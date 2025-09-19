import os
import re
import csv
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import requests
import numpy as np
from scipy.spatial import cKDTree
from Bio import SeqIO
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
from tqdm import tqdm

warnings.simplefilter("ignore", PDBConstructionWarning)

ALPHAFOLD_PATTERNS = [
    "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v3.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v2.pdb",
]

AA_CAT = {
    "A":"hyd","V":"hyd","I":"hyd","L":"hyd","M":"hyd",   # Hydrophobic
    "F":"aro","Y":"aro","W":"aro",                      # Aromatic
    "S":"pol","T":"pol","N":"pol","Q":"pol","C":"pol",  # Polar
    "K":"pos","R":"pos","H":"pos",                      # Positive charge
    "D":"neg","E":"neg",                                # Negative charge
    "G":"spc","P":"spc",                                # Special
}

def parse_args():
    ap = argparse.ArgumentParser(description="Pocket mapping + pLDDT QC + pocket 3D scoring")
    ap.add_argument("--core-fasta", required=True, help="Core/reference sequence FASTA file (only one sequence)")
    ap.add_argument("--msa", default="result/msa.fa", help="MAFFT result FASTA (e.g., result/msa.fa)")
    ap.add_argument("--pocket", required=True, help="Core pocket residue numbers (1-based), one per line")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--radius", type=float, default=3.0, help="ShapeOverlap nearest-neighbor radius (Å), default 3.0")
    ap.add_argument("--motif-file", default="data/motif.txt", help="File whose FIRST line is used as the anchor motif (default: data/motif.txt)")
    return ap.parse_args()

def guess_accession_from_header(header: str) -> Optional[str]:
    # Typical: sp|Q9XXXX|NAME or tr|A0A0...| or Q9XXXX
    m = re.search(r"\|(A0A[A-Z0-9]+|[OPQ][0-9][A-Z0-9]{3}[0-9])\|", header)
    if m:
        return m.group(1)
    m2 = re.search(r"\b([A-NR-Z0-9]{6,10})\b", header)
    if m2:
        return m2.group(1)
    return None

def read_pocket_positions(txt_path: str) -> List[int]:
    pos = []
    with open(txt_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if not s.isdigit():
                raise ValueError(f"Non-numeric line in pocket file: {s}")
            pos.append(int(s))
    return sorted(set(pos))

def read_first_motif(motif_path: str) -> str:
    """Read the FIRST non-empty line from motif file and return it as the motif string."""
    with open(motif_path, "r", encoding="utf-8") as f:
        for line in f:
            mot = line.strip()
            if mot:
                return mot
    raise ValueError(f"No non-empty motif line found in: {motif_path}")

def parse_msa(msa_path: str) -> List[SeqIO.SeqRecord]:
    return list(SeqIO.parse(msa_path, "fasta"))

def build_pos_maps(aln_seq: str) -> Tuple[Dict[int,int], Dict[int,int]]:
    """Build mapping between ungapped sequence positions (1-based) and alignment indices (0-based)."""
    map_seq2aln: Dict[int,int] = {}
    map_aln2seq: Dict[int,int] = {}
    seq_pos = 0
    for aln_i, ch in enumerate(aln_seq):
        if ch != "-":
            seq_pos += 1
            map_seq2aln[seq_pos] = aln_i
            map_aln2seq[aln_i] = seq_pos
        else:
            map_aln2seq[aln_i] = 0
    return map_seq2aln, map_aln2seq

def map_core_pocket_to_target(core_aln: str, tgt_aln: str, core_pocket_1b: List[int]) -> List[int]:
    """Map core pocket positions (1-based on ungapped core) to target ungapped positions via alignment columns."""
    core_s2a, _ = build_pos_maps(core_aln)
    _, tgt_a2s   = build_pos_maps(tgt_aln)
    out = []
    for p in core_pocket_1b:
        aln_idx = core_s2a.get(p)
        if aln_idx is None:
            continue
        tgt_pos = tgt_a2s.get(aln_idx, 0)
        if tgt_pos > 0:
            out.append(tgt_pos)
    return sorted(set(out))

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def download_alphafold_pdb(acc: str, out_path: Path) -> bool:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for url in ALPHAFOLD_PATTERNS:
        u = url.format(acc=acc)
        try:
            r = requests.get(u, timeout=30)
            if r.status_code == 200 and len(r.content) > 1000:
                out_path.write_bytes(r.content)
                return True
        except requests.RequestException:
            pass
    return False

def load_structure_coords_and_plddt(pdb_path: Path, chain_prefer: Optional[str]=None):
    """Return (seq_1letter, CA_xyz[N,3], pLDDT_per_res[N], seqIndex(1-based)->listIndex)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af", str(pdb_path))
    model = next(structure.get_models())

    chain = None
    if chain_prefer:
        for c in model:
            if c.id == chain_prefer:
                chain = c; break
    if chain is None:
        chain = next(model.get_chains())

    seq = []
    ca_xyz = []
    res_plddt = []
    seqidx_to_listi = {}

    i_seq = 0
    for res in chain:
        hetfield, resseq, icode = res.id
        if hetfield != " ":
            continue
        resname = res.get_resname()
        try:
            aa = three_to_one(resname)
        except Exception:
            continue
        if "CA" not in res:
            continue

        i_seq += 1
        seq.append(aa)
        ca = res["CA"].get_coord()
        ca_xyz.append(ca)

        bf = [atom.get_bfactor() for atom in res.get_atoms()]
        plddt = float(np.mean(bf)) if bf else np.nan
        res_plddt.append(plddt)

        seqidx_to_listi[i_seq] = len(ca_xyz)-1

    return "".join(seq), np.array(ca_xyz, dtype=float), np.array(res_plddt, dtype=float), seqidx_to_listi

def three_to_one(resname: str) -> str:
    resname = resname.upper().strip()
    table = {
        "ALA":"A","VAL":"V","ILE":"I","LEU":"L","MET":"M","PHE":"F","TYR":"Y","TRP":"W",
        "SER":"S","THR":"T","ASN":"N","GLN":"Q","CYS":"C","LYS":"K","ARG":"R","HIS":"H",
        "ASP":"D","GLU":"E","GLY":"G","PRO":"P",
    }
    if resname in table:
        return table[resname]
    raise KeyError(resname)

# ====== Scoring constants (tunable) ======
ANCHOR_WEIGHT  = 2.0       # multiplier for anchor positions
PLDDT_GAMMA    = 1.5       # exponent for pLDDT-based weight; >1 emphasizes high-confidence sites

def find_core_anchor_positions(core_seq: str, motif: Optional[str] = None) -> set:
    """Find 1-based indices of the anchor motif in the GAP-REMOVED core sequence."""
    m = (motif or ANCHOR_MOTIF).upper()
    s = core_seq.upper()
    pos = set(); i = 0
    while True:
        j = s.find(m, i)
        if j == -1:
            break
        for k in range(j+1, j+1+len(m)):
            pos.add(k)
        i = j + 1
    return pos

# ====== Side-chain chemistry similarity ======
_CONSERVATIVE = {
    ("I","L"),("I","V"),("L","V"),
    ("F","Y"),("F","W"),("Y","W"),
    ("K","R"),("K","H"),("R","H"),
    ("S","T"),
    ("N","Q"),
    ("D","E"),("D","N"),("E","Q"),
}
def _cons(a: str, b: str) -> bool:
    a=a.upper(); b=b.upper()
    return (a,b) in _CONSERVATIVE or (b,a) in _CONSERVATIVE

def chem_similarity(a: str, b: str) -> float:
    if a == b and a in AA_CAT: return 1.0
    ca = AA_CAT.get(a); cb = AA_CAT.get(b)
    if (ca is not None) and (cb is not None) and (ca == cb): return 0.7
    if _cons(a,b): return 0.5
    return 0.0

def compute_pair_weights(core_idxmap, tgt_idxmap, core_plddt, tgt_plddt,
                         core_pos_valid: List[int], tgt_pos_valid: List[int],
                         anchor_pos_set: set) -> Tuple[np.ndarray, float]:
    """Weight per paired site: (min pLDDT)^gamma × (anchor multiplier). Also return mean weight as 'quality'."""
    w = []
    for cp, tp in zip(core_pos_valid, tgt_pos_valid):
        ic = core_idxmap.get(cp); it = tgt_idxmap.get(tp)
        if ic is None or it is None:
            w.append(0.0); continue
        qc = float(core_plddt[ic]) if np.isfinite(core_plddt[ic]) else 0.0
        qt = float(tgt_plddt[it])  if np.isfinite(tgt_plddt[it])  else 0.0
        base = max(0.0, min(qc, qt)) / 100.0
        base = base ** PLDDT_GAMMA
        mult = ANCHOR_WEIGHT if (cp in anchor_pos_set) else 1.0
        w.append(base * mult)
    w = np.array(w, dtype=float)
    quality = float(np.mean(w)) if np.any(w > 0) else 0.0
    return w, quality

def weighted_kabsch(P: np.ndarray, Q: np.ndarray, w: np.ndarray) -> Tuple[np.ndarray, float]:
    """Weighted rigid alignment (Kabsch) and weighted RMSD."""
    idx = np.where(w > 0)[0]
    if idx.size < 3:
        raise ValueError("Insufficient weighted points for rigid fitting")
    P2, Q2, w2 = P[idx], Q[idx], w[idx]
    W = w2 / np.sum(w2)
    muP = np.sum(P2 * W[:,None], axis=0)
    muQ = np.sum(Q2 * W[:,None], axis=0)
    Pc, Qc = P2 - muP, Q2 - muQ
    C = Qc.T @ (W[:,None] * Pc)
    V,S,Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1,1,d]) @ Wt
    Q_fit_full = (Q - muQ) @ R + muP
    diff = P2 - ((Q2 - muQ) @ R + muP)
    wrmsd = float(np.sqrt(np.sum(w2 * np.sum(diff**2, axis=1)) / np.sum(w2)))
    return Q_fit_full, wrmsd

def weighted_shape_overlap(P: np.ndarray, Q: np.ndarray,
                           wP: np.ndarray, wQ: np.ndarray, radius: float=3.0) -> float:
    """Weighted bidirectional coverage within a radius; average of (P covers Q) and (Q covers P)."""
    if len(P)==0 or len(Q)==0: return 0.0
    tP = cKDTree(P); tQ = cKDTree(Q)
    dQP,_ = tP.query(Q, k=1)  # Q covered by P
    dPQ,_ = tQ.query(P, k=1)  # P covered by Q
    covQ = float(np.sum((dQP <= radius) * wQ) / (np.sum(wQ) + 1e-12))
    covP = float(np.sum((dPQ <= radius) * wP) / (np.sum(wP) + 1e-12))
    return (covP + covQ) / 2.0

def weighted_chem_conservation(core_seq3d: str, tgt_seq3d: str,
                               core_pos_valid: List[int], tgt_pos_valid: List[int],
                               w: np.ndarray) -> float:
    """Weighted side-chain chemistry conservation over paired pocket positions."""
    num = 0.0; den = 0.0
    for i,(cp,tp) in enumerate(zip(core_pos_valid, tgt_pos_valid)):
        wi = float(w[i])
        if wi <= 0: continue
        a = core_seq3d[cp-1] if cp-1 < len(core_seq3d) else "X"
        b = tgt_seq3d[tp-1]  if tp-1  < len(tgt_seq3d)  else "X"
        num += chem_similarity(a,b) * wi
        den += wi
    return float(num / (den + 1e-12))

def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, float]:
    """Unweighted Kabsch for reference; returns aligned Q and RMSD."""
    assert P.shape == Q.shape and P.shape[1] == 3
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Qc.T @ Pc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1,1,d])
    R = V @ D @ Wt
    Q_fit = Qc @ R + Pc.mean(axis=0)
    rmsd = float(np.sqrt(np.mean(np.sum((P - Q_fit)**2, axis=1))))
    return Q_fit, rmsd

def shape_overlap(P: np.ndarray, Q: np.ndarray, radius: float=2.0) -> float:
    """Unweighted shape overlap (for reference)."""
    if len(P)==0 or len(Q)==0:
        return 0.0
    tP = cKDTree(P); tQ = cKDTree(Q)
    covP = np.mean(tP.query(Q, k=1)[0] <= radius)
    covQ = np.mean(tQ.query(P, k=1)[0] <= radius)
    return float((covP + covQ)/2.0)

def sidechain_category_match(core_seq: str, tgt_seq: str, core_pos_1b: List[int], tgt_pos_1b: List[int]) -> float:
    """Unweighted categorical match of side-chain types over paired positions (for reference)."""
    L = min(len(core_pos_1b), len(tgt_pos_1b))
    if L == 0:
        return 0.0
    agree = 0
    for i in range(L):
        a = core_seq[ core_pos_1b[i]-1 ] if core_pos_1b[i]-1 < len(core_seq) else "X"
        b = tgt_seq[  tgt_pos_1b[i]-1 ] if tgt_pos_1b[i]-1 < len(tgt_seq)  else "X"
        ca = AA_CAT.get(a, "unk")
        cb = AA_CAT.get(b, "unk")
        if ca == cb and ca != "unk":
            agree += 1
    return float(agree / L)

# ---------- Read core header and locate the corresponding MSA record ----------
def read_core_header(core_fasta: str) -> Tuple[str, str]:
    recs = list(SeqIO.parse(core_fasta, "fasta"))
    if len(recs) == 0:
        raise ValueError(f"core-fasta is empty: {core_fasta}")
    if len(recs) > 1:
        print(f"[WARN] Multiple sequences in core-fasta ({len(recs)}); using the first: {recs[0].id}")
    r = recs[0]
    return r.id, (r.description or r.id)

def find_core_record_auto(records: List[SeqIO.SeqRecord], core_header: str, core_id: str) -> SeqIO.SeqRecord:
    # 1) Prefer exact id match
    for r in records:
        if r.id == core_id:
            return r
    # 2) Match first token of description (substring)
    head_token = core_header.split()[0]
    cand = [r for r in records if (head_token in r.id) or (head_token in (r.description or ""))]
    if len(cand) == 1:
        return cand[0]
    if len(cand) > 1:
        print(f"[WARN] Multiple hits by header token; using the first: {cand[0].id}")
        return cand[0]
    # 3) Try UniProt accession from header
    acc = guess_accession_from_header(core_header)
    if acc:
        cand2 = [r for r in records if (acc in r.id) or (acc in (r.description or ""))]
        if len(cand2) >= 1:
            if len(cand2) > 1:
                print(f"[WARN] Multiple hits by accession; using the first: {cand2[0].id}")
            return cand2[0]
    # 4) Fallback: substring of core_id
    cand3 = [r for r in records if (core_id in r.id) or (core_id in (r.description or ""))]
    if len(cand3) >= 1:
        if len(cand3) > 1:
            print(f"[WARN] Multiple hits by core_id substring; using the first: {cand3[0].id}")
        return cand3[0]
    raise ValueError("Cannot locate the core sequence in MSA (id/description/accession all failed).")

# ------------- Main --------------
def main():
    args = parse_args()
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    struct_dir = outdir / "structures"
    ensure_dir(struct_dir)

    # Read anchor motif from file FIRST LINE and override the global ANCHOR_MOTIF
    global ANCHOR_MOTIF
    try:
        ANCHOR_MOTIF = read_first_motif(args.motif_file)
        print(f"[INFO] Anchor motif loaded from {args.motif_file}: {ANCHOR_MOTIF}")
    except Exception as e:
        print(f"[WARN] Failed to read motif from {args.motif_file}, fallback to default '{ANCHOR_MOTIF}': {e}")

    # 1) Read MSA and locate the core record
    core_id_token, core_header = read_core_header(args.core_fasta)
    records = parse_msa(args.msa)
    core_rec = find_core_record_auto(records, core_header, core_id_token)
    core_aln = str(core_rec.seq).upper()
    core_seq = core_aln.replace("-", "")

    # Build MSA dicts
    id2rec = {r.id: r for r in records}
    all_ids = [r.id for r in records]

    # 2) Read core pocket positions (1-based on gap-removed core sequence)
    core_pocket_pos = read_pocket_positions(args.pocket)

    # 3) Accession map (if needed for overrides; currently empty)
    acc_map = {}

    # 4) Map core pocket to each sequence via alignment columns
    pocket_map: Dict[str, List[int]] = {}
    pocket_map["core"] = core_pocket_pos

    print("[INFO] Mapping core pocket positions to candidate sequences...")
    for pid in all_ids:
        aln = str(id2rec[pid].seq).upper()
        if pid == core_rec.id:
            pocket_map[pid] = core_pocket_pos
        else:
            tgt_pos = map_core_pocket_to_target(core_aln, aln, core_pocket_pos)
            pocket_map[pid] = tgt_pos

    with open(outdir/"pocket_residues.json", "w", encoding="utf-8") as f:
        json.dump(pocket_map, f, indent=2, ensure_ascii=False)

    # 5) Download/read AlphaFold structures; extract seq, CA coords, pLDDT
    qc_rows = []
    score_rows = []

    def get_acc(pid: str) -> Optional[str]:
        if pid in acc_map and acc_map[pid]:
            return acc_map[pid]
        header = id2rec[pid].description or id2rec[pid].id
        return guess_accession_from_header(header)

    # Load core structure
    core_acc = get_acc(core_rec.id)
    if not core_acc:
        raise RuntimeError(f"Cannot infer accession for core={core_rec.id}; please provide an accession mapping")
    core_pdb = struct_dir / f"{core_acc}.pdb"
    if not core_pdb.exists():
        ok = download_alphafold_pdb(core_acc, core_pdb)
        if not ok:
            raise RuntimeError(f"AlphaFold download failed for core {core_acc}. Provide PDB manually at: {core_pdb}")

    core_seq3d, core_ca, core_plddt, core_idxmap = load_structure_coords_and_plddt(core_pdb)

    # Restrict core pocket to CA-present residues in structure
    core_pocket_3d_idx = []
    for p in core_pocket_pos:
        li = core_idxmap.get(p)
        if li is not None:
            core_pocket_3d_idx.append(li)
    core_pocket_ca = core_ca[core_pocket_3d_idx, :]

    print("[INFO] Structure download + pLDDT stats + pocket scoring ...")
    for pid in tqdm(all_ids):
        acc = get_acc(pid)
        if not acc:
            print(f"[WARN] {pid} has no inferred accession; skip")
            continue
        pdb_path = struct_dir / f"{acc}.pdb"
        if not pdb_path.exists():
            ok = download_alphafold_pdb(acc, pdb_path)
            if not ok:
                print(f"[WARN] {pid} ({acc}) AlphaFold download failed; skip")
                continue

        try:
            seq3d, ca_xyz, plddt_arr, idxmap = load_structure_coords_and_plddt(pdb_path)
        except Exception as e:
            print(f"[WARN] Failed to read structure for {pid} ({acc}): {e}")
            continue

        tgt_pocket_pos = pocket_map.get(pid, [])
        pocket_list_idx = [idxmap[p] for p in tgt_pocket_pos if p in idxmap]
        if len(pocket_list_idx) == 0:
            print(f"[WARN] {pid} pocket mapping empty or missing Cα; skip scoring")
            continue

        pocket_plddt = plddt_arr[pocket_list_idx]
        pocket_mean = float(np.nanmean(pocket_plddt)) if pocket_plddt.size else float("nan")
        low70_pct = float(np.mean(pocket_plddt < 70.0)) if pocket_plddt.size else 1.0
        qc_pass = (low70_pct <= 0.30)

        qc_rows.append({
            "protein_id": pid,
            "pocket_pLDDT_mean": f"{pocket_mean:.2f}",
            "pocket_low70%": f"{100.0*low70_pct:.2f}",
            "pass": "pass" if qc_pass else "fail",
        })

        if not qc_pass:
            continue

        try:
            # Pair strictly by the SAME alignment columns to avoid mis-zip
            core_s2a, core_a2s = build_pos_maps(core_aln)               # core alignment (with '-')
            tgt_s2a,  tgt_a2s  = build_pos_maps(str(id2rec[pid].seq).upper())

            core_pos_set = set(pocket_map["core"])      # 1-based on GAP-REMOVED core
            tgt_pos_set  = set(pocket_map[pid])         # 1-based on GAP-REMOVED target

            paired = []
            core_pos_valid = []
            tgt_pos_valid  = []

            # Iterate over alignment columns; pair when both core and target positions exist and are in pocket sets
            for col in range(len(core_aln)):
                cpos = core_a2s.get(col, 0)   # 0 if gap
                tpos = tgt_a2s.get(col, 0)
                if cpos in core_pos_set and tpos in tgt_pos_set:
                    li_core = core_idxmap.get(cpos)
                    li_tgt  = idxmap.get(tpos)
                    if (li_core is not None) and (li_tgt is not None):
                        core_pos_valid.append(cpos)
                        tgt_pos_valid.append(tpos)
                        paired.append((li_core, li_tgt))

            if len(paired) < 3:
                print(f"[WARN] {pid} has too few common pocket positions ({len(paired)}); skip scoring")
                continue

            P = core_ca[[i for i,_ in paired], :]
            Q = ca_xyz[[j for _,j in paired], :]

            # --- Weighted scoring ---
            # 1) Per-pair weight: (pLDDT^gamma) × (anchor multiplier only for anchor-in-pocket)
            anchor_core_pos  = find_core_anchor_positions(core_seq, motif=ANCHOR_MOTIF)  # core GAP-REMOVED
            anchor_in_pocket = anchor_core_pos & set(core_pos_valid)
            w_pair, quality  = compute_pair_weights(
                core_idxmap, idxmap, core_plddt, plddt_arr,
                core_pos_valid, tgt_pos_valid,
                anchor_in_pocket
            )
            if not np.any(w_pair > 0):
                w_pair = np.ones(len(paired), dtype=float)
                quality = 0.0

            # 2) Weighted Kabsch + weighted RMSD
            Q_fit, wrmsd = weighted_kabsch(P, Q, w_pair)

            # 3) Weighted shape overlap (each side weighted by w_pair)
            sh_ov = weighted_shape_overlap(P, Q_fit, w_pair, w_pair, radius=args.radius)

            # 4) Weighted chemical conservation
            sc_match = weighted_chem_conservation(core_seq3d, seq3d, core_pos_valid, tgt_pos_valid, w_pair)

            # 5) Map wRMSD to [0,1] robustly
            norm_wrmsd = 1.0 / (1.0 + (wrmsd/3.0)**2)

            # 6) Final composite score (geometry > RMSD > chemistry > quality)
            pocket3d = 0.45*sh_ov + 0.30*norm_wrmsd + 0.20*sc_match + 0.05*quality

            score_rows.append({
                "protein_id": pid,
                "RMSD": f"{wrmsd:.3f}",
                "ShapeOverlap": f"{sh_ov:.3f}",
                "SidechainChem": f"{sc_match:.3f}",
                "Pocket3D_Score": f"{pocket3d:.3f}",
            })
        except Exception as e:
            print(f"[WARN] Scoring failed for {pid}: {e}")
            continue

    qc_csv = outdir/"qc_plddt.csv"
    with open(qc_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["protein_id","pocket_pLDDT_mean","pocket_low70%","pass"])
        w.writeheader()
        for row in qc_rows:
            w.writerow(row)

    score_csv = outdir/"pocket_scores.csv"
    with open(score_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["protein_id","RMSD","ShapeOverlap","SidechainChem","Pocket3D_Score"])
        w.writeheader()
        for row in score_rows:
            w.writerow(row)

    print(f"[DONE] pocket_residues.json  → {outdir/'pocket_residues.json'}")
    print(f"[DONE] qc_plddt.csv         → {qc_csv}")
    print(f"[DONE] pocket_scores.csv    → {score_csv}")

    # ============ Final aggregation: result/final_result.csv ============
    # Read three CSVs from fixed paths and write the intersection to result/final_result.csv
    res_dir = Path("result")
    res_dir.mkdir(parents=True, exist_ok=True)

    path_core   = res_dir / "core_vs_all.csv"
    path_qc     = res_dir / "qc_plddt.csv"
    path_score  = res_dir / "pocket_scores.csv"
    out_final   = res_dir / "final_result.csv"

    # core_vs_all: keep rows where fail_tag is empty
    core_pass = {}
    core_fields = []
    with open(path_core, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        core_fields = r.fieldnames or []
        for row in r:
            fail = (row.get("fail_tag") or "").strip()
            if fail == "":
                pid = (row.get("protein_id") or "").strip()
                if pid:
                    core_pass[pid] = row

    # qc_plddt: keep rows where pass == "pass"
    qc_pass = {}
    qc_fields = []
    with open(path_qc, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        qc_fields = r.fieldnames or []
        for row in r:
            if (row.get("pass") or "").strip().lower() == "pass":
                pid = (row.get("protein_id") or "").strip()
                if pid:
                    qc_pass[pid] = row

    # pocket_scores: keep all rows (map by protein_id)
    score_rows_map = {}
    score_fields = []
    with open(path_score, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        score_fields = r.fieldnames or []
        for row in r:
            pid = (row.get("protein_id") or "").strip()
            if pid:
                score_rows_map[pid] = row

    # Intersection: protein_id present in core_pass AND qc_pass AND score_rows_map
    keep_ids = [pid for pid in core_pass.keys() if pid in qc_pass and pid in score_rows_map]

    # Compose output header: all core_vs_all fields (+ protein_id if missing) + selected qc + selected score fields
    cols_core  = core_fields[:] if "protein_id" in core_fields else ["protein_id"] + core_fields
    cols_qc    = [c for c in ["pocket_pLDDT_mean","pocket_low70%","pass"] if c in qc_fields]
    cols_score = [c for c in ["RMSD","ShapeOverlap","SidechainChem","Pocket3D_Score"] if c in score_fields]

    final_header = []
    seen = set()
    for c in cols_core + [c for c in cols_qc + cols_score if c not in cols_core]:
        if c not in seen:
            final_header.append(c); seen.add(c)

    with open(out_final, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=final_header)
        w.writeheader()
        for pid in keep_ids:
            row = {}
            row.update(core_pass[pid])         # base info from core_vs_all
            row["protein_id"] = pid            # ensure protein_id exists
            for c in cols_qc:
                row[c] = qc_pass[pid].get(c, "")
            for c in cols_score:
                row[c] = score_rows_map[pid].get(c, "")
            w.writerow(row)

    print(f"[DONE] final_result.csv    → {out_final}")

if __name__ == "__main__":
    main()
