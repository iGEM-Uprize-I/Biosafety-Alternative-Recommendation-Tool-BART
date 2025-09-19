import re, csv, subprocess, argparse, sys
from pathlib import Path


parser = argparse.ArgumentParser(description="StepA: core vs all (mmseqs2) + mafft + motif QC")
parser.add_argument("--data", required=True, help="数据目录 (包含 all.fasta, core.fasta, motif.txt)")
parser.add_argument("--outdir", required=True, help="结果输出目录")
args = parser.parse_args()

DATA_DIR = Path(args.data)
OUTDIR   = Path(args.outdir)

FASTA      = DATA_DIR/"all.fasta"
CORE_FASTA = DATA_DIR/"core.fasta"
MOTIF_FILE = DATA_DIR/"motif.txt"

MMSEQS_DIR = OUTDIR/"mmseqs_work"
MMSEQS_M8  = MMSEQS_DIR/"core_vs_all.m8"
MSA_FA     = OUTDIR/"msa.fa"
OUT_CSV    = OUTDIR/"core_vs_all.csv"

OUTDIR.mkdir(parents=True, exist_ok=True)
MMSEQS_DIR.mkdir(parents=True, exist_ok=True)


def read_fasta_to_dict(path: Path):
    seqs = {}
    with open(path) as f:
        sid, buf = None, []
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if sid: seqs[sid] = "".join(buf)
                sid = s[1:].split()[0]  
                buf = []
            else:
                buf.append(s)
        if sid: seqs[sid] = "".join(buf)
    return seqs

def read_fasta_ids(path: Path):
    return set(read_fasta_to_dict(path).keys())

def load_motifs(path: Path):
    motifs, names = [], []
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            motifs.append(re.compile(line))
            names.append(line)
    if not motifs:
        sys.exit("[ERR] motif.txt is empty. Please write at least one line of regular expression")
    return motifs, names


def motifs_hits(seq_aln: str, motifs):
    raw = seq_aln.replace("-", "")
    return [bool(m.search(raw)) for m in motifs]

def parse_mmseqs(m8_path: Path, fasta_ids=None):
    """
    Parse the output of mmseqs2 convertalis.
    Return {full_target_id: (pident, qcov, alnlen, evalue, bits)}
    """
    hits = {}
    if not m8_path.exists():
        return hits

    # accession → full fasta id 映射
    acc2full = {}
    if fasta_ids:
        for fid in fasta_ids:
            if "|" in fid:
                parts = fid.split("|")
                for p in parts:
                    if p and p.isalnum() and len(p) >= 5:
                        acc2full[p] = fid
            else:
                acc2full[fid] = fid

    with open(m8_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            _, t = parts[0], parts[1]  
            pid, qcov, alnlen, evalue, bits = parts[2:7]
            full_t = acc2full.get(t, t)
            hits[full_t] = (
                float(pid),
                float(qcov),
                int(alnlen),
                float(evalue.replace("E", "e")),
                float(bits),
            )
    return hits


# ============ 0. Read the core ID & motif ============
if not FASTA.exists() or not CORE_FASTA.exists() or not MOTIF_FILE.exists():
    sys.exit("[ERR] Data directory must contain all.fasta, core.fasta, motif.txt")

CORE_IDS = read_fasta_ids(CORE_FASTA)
print(f"[INFO] Number of core sequences: {len(CORE_IDS)}")

MOTIFS, MOTIF_NAMES = load_motifs(MOTIF_FILE)
print(f"[INFO] Number of motif entries: {len(MOTIFS)} ({', '.join(MOTIF_NAMES)})")

# ============ 1) mmseqs2 core vs all ============
if not MMSEQS_M8.exists():
    print("[RUN] mmseqs2 core vs all...")
    coreDB = MMSEQS_DIR/"coreDB"
    allDB  = MMSEQS_DIR/"allDB"
    RES    = MMSEQS_DIR/"res"
    TMP    = MMSEQS_DIR/"tmp"

    subprocess.run(["mmseqs","createdb",str(CORE_FASTA),str(coreDB)], check=True)
    subprocess.run(["mmseqs","createdb",str(FASTA),str(allDB)], check=True)

    subprocess.run(["mmseqs","search",str(coreDB),str(allDB),str(RES),str(TMP),
                    "--cov-mode","1","-c","0.7","--min-seq-id","0.25"], check=True)

    subprocess.run([
        "mmseqs","convertalis",str(coreDB),str(allDB),str(RES),str(MMSEQS_M8),
        "--format-output","query,target,pident,qcov,alnlen,evalue,bits"
    ], check=True)
    print("[OK] mmseqs2 complete")
else:
    print("[SKIP] mmseqs2 existed")

# ============ 2) mafft ============
if not MSA_FA.exists():
    print("[RUN] mafft...")
    with open(OUTDIR/"mafft.log","w") as log, open(MSA_FA,"w") as out:
        subprocess.run(["mafft","--auto",str(FASTA)], stdout=out, stderr=log, check=True)
    print("[OK] mafft complete")
else:
    print("[SKIP] mafft existed")

# ============ 3) parse & merge ============
msa  = read_fasta_to_dict(MSA_FA)
hits = parse_mmseqs(MMSEQS_M8)
print(hits)

missing_cores = [cid for cid in CORE_IDS if cid not in msa]
if missing_cores:
    print(f"[WARN] The following core ids were not found in MSA (check ID consistency) : {', '.join(missing_cores)}")

rows=[]
for pid, aln in msa.items():
    species = pid.split("|")[1] if "|" in pid else "."
    is_core = 1 if pid in CORE_IDS else 0

    hit_info = hits.get(species)
    if hit_info:
        pidval, covval, alnlen, evalue, bits = hit_info
        hit_ok = 1
    else:
        pidval, covval, alnlen, evalue, bits = ".", ".", ".", ".", "."
        hit_ok = 0

    motif_hit_vec = motifs_hits(aln, MOTIFS)
    motifs_all = 1 if all(motif_hit_vec) else 0

    fail = ""
    if hit_ok==0 or motifs_all==0:
        fail = "Fail-A1"

    rows.append([species, pid, is_core, hit_ok, pidval, covval, evalue, bits] 
                + motif_hit_vec + [motifs_all, fail])

# ============ 4. Output CSV ============
with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    header = ["accession","protein_id","is_core","hit_core(0/1)",
              "pid(%)","aln_cov(%)","evalue","bitscore"] \
             + MOTIF_NAMES + ["motifs_ok(0/1)","fail_tag"]
    w.writerow(header)
    w.writerows(rows)

print(f"[OK] Write the {OUT_CSV} ({len(rows)} entries)")
