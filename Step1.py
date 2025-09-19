import csv
import os
import re
import time
import argparse
from typing import List, Tuple, Optional, Dict, Iterator
import requests

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"

def slugify(name: str) -> str:
    name = name.strip().replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    return name or "unknown"

def read_species_csv(path: str) -> List[Tuple[str, str]]:
    """
    Read the species name and taxon id from the CSV.
    Automatically search for the column named 'taxon' as the taxid
    By default, the first column is used as the species name (if there is a 'species' column, it is preferred).
    """
    items: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if "Taxon" not in reader.fieldnames:
            raise SystemExit(f"CSV is missing the 'Taxon' column. The column names are: {reader.fieldnames}")

        species_col = "species" if "species" in reader.fieldnames else reader.fieldnames[0]

        for row in reader:
            tax = row["Taxon"].strip()
            name = row[species_col].strip()
            if not tax.isdigit():
                continue
            items.append((name, tax))
    return items

def quote_phrase(s: str) -> str:
    return '"' + s.replace('"', r'\"') + '"'

def build_query_taxonomy(taxid: str, searchword: str,
                         reviewed: bool = True, no_fragment: bool = False) -> str:
    protein_term = f'protein_name:{quote_phrase(searchword)}'
    parts = [f"({protein_term})", f"(taxonomy_id:{taxid})"]
    if reviewed:
        parts.append("(reviewed:true)")
    if no_fragment:
        parts.append("(fragment:false)")
    return " AND ".join(parts)

def fetch_fasta(query: str, timeout: int = 60) -> Optional[bytes]:
    params = {"format": "fasta", "query": query}
    req = requests.Request("GET", UNIPROT_STREAM_URL, params=params).prepare()
    print(f"  URL: {req.url}")
    with requests.Session() as s:
        with s.send(req, stream=True, timeout=timeout) as r:
            if r.status_code == 204:
                return b""
            r.raise_for_status()
            return b"".join(r.iter_content(chunk_size=8192))

def iter_fasta(text: str) -> Iterator[Tuple[str, str]]:
    header, buf = None, []
    for line in text.splitlines():
        line = line.rstrip("\r\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(buf)
            header = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if header is not None:
        yield header, "".join(buf)

def wrap_seq(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def parse_uniprot_header(h: str) -> Dict[str, Optional[str]]:
    m = re.match(r"^(?P<db>sp|tr)\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<rest>.*)$", h)
    if not m:
        return {"db": None, "accession": None, "entry": None,
                "description": h, "OS": None, "OX": None, "GN": None}
    db, acc, entry, rest = m.group("db"), m.group("acc"), m.group("entry"), m.group("rest")

    fields = {}
    for tag in (" OS=", " OX=", " GN=", " PE=", " SV="):
        idx = rest.find(tag)
        if idx >= 0:
            fields[tag.strip()] = idx
    tag_positions = sorted(fields.items(), key=lambda x: x[1])

    desc_end = tag_positions[0][1] if tag_positions else len(rest)
    description = rest[:desc_end].strip()

    tag_values: Dict[str, Optional[str]] = {"OS": None, "OX": None, "GN": None}
    for i, (tag, start) in enumerate(tag_positions):
        val_start = start + len(tag)
        val_end = len(rest) if i == len(tag_positions) - 1 else tag_positions[i+1][1]
        key = tag.strip().replace("=", "")
        if key in tag_values:
            tag_values[key] = rest[val_start:val_end].strip()

    return {"db": db, "accession": acc, "entry": entry,
            "description": description if description else None,
            "OS": tag_values["OS"], "OX": tag_values["OX"], "GN": tag_values["GN"]}


def main():
    ap = argparse.ArgumentParser(description="Download the UniProt protein and combine core.fasta")
    ap.add_argument("--data_dir", help="The data directory should include whitelist.csv and core.fasta")
    ap.add_argument("--searchword", required=True, help="Keyword for searching, eg. \"Chitin synthase\"")
    ap.add_argument("--no-reviewed", dest="reviewed", action="store_false")
    ap.set_defaults(reviewed=True)
    ap.add_argument("--no-fragment", action="store_true")
    ap.add_argument("--sleep", type=float, default=0.5)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--dedup", choices=["accession", "none"], default="accession")
    args = ap.parse_args()

    data_dir = args.data_dir
    whitelist_csv = os.path.join(data_dir, "whitelist.csv")
    core_fasta = os.path.join(data_dir, "core.fasta")
    merge_fasta = os.path.join(data_dir, "merge.fasta")
    meta_csv = os.path.join(data_dir, "meta.csv")
    all_fasta = os.path.join(data_dir, "all.fasta")

    species_list = read_species_csv(whitelist_csv)
    if not species_list:
        raise SystemExit("The whitelist.csv format is incorrect or empty.")

    seen_accessions = set()
    fasta_fh = open(merge_fasta, "w", encoding="utf-8")
    meta_fh = open(meta_csv, "w", encoding="utf-8-sig", newline="")
    meta_writer = csv.writer(meta_fh)
    meta_writer.writerow([
        "input_species_name","input_taxid","searchword",
        "db","accession","entry","gene","description",
        "organism_OS","organism_taxid_OX","length_aa"
    ])

    try:
        for i,(sp_name,taxid) in enumerate(species_list,1):
            query = build_query_taxonomy(taxid, args.searchword,
                                         reviewed=args.reviewed,
                                         no_fragment=args.no_fragment)
            print(f"[{i}/{len(species_list)}] {sp_name} (taxid={taxid}) ...")
            for attempt in range(1, args.retries+1):
                try:
                    fasta_bytes = fetch_fasta(query)
                    if fasta_bytes is None: raise RuntimeError("Request failed")
                    if not fasta_bytes:
                        print("  -> No result")
                        break
                    text = fasta_bytes.decode("utf-8","replace")
                    for header,seq in iter_fasta(text):
                        meta = parse_uniprot_header(header)
                        acc = meta.get("accession")
                        if args.dedup=="accession" and acc and acc in seen_accessions:
                            continue
                        if acc: seen_accessions.add(acc)
                        fasta_fh.write(f">{header}\n{wrap_seq(seq)}\n")
                        meta_writer.writerow([
                            sp_name,taxid,args.searchword,
                            meta["db"] or "", acc or "", meta["entry"] or "", meta["GN"] or "",
                            meta["description"] or "", meta["OS"] or "", meta["OX"] or "",
                            len(seq)
                        ])
                    break
                except Exception as e:
                    print(f"  ERROR：{e}")
                    time.sleep(args.sleep*attempt)
            time.sleep(args.sleep)
    finally:
        fasta_fh.close()
        meta_fh.close()

    # 合并 core.fasta
    with open(all_fasta,"w",encoding="utf-8") as out_f:
        with open(merge_fasta,"r",encoding="utf-8") as f_in:
            out_f.writelines(f_in)
        if os.path.exists(core_fasta):
            with open(core_fasta,"r",encoding="utf-8") as f_in:
                out_f.writelines(f_in)
            print(f"[INFO] core.fasta has been added")
        else:
            print(f"[WARN] core.fasta was not found")

    os.remove(merge_fasta)
    print(f"[DONE] all.fasta: {all_fasta}")

if __name__ == "__main__":
    main()
