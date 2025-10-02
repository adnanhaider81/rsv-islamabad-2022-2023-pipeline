#!/usr/bin/env python3
import os, sys, argparse
from Bio import Entrez, SeqIO

def main():
    p = argparse.ArgumentParser(description="Fetch a GenBank or RefSeq record by accession and write FASTA")
    p.add_argument("--acc", required=True, help="Accession, for example NC_038235.1")
    p.add_argument("--out", required=True, help="Output FASTA path")
    p.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""), help="Contact email for NCBI E-utilities")
    args = p.parse_args()

    if not args.email:
        print("Error: set --email or NCBI_EMAIL env var", file=sys.stderr)
        sys.exit(2)

    Entrez.email = args.email
    with Entrez.efetch(db="nucleotide", id=args.acc, rettype="fasta", retmode="text") as handle:
        seq_records = list(SeqIO.parse(handle, "fasta"))
        if not seq_records:
            print("No sequence returned for", args.acc, file=sys.stderr)
            sys.exit(1)
        with open(args.out, "w") as out_fh:
            SeqIO.write(seq_records, out_fh, "fasta")

if __name__ == "__main__":
    main()
