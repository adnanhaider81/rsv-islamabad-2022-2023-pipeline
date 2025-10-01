
import sys, pathlib

consensus = pathlib.Path("results/consensus/TestA_01.consensus.fasta")
if not consensus.exists():
    print("Missing consensus:", consensus, file=sys.stderr)
    sys.exit(1)
txt = consensus.read_text()
lines = [l.strip() for l in txt.splitlines() if l.strip()]
if not lines or not lines[0].startswith(">"):
    print("Consensus header missing or malformed", file=sys.stderr); sys.exit(2)
seq = "".join(l for l in lines[1:] if not l.startswith(">"))
if len(seq) < 50:
    print("Consensus too short:", len(seq), file=sys.stderr); sys.exit(3)
print("OK consensus length", len(seq))
