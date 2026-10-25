#!/usr/bin/env python3
"""Six-frame stop-to-stop ORF extraction and translation for the Pfam
filter.

The pfam_scan 1.6 biocontainer (hdfd78af_5) does NOT ship the standalone
'translate' binary pfam_scan.pl -translate requires (verified 2026-10-26:
'Command translate not found' on a real container run), so the module
translates here instead and pfam_scan consumes protein sequences directly.

ORF ids are '<original-id>|orf<N>' so the hit parser can map a domain hit
back to its transcript with an unambiguous separator.

Usage:
    six_frame_orfs.py <input.fasta> <orfs.fasta> [--min-aa 20]
"""
import argparse
import sys

CODONS = {}
_B = "TCAG"
_AA = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
_i = 0
for b1 in _B:
    for b2 in _B:
        for b3 in _B:
            CODONS[b1 + b2 + b3] = _AA[_i]
            _i += 1


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def orfs(frame_seq):
    """Stop-to-stop ORF peptides in one frame (leading M kept as-is)."""
    out = []
    run = []
    for i in range(0, len(frame_seq) - 2, 3):
        aa = CODONS.get(frame_seq[i:i + 3].upper(), "X")
        if aa == "*":
            if run:
                out.append("".join(run))
                run = []
        else:
            run.append(aa)
    if run:
        out.append("".join(run))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("fasta")
    ap.add_argument("out")
    ap.add_argument("--min-aa", type=int, default=20)
    args = ap.parse_args()

    name, seqs = None, []
    records = []

    def flush():
        if name is not None:
            records.append((name, "".join(seqs)))

    with open(args.fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line.strip())
    flush()

    written = 0
    with open(args.out, "w") as out:
        for name, seq in records:
            frames = [seq[i:] for i in range(3)] + [revcomp(seq)[i:] for i in range(3)]
            n = 0
            for frame_seq in frames:
                for pep in orfs(frame_seq):
                    if len(pep) >= args.min_aa:
                        n += 1
                        out.write(f">{name}|orf{n}\n")
                        for i in range(0, len(pep), 60):
                            out.write(pep[i:i + 60] + "\n")
                        written += 1
    print(f"six-frame ORFs (>= {args.min_aa} aa): {written} from {len(records)} sequences")


if __name__ == "__main__":
    sys.exit(main())
