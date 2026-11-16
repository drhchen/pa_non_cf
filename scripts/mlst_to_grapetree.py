#!/usr/bin/env python3
"""Convert an MLST table (name, ST, locus1, locus2, ...) into a GrapeTree profile.

GrapeTree's profile rules, as implemented in grapetree/module/MSTrees.py:
  * the header line must begin with '#'            -> 'name' becomes '#Name'
  * the file is tab-delimited, one isolate per row
  * a column called 'ST' (or 'st_id') is allowed; GrapeTree keeps it but
    excludes it from the allele comparison
  * missing alleles must be written as '0', '-' or 'N'. Anything else,
    including 'NA', is read as a real allele and will make isolates match.

The mlst tool (Seemann) additionally emits:
  * 'new'        -> a novel allele not present in the MLST database
  * '5;5'/'4,4'  -> two allele calls for a duplicated locus
Both are handled below.

Usage:
  python mlst_to_grapetree.py mlst/mlst.txt grapetree/mlst.profile
"""
import argparse
import re
import sys

# Tokens that mean "no allele call".
MISSING = {"NA", "N/A", "na", "n/a", "nan", "NaN", "None", "none", "", "-", "?", "."}

# '5;5' or '4,4' -> keep the first call
MULTI = re.compile(r"^\s*(\d+)\s*[;,]")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--missing", default="0",
                    help="code written for missing alleles (default: 0)")
    ap.add_argument("--keep-new", action="store_true",
                    help="keep the 'new' token instead of coding it as missing")
    args = ap.parse_args(argv)

    with open(args.infile) as fh:
        rows = [ln.rstrip("\n") for ln in fh if ln.strip()]

    header = rows[0].split("\t")
    header[0] = "#Name"

    n_na = n_new = n_multi = n_dropped = 0
    out = ["\t".join(header)]
    for line in rows[1:]:
        f = line.split("\t")
        # column 1 is the ST: GrapeTree ignores it, so leave it untouched
        for i in range(2, len(f)):
            v = f[i]
            m = MULTI.match(v)
            if m:
                f[i] = m.group(1)
                n_multi += 1
            elif v in MISSING:
                f[i] = args.missing
                n_na += 1
            elif v.lower() == "new":
                n_new += 1
                if not args.keep_new:
                    f[i] = args.missing
        if all(x == args.missing for x in f[2:]):
            n_dropped += 1
        out.append("\t".join(f))

    with open(args.outfile, "w") as fh:
        fh.write("\n".join(out) + "\n")

    print(f"wrote {args.outfile}: {len(out) - 1} isolates x {len(header) - 1} loci "
          f"(+ name)", file=sys.stderr)
    print(f"  'NA'  -> {args.missing}: {n_na} cells", file=sys.stderr)
    print(f"  'new' -> {'kept' if args.keep_new else args.missing}: {n_new} cells", file=sys.stderr)
    print(f"  'a;b' -> first call  : {n_multi} cells", file=sys.stderr)
    print(f"  isolates with no allele left (GrapeTree drops these): {n_dropped}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
