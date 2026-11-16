#!/usr/bin/env python3
"""Add the MLST sequence type (ST) from mlst/mlst.txt to metadata/meta_gene_data.csv.

Joins on the isolate name (first column of both files) and appends a single
'ST' column to the metadata CSV. Missing STs are written as 'NA', matching the
convention already used for Year/Country in that file.

Usage:
  python add_mlst_to_meta.py
"""
import csv

META = "metadata/meta_gene_data.csv"
MLST = "mlst/mlst.txt"


def main():
    st = {}
    with open(MLST) as fh:
        fh.readline()  # header: name, ST, acsA, ...
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            st[p[0]] = p[1] if p[1] else "NA"

    with open(META) as fh:
        rows = list(csv.reader(fh))
    header, body = rows[0], rows[1:]

    if "ST" not in header:
        header.append("ST")

    n_missing = 0
    for r in body:
        v = st.get(r[0], "NA")
        if v == "NA":
            n_missing += 1
        r.append(v)

    with open(META, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(body)

    print(f"{META}: {len(body)} isolates, added column 'ST' "
          f"({n_missing} with no ST)")


if __name__ == "__main__":
    main()
