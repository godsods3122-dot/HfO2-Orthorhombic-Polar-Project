#!/usr/bin/env python3
"""Rewrite pn.in's KPATH_BULK to a specific label sequence, using each
structure's own seekpath-derived point coordinates, and turn on BulkBand_calc."""
import re
import sys

base_pnin, outfile = sys.argv[1], sys.argv[2]
SEQ = ["G", "X", "S", "Y", "G", "Z", "U", "R", "T", "Z"]

text = open(base_pnin).read()

m = re.search(r"KPATH_BULK\n(\d+)\n((?:.*\n)+?)\n", text)
body = m.group(2)

points = {}
for line in body.strip().split("\n"):
    parts = line.split()
    # label k1 k2 k3 label2 k1_2 k2_2 k3_2
    points[parts[0]] = tuple(parts[1:4])
    points[parts[4]] = tuple(parts[5:8])

segs = list(zip(SEQ[:-1], SEQ[1:]))
new_body = ""
for a, b in segs:
    ka, kb = points[a], points[b]
    new_body += f"{a:4s}    {ka[0]:>20s}    {ka[1]:>20s}    {ka[2]:>20s}    {b:4s}    {kb[0]:>20s}    {kb[1]:>20s}    {kb[2]:>20s}\n"

new_block = f"KPATH_BULK\n{len(segs)}\n{new_body}\n"
text = text[:m.start()] + new_block + text[m.end():]

text = text.replace("BulkBand_calc         = F", "BulkBand_calc         = T")
text = text.replace(
    "  Gap_threshold = 0.00001 !>threshold for GapCube output \n/",
    "  Gap_threshold = 0.00001 !>threshold for GapCube output \n"
    "  LOTO_method = 'phonopy'\n/",
)

open(outfile, "w").write(text)
print(f"wrote {outfile}")
for a, b in segs:
    print(f"  {a} {points[a]} -> {b} {points[b]}")
