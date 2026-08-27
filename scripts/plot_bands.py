#!/usr/bin/env python3
import sys
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def parse_bulkek(fname):
    """gnuplot-block format: each band is a block of k,E pairs separated by blank lines."""
    bands = []
    cur = []
    for line in open(fname):
        line = line.strip()
        if line.startswith("#"):
            continue
        if line == "":
            if cur:
                bands.append(np.array(cur))
                cur = []
            continue
        parts = line.split()
        cur.append((float(parts[0]), float(parts[1])))
    if cur:
        bands.append(np.array(cur))
    return bands

def parse_gnu_ticks(fname):
    text = open(fname).read()
    m = re.search(r'set xtics \((.*?)\)\n', text)
    entries = re.findall(r'"([^"]*)"\s*([\d.]+)', m.group(1))
    labels = [e[0].strip() for e in entries]
    labels = ["$\\Gamma$" if l == "G" else l for l in labels]
    pos = [float(e[1]) for e in entries]
    return labels, pos

def plot_one(dat_file, gnu_file, title, outfile, highlight_bands=(17, 18)):
    bands = parse_bulkek(dat_file)
    labels, ticks = parse_gnu_ticks(gnu_file)

    fig, ax = plt.subplots(figsize=(9, 6))
    for i, b in enumerate(bands):
        n = i + 1  # 1-indexed band number
        if n in highlight_bands:
            color = "tab:red" if n == highlight_bands[0] else "tab:orange"
            ax.plot(b[:, 0], b[:, 1], lw=1.6, color=color, zorder=5,
                    label=f"band {n}")
        else:
            ax.plot(b[:, 0], b[:, 1], lw=0.6, color="0.35", zorder=2)

    for t in ticks:
        ax.axvline(t, color="0.7", lw=0.6, zorder=1)
    ax.axhline(0, color="0.5", lw=0.5, ls="--", zorder=1)

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlim(ticks[0], ticks[-1])
    ax.set_ylabel("Frequency (THz)")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    print("wrote", outfile)

if __name__ == "__main__":
    dat, gnu, title, out = sys.argv[1:5]
    plot_one(dat, gnu, title, out)
