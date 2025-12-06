#!/usr/bin/env python3
import argparse
import glob
import re
import numpy as np
import matplotlib.pyplot as plt


def read_xvg_fast(filename):
    """
    Fast reader for 2-column XVG files.
    Extracts subtitle → title → filename as legend label.
    """
    title = None
    subtitle = None

    with open(filename, "r") as f:
        for line in f:
            if not line.startswith("@"):
                break
            if "subtitle" in line:
                m = re.search(r'"(.+)"', line)
                if m:
                    subtitle = m.group(1)
            if "title" in line:
                m = re.search(r'"(.+)"', line)
                if m:
                    title = m.group(1)

    label = subtitle if subtitle else (title if title else filename)

    data = np.loadtxt(filename, comments=("@", "#"))
    x = data[:, 0]
    y = data[:, 1]

    return x, y, label


def main():
    parser = argparse.ArgumentParser(description="Fast multi-XVG plotter with optional formatting.")
    parser.add_argument("-i", "--input", required=True,
                        help="Input glob pattern, e.g. 'rmsd_*g1.xvg'")
    parser.add_argument("-o", "--output", default="plot.png",
                        help="Output image file name (default: plot.png)")
    parser.add_argument("--stride", type=int, default=1,
                        help="Downsample by keeping every Nth point (default: 1)")

    # Optional formatting
    parser.add_argument("--xlabel", type=str, help="Override x-axis label")
    parser.add_argument("--ylabel", type=str, help="Override y-axis label")
    parser.add_argument("--yscale", type=str, choices=["linear", "log"],
                        help="Y-axis scale (linear or log)")
    parser.add_argument("--xlim", nargs=2, type=float,
                        help="Set x-limits: --xlim xmin xmax")
    parser.add_argument("--ylim", nargs=2, type=float,
                        help="Set y-limits: --ylim ymin ymax")

    args = parser.parse_args()

    files = sorted(glob.glob(args.input))
    if not files:
        print(f"No files found for pattern: {args.input}")
        return

    print(f"Plotting {len(files)} files:")
    for f in files:
        print("  ", f)

    # ----- Plot -----
    plt.figure(figsize=(10, 6), dpi=200)

    for f in files:
        x, y, label = read_xvg_fast(f)

        if args.stride > 1:
            x = x[::args.stride]
            y = y[::args.stride]

        plt.plot(x, y, label=label, linewidth=1.5)

    # Labels
    plt.xlabel(args.xlabel if args.xlabel else "x")
    plt.ylabel(args.ylabel if args.ylabel else "y")
    plt.title("Multi-XVG Plot")

    # Axis limits
    if args.xlim:
        plt.xlim(args.xlim)
    if args.ylim:
        plt.ylim(args.ylim)

    # Scale
    if args.yscale:
        plt.yscale(args.yscale)

    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(args.output, dpi=300)
    print(f"Saved figure → {args.output}")


if __name__ == "__main__":
    main()

