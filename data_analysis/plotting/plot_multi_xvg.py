#!/usr/bin/env python3
import argparse
import glob
import re
import numpy as np
import matplotlib.pyplot as plt


def read_xvg_fast(file):
    """
    Fast reader for 2-column XVG files.
    Extracts subtitle → title → filename as label.
    """
    title = None
    subtitle = None

    # Read metadata lines at the top of the file
    with open(file, "r") as f:
        for line in f:
            if not line.startswith("@"):
                break

            # Title
            if "title" in line:
                m = re.search(r'"(.+)"', line)
                if m:
                    title = m.group(1)

            # Subtitle
            if "subtitle" in line:
                m = re.search(r'"(.+)"', line)
                if m:
                    subtitle = m.group(1)

    # Priority for label
    label = subtitle if subtitle else (title if title else file)

    # Fast numeric read
    data = np.loadtxt(file, comments=("@", "#"))

    x = data[:, 0]
    y = data[:, 1]

    return x, y, label


def main():
    parser = argparse.ArgumentParser(description="Fast multi-XVG plotter with subtitle legend labels.")
    parser.add_argument("-i", "--input", required=True,
                        help="Input glob pattern, e.g. 'rmsd_*g1.xvg'")
    parser.add_argument("-o", "--output", default="plot.png",
                        help="Output image file name (default: plot.png)")
    parser.add_argument("--stride", type=int, default=1,
                        help="Downsample by keeping every Nth point (default: 1)")

    args = parser.parse_args()

    # Find matching XVG files
    files = sorted(glob.glob(args.input))
    if not files:
        print(f"No files found for pattern: {args.input}")
        return

    print(f"Found {len(files)} files:")
    for f in files:
        print("   ", f)

    # ----- Plotting -----
    plt.figure(figsize=(10, 6), dpi=200)

    for f in files:
        x, y, label = read_xvg_fast(f)

        # Downsampling if needed
        if args.stride > 1:
            x = x[::args.stride]
            y = y[::args.stride]

        plt.plot(x, y, label=label, linewidth=1.5)

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Multi-XVG Plot")
    plt.legend()
    plt.tight_layout()

    # Save figure
    plt.savefig(args.output, dpi=300)
    print(f"Saved plot to {args.output}")

    plt.show()


if __name__ == "__main__":
    main()

