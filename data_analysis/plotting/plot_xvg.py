#!/usr/bin/env python3
import argparse
import glob
import matplotlib.pyplot as plt
import os

def parse_xvg_metadata(filename):
    """Extract title, subtitle, xlabel, ylabel from XVG header."""
    title = ""
    subtitle = ""
    xlabel = "X"
    ylabel = "Y"

    with open(filename, "r") as f:
        for line in f:
            if line.startswith('@'):
                if "title" in line:
                    title = line.split('"')[1]
                elif "subtitle" in line:
                    subtitle = line.split('"')[1]
                elif "xaxis" in line and "label" in line:
                    xlabel = line.split('"')[1]
                elif "yaxis" in line and "label" in line:
                    ylabel = line.split('"')[1]
            elif not line.startswith(('#', '@')):
                break

    return title, subtitle, xlabel, ylabel


def load_xvg_data(filename):
    """Load numeric columns from XVG (skip comments)."""
    x, y = [], []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(('#', '@')) or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return x, y


def plot_xvg(filename, args):
    """Read, plot, and save XVG as JPEG."""
    title, subtitle, xlabel, ylabel = parse_xvg_metadata(filename)
    x, y = load_xvg_data(filename)

    plt.figure(figsize=(7, 5))
    plt.plot(x, y, linewidth=1.5)

    # Title
    plt.title(f"{title}\n{subtitle}" if subtitle else title)

    # Labels (override if provided)
    plt.xlabel(args.xlabel if args.xlabel else xlabel)
    plt.ylabel(args.ylabel if args.ylabel else ylabel)

    # Axis limits
    if args.xlim:
        plt.xlim(args.xlim)
    if args.ylim:
        plt.ylim(args.ylim)

    # Y-scale (linear/log)
    if args.yscale:
        plt.yscale(args.yscale)

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    outname = os.path.splitext(filename)[0] + ".jpg"
    plt.savefig(outname, dpi=300)
    plt.close()

    print(f"Saved plot → {outname}")


def main():
    parser = argparse.ArgumentParser(description="Plot XVG files as JPEG with optional formatting.")
    parser.add_argument("-i", "--input", required=True,
                        help="Input XVG glob pattern, e.g. 'rmsd*.xvg'")
    
    parser.add_argument("--xlim", nargs=2, type=float,
                        help="Set x-axis limits: --xlim xmin xmax")
    parser.add_argument("--ylim", nargs=2, type=float,
                        help="Set y-axis limits: --ylim ymin ymax")
    parser.add_argument("--xlabel", type=str,
                        help="Override x-axis label")
    parser.add_argument("--ylabel", type=str,
                        help="Override y-axis label")
    parser.add_argument("--yscale", type=str, choices=["linear", "log"],
                        help="Set y-axis scale (linear or log)")

    args = parser.parse_args()
    files = sorted(glob.glob(args.input))

    if not files:
        print("No files found for pattern:", args.input)
        return

    for f in files:
        plot_xvg(f, args)


if __name__ == "__main__":
    main()

