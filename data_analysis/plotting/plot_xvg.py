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
                # First data line reached → stop parsing metadata
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


def plot_xvg(filename):
    """Read, plot, and save XVG as JPEG."""
    title, subtitle, xlabel, ylabel = parse_xvg_metadata(filename)
    x, y = load_xvg_data(filename)

    plt.figure(figsize=(7, 5))
    plt.plot(x, y, linewidth=1.5)

    plt.title(f"{title}\n{subtitle}" if subtitle else title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)

    outname = os.path.splitext(filename)[0] + ".jpg"
    plt.tight_layout()
    plt.savefig(outname, dpi=300)
    plt.close()

    print(f"Saved plot → {outname}")


def main():
    parser = argparse.ArgumentParser(description="Plot XVG files as JPEG")
    parser.add_argument(
        "-i", "--input", required=True, 
        help="Input XVG glob pattern, e.g. 'rmsd*.xvg'"
    )
    
    args = parser.parse_args()
    files = sorted(glob.glob(args.input))

    if not files:
        print("No files found for pattern:", args.input)
        return

    for f in files:
        plot_xvg(f)


if __name__ == "__main__":
    main()

