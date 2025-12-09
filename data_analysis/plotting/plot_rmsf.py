#!/usr/bin/env python3
import argparse
import glob
import matplotlib.pyplot as plt
import os

def parse_xvg_metadata(filename):
    """Extract title, subtitle, xlabel, ylabel from XVG header."""
    title = ""
    subtitle = ""
    xlabel = "Residue Index"
    ylabel = "RMSF (nm)"

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


def load_rmsf_xvg(filename):
    """Load RMSF data (residue index, RMSF value) from XVG."""
    residues = []
    rmsf = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(('#', '@')) or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                residues.append(int(parts[0]))
                rmsf.append(float(parts[1]))
    return residues, rmsf


def plot_rmsf_bar(filename):
    """Plot RMSF as a bar chart."""
    title, subtitle, xlabel, ylabel = parse_xvg_metadata(filename)
    residues, rmsf_values = load_rmsf_xvg(filename)

    plt.figure(figsize=(12, 4))
    plt.bar(residues, rmsf_values, width=0.9, color="steelblue", edgecolor="black")

    plt.title(f"{title}\n{subtitle}" if subtitle else title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.tight_layout()

    outname = os.path.splitext(filename)[0] + "_bar.jpg"
    plt.savefig(outname, dpi=300)
    plt.close()

    print(f"Saved bar plot → {outname}")


def main():
    parser = argparse.ArgumentParser(description="Plot RMSF XVG files as bar plots")
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input XVG glob pattern, e.g. 'rmsf*.xvg'"
    )

    args = parser.parse_args()
    files = sorted(glob.glob(args.input))

    if not files:
        print("No files found for pattern:", args.input)
        return

    for f in files:
        plot_rmsf_bar(f)


if __name__ == "__main__":
    main()

