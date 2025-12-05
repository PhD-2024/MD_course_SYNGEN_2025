#!/bin/bash

# List of EDR files to process
edr_list=("npt_ber.edr" "npt_ber_posres.edr" "npt.edr" "nvt.edr")

echo "=== Extracting Pressure, Density, Total Energy from EDR files ==="

for edr in "${edr_list[@]}"; do
    base="${edr%.edr}"
    echo "Processing: $edr"

    # --- Pressure ---
    echo "Pressure" | gmx energy -f "$edr" -o "${base}_pressure.xvg" >/dev/null 2>&1

    # --- Density ---
    echo "Density" | gmx energy -f "$edr" -o "${base}_density.xvg" >/dev/null 2>&1

    # --- Total Energy or Potential (whichever exists) ---
    # First try Total Energy
    echo "Total Energy" | gmx energy -f "$edr" -o "${base}_total_energy.xvg" >/dev/null 2>&1

    # If Total Energy wasn't found, use Potential instead
    if [ ! -s "${base}_total_energy.xvg" ]; then
        echo "Potential" | gmx energy -f "$edr" -o "${base}_total_energy.xvg" >/dev/null 2>&1
    fi

    echo " → Output: ${base}_pressure.xvg, ${base}_density.xvg, ${base}_total_energy.xvg"
    echo "--------------------------------------------------------------"

done

echo "Energy extraction complete!"

