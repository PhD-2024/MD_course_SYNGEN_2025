#!/usr/bin/env python3
import matplotlib.pyplot as plt

def load_xvg(filename):
    """Load XVG file and return x, y1, y2 lists."""
    x, y1, y2 = [], [], []
    with open(filename) as f:
        for line in f:
            if line.startswith(('#','@')) or not line.strip():
                continue
            cols = line.split()
            # lifetime.xvg always has 3 numeric columns: t, p(t), t·p(t)
            x.append(float(cols[0]))
            y1.append(float(cols[1]))
            if len(cols) > 2:
                y2.append(float(cols[2]))
    return x, y1, y2


# -------------------- Load data --------------------
x, p_t, tp_t = load_xvg("hb_lifetime.xvg")


# -------------------- Plot --------------------
plt.figure(figsize=(7,5))

plt.plot(x, p_t, 'bo', label="p(t) — autocorrelation")
plt.plot(x, tp_t, 'rd', lw=2, label="t·p(t)")

plt.yscale("log")                         # LOG SCALE
plt.xlabel("Time (ps or ns)")
plt.ylabel("Autocorrelation C(t)")
plt.title("H-Bond Lifetime Autocorrelation")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("hb_lifetime_logscale.jpg", dpi=300)
plt.close()

print("Saved: hb_lifetime_logscale.jpg")

