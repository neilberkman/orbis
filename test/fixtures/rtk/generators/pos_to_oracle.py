#!/usr/bin/env python3
"""Convert an RTKLIB rnx2rtkp ENU .pos solution into the vendored RTK oracle
JSON shape (matching test/fixtures/rtk/wtzr_wtzz_rtklib_oracle.json).

Usage:
    pos_to_oracle.py POS CONF LABEL DESCRIPTION OUT.json

Truth is the WTZR/WTZZ static antenna baseline (receivers do not move), copied
verbatim from the existing Wettzell oracle, so kinematic-mode and multi-GNSS
oracles share the same physical truth — only the RTKLIB processing changes.
"""
import json
import sys
from datetime import datetime

# Physical truth for the co-located WTZR(base)/WTZZ(rover) pair, identical to
# test/fixtures/rtk/wtzr_wtzz_rtklib_oracle.json ("truth").
TRUTH = {
    "frame": "ENU at WTZR ARP, metres",
    "base_marker_ecef_m": {"x": 4075580.3111, "y": 931854.0543, "z": 4801568.2808},
    "rover_marker_ecef_m": {"x": 4075579.1913, "y": 931853.3696, "z": 4801569.1897},
    "base_antenna_height_m": 0.071,
    "rover_antenna_height_m": 0.284,
    "antenna_baseline_enu_m": {
        "east": -0.41788146461250397,
        "north": 1.5352286147033802,
        "up": 0.08141828505054118,
    },
}
TE = TRUTH["antenna_baseline_enu_m"]["east"]
TN = TRUTH["antenna_baseline_enu_m"]["north"]
TU = TRUTH["antenna_baseline_enu_m"]["up"]

Q_STATUS = {1: "fixed", 2: "float", 3: "sbas", 4: "dgps", 5: "single", 6: "ppp"}

RTKLIB = {"program": "rnx2rtkp", "version": "v2.4.2-p13", "commit": "71db0ff"}


def parse_pos(path):
    epochs = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("2020"):
                continue
            f = line.split()
            # date time e n u Q ns sde sdn sdu sden sdnu sdue age ratio
            d, t = f[0], f[1]
            iso = datetime.strptime(d + " " + t, "%Y/%m/%d %H:%M:%S.%f").isoformat()
            e, n, u = float(f[2]), float(f[3]), float(f[4])
            q, ns, ratio = int(f[5]), int(f[6]), float(f[14])
            epochs.append(
                {
                    "time": iso,
                    "fix_status": Q_STATUS.get(q, str(q)),
                    "q": q,
                    "satellites": ns,
                    "baseline_enu_m": {"east": e, "north": n, "up": u},
                    "ratio": ratio,
                }
            )
    return epochs


def truth_err(ep):
    b = ep["baseline_enu_m"]
    return ((b["east"] - TE) ** 2 + (b["north"] - TN) ** 2 + (b["up"] - TU) ** 2) ** 0.5


def main():
    pos, conf, label, desc, out = sys.argv[1:6]
    epochs = parse_pos(pos)
    fixed = [i for i, e in enumerate(epochs) if e["q"] == 1]
    first_fix = fixed[0] if fixed else None
    last = epochs[-1]
    ref = {
        "label": label,
        "config": conf,
        "source_pos": pos.split("/")[-1],
        "epochs": len(epochs),
        "fixed_epochs": len(fixed),
        "first_fixed_index": first_fix,
        "first_fixed_time": epochs[first_fix]["time"] if first_fix is not None else None,
        "final_status": last["fix_status"],
        "final_ratio": last["ratio"],
        "final_baseline_enu_m": last["baseline_enu_m"],
        "final_truth_error_m": round(truth_err(last), 12),
        "mean_truth_error_m": round(sum(truth_err(e) for e in epochs) / len(epochs), 12),
        "max_truth_error_m": round(max(truth_err(e) for e in epochs), 12),
    }
    doc = {
        "version": "1",
        "description": desc,
        "generator": {"rtklib": RTKLIB, "config": conf},
        "truth": TRUTH,
        "reference": ref,
        "per_epoch": epochs,
    }
    with open(out, "w") as fh:
        json.dump(doc, fh, indent=1)
        fh.write("\n")
    print(
        f"{out}: {ref['fixed_epochs']}/{ref['epochs']} fixed, "
        f"first_fix@{first_fix}, final_err={ref['final_truth_error_m']*1000:.1f}mm"
    )


if __name__ == "__main__":
    main()
