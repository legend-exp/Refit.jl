"""Write the fake input files of the tutorial in docs/src/tutorial_fccd.md.

- templates.lh5 has the generalised histogram pdf/geds, the spectra of a 60 keV gamma
  line simulated on a grid of FCCD and transition layer (TL) thickness, and the number
  of simulated primaries,
- data.lh5 has the measured histogram hist/geds.

Each gamma deposits its 60 keV at a depth drawn uniformly from 0 to 2 mm below the
surface of the detector. The charge collected from that depth is 0 in the dead layer
(depth below FCCD - TL), rises linearly through the transition layer and is 1 below
the FCCD. So the spectrum is a peak at 60 keV with probability (2 mm - FCCD) / 2 mm and
a flat continuum from 0 to 60 keV with probability TL / 2 mm. Events from the dead
layer have no energy and are not recorded.

The files are written in Python with lgdo, as legend-simflow writes its PDFs. Run with:
python make_files.py
"""

from pathlib import Path

import lh5
import numpy as np
from lgdo import Histogram, Scalar, Struct

out_dir = Path(__file__).parent  # the files are written next to this script
edges = np.arange(0.0, 71.0, 1.0)  # 1 keV bins from 0 to 70 keV


def activeness(depth, fccd, tl):
    """Fraction of the charge collected from the depth, all lengths in mm."""
    if tl == 0:
        return (depth >= fccd).astype(float)
    return np.clip((depth - (fccd - tl)) / tl, 0, 1)


def simulate(rng, n, fccd, tl):
    """The spectrum of n gammas of 60 keV, as a histogram in 1 keV bins."""
    energy = 60 * activeness(2 * rng.random(n), fccd, tl)
    counts, _ = np.histogram(energy[energy > 0], bins=edges)
    return Histogram(counts.astype(float), [(0.0, 70.0, 1.0)])


def mm(value):
    return Scalar(float(value), attrs={"units": "mm"})


# the templates, on a grid of 5 FCCD values times 6 TL values, in mm
rng = np.random.default_rng(1)
n_primaries = 1_000_000
grid = {"fccd": (0.6, 1.4, 0.2), "tl": (0.0, 0.5, 0.1)}  # first, last, step
values = {
    name: np.linspace(first, last, round((last - first) / step) + 1)
    for name, (first, last, step) in grid.items()
}

grid_info = Struct(
    {
        f"axis_{d}": Struct(
            {
                "name": Scalar(name),
                "values": Struct(
                    {"first": mm(first), "last": mm(last), "step": mm(step)}
                ),
            }
        )
        for d, (name, (first, last, step)) in enumerate(grid.items())
    }
)
points = Struct(
    {
        f"fccd_{i}": Struct(
            {
                f"tl_{j}": simulate(rng, n_primaries, fccd, tl)
                for j, tl in enumerate(values["tl"])
            }
        )
        for i, fccd in enumerate(values["fccd"])
    }
)

pdf = Struct({"geds": Struct({"grid_info": grid_info, "grid": points})})
lh5.write(pdf, "pdf", out_dir / "templates.lh5", wo_mode="of")
lh5.write(Scalar(n_primaries), "nr_sim_events", out_dir / "templates.lh5", wo_mode="a")

# the measurement: 20000 gammas, a number known from the activity of the source and
# the livetime
true_fccd, true_tl = 1.05, 0.27  # mm
data = simulate(np.random.default_rng(2), 20_000, true_fccd, true_tl)
lh5.write(Struct({"geds": data}), "hist", out_dir / "data.lh5", wo_mode="of")
