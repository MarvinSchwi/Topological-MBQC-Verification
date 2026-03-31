# Topological-MBQC-Verification

This repository contains source code for the algorithmic verification of topological MBQC circuits, as described in the paper "[A 3D lattice defect and efficient computations in topological MBQC](https://quantum-journal.org/papers/q-2026-02-06-1997/)". Please cite this paper as:

`Gabrielle Tournaire, Marvin Schwiering, Robert Raussendorf, Sven Bachmann, Quantum 10, 1997 (2026)`

The corresponding BibTex entry is:

```bibtex
@article{TournaireSchwieringRaussendorfBachmann2026,
  author    = {Tournaire, Gabrielle and Schwiering, Marvin and Raussendorf, Robert and Bachmann, Sven},
  year      = {2026},
  title     = {A {{3D}} Lattice Defect and Efficient Computations in Topological {{MBQC}}},
  journal   = {Quantum},
  volume    = {10},
  pages     = {1997},
  publisher = {Verein zur F{\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften},
  doi       = {10.22331/q-2026-02-06-1997},
  url       = {https://quantum-journal.org/papers/q-2026-02-06-1997/}
}
```

## License

This project is licensed under the GNU General Public License v3 (GPLv3). See the [LICENSE](./LICENSE) file for details.

## Installation

To execute the verifier, install conda and run `conda env create [-n CUSTOMENVNAME] -f environment.yml` where `CUSTOMENVNAME` (if specified) is replaced with your desired environment name. To uninstall, run `conda remove -n CUSTOMENVNAME --all`.

## Code Overview

Running the visualizer is done via `python app.py possibly/relative/path/to/save.json [--fullscreen]`. The verification is currently not accessible via the GUI but only via the command line (see below) or python imports.

## Visualizer Controls

Camera movement:
- Hold the left mouse button to move the lattice in the plane orthogonal to your line of sight.
- Hold the middle mouse button (alternatively: both the left and right mouse button) to rotate the lattice.
- Hold the right mouse button to move the lattice along your line of sight.

Further key bindings:
- `L`: Toggle the visibility of the lattice (primal/dual/none).
- `D`: Toggle the visibility of the defects (all/none/primal/dual).
- `T`: Toggle the visibility of the target vector (none/1st target vector (if given)/2nd target vector (if given)/...).
- `S`: Toggle the visibility of the correlation surfaces for the shown target vector, if specified (on/off).
- `O`: Open another lattice save.

## Verification

Verification can be accessed via the `verification` method of a `Lattice` object. `return_each = True` will produce a verification for each target vector.

```python
>>> from main import Lattice
>>> lattice = Lattice.load(save_name = "CNOT_Correct") # A Valid Topological Circuit
>>> lattice.verification() # Will Return True
True
>>> lattice = Lattice.load(save_name = "IdentityPrimal_Straight_IncorrectUpDown") # A (Partially) Invalid Topological Circuit
>>> lattice.verification() # Will Return False
False
>>> lattice.verification(return_each = True) # Will Show Individual Results for Target Vectors
{'tv00': True, 'tv01': False}
```
