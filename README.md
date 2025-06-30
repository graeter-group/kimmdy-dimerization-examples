# kimmdy-dimerization-examples
Examples demonstrating the usage of the kimmdy-dimerization plugin including analysis scripts for quantum yield determination.

- `data` contains the extracted data such as rates, distances, angles (from KIMMDY or gmx angle/distance) or grid search results for parameter fitting
- `graphics` contains Images obtained from the raw data used in the manuscript (or the SI) 
- `input_structures` contains the pdb input files for all structures as well as the equilibrated structures that were used as Input for the KIMMDY runs
- `parametrization` contains everything regarding charge and residue parametrization of the additional dimer residue - including the modified forcefield (which is based on OL21).
- `scripts_and_setup` contains all scripts / input files (besides structures) used to setup GROMACS or KIMMDY simulations and analyize them


Sample names different from the manuscript: `ds_left` and `ds_right` stand for `left only` or `right only`. `ds_both` in equilibration is the precursor for `left 1st`, `left 2nd`, `right 1st` and `right 2nd`. `ds_both_pre` in KIMMDY runs contains `left 1st` and `right 1st` (measure simulaniously, see SI) and ds_both_post in KIMMDY runs contains `left 2nd` and `right 2nd` (each simulation 'decided' for one of them, see SI). 