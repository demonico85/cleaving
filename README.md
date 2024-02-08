# CLEAVING

[![DOI](https://joss.theoj.org/papers/10.21105/joss.05886/status.svg)](https://doi.org/10.21105/joss.05886)

The LAMMPS package CLEAVING contains all the modified potentials and additional computes/fixes needed by the cleaving model to use Thermodynamic Integration to calculate the Surface Free Energy of a given system.

This package is presented as part of the consortium Simulation Interfacial Free Energy Techniques (SIFT) which include the packages:

  - [MOLD](https://github.com/AndresRTejedor/Mold)
  - [Einstein IFE](https://github.com/syeandel/Einstein_IFE)

![SIFT group](./docs/figs/final_logo.png)

## Documentation

The documentation can be browsed online at [this link](https://demonico85.github.io/cleaving/). The documentation can also be generated locally through the following procedure, which requires a working python3 installation and an active internet connection

```bash
cd docs
pip3 install -r docs_requirements.txt # use pip3 or pip depending on your local setup
make html
```

At the end of the generation point your browser to `build/html/index.html` to browse the docs.

NOTE: The repository contains a directory `/utils` which includes a miscellanea of programs useful for  cleaving calculations (e.g., extracting and averaging useful quantities). These small programs are given without description and without warranty  

## Contributing

If you want to add functionalities (or fix issues) to CLEAVING, you are more than welcome to do so. Feel free to fork the repo and open a pull request: we will do our best to review it. Note that if we accept your changes you will be asked to update the documentation to reflect what you have done.

If you have any questions, find a bug or an issue or want to propose a new feature, please do so by creating a [new issue](https://github.com/demonico85/cleaving/issues/new/choose). Pick the issue template that makes the most sense and edit it to add as much detail as possible.

## Citing CLEAVING

[![DOI](https://joss.theoj.org/papers/10.21105/joss.05886/status.svg)](https://doi.org/10.21105/joss.05886)

Please cite the following paper for any work that uses the CLEAVING package:

* Pasquale et al., (2024). *CLEAVING: a LAMMPS package to compute surface free energies*. Journal of Open Source Software, 9(94), 5886, [https://doi.org/10.21105/joss.05886](https://doi.org/10.21105/joss.05886)
