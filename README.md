# cleaving

Lammps Package for cleaving calculations USER-CLEAVING

The package contains all the modified potentials and additional computes/fixes needed by the cleaving model to calculate the Surface Free Energy in a system using Thermodynamic Integration.

This package is presented as part of the consortium Simulation Interfacial Free Energy Techniques (SIFT) which include the packages:

  - [MOLD](https://github.com/AndresRTejedor/Mold)

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
