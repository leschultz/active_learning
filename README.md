# Active learning using VASP, LAMMPS, and MLIP-2

## Dependencies

```
ONEAPI from Intel
mlip-2
VASP.5*
LAMMPS
LAMMPS mlip-2 interface
``` 

Python 3.6.8 was used with the packages listed under python_requirements.txt. Not all package are needed for functionality, but installing them ensures functionality.

## To run
Execute the following command:

```
./active_md.sh
```

## Tool description
The files along with their description follow the convention of \<name: \> \<description\> below. Make sure to edit files in the examples directory in your own run. All files in the examles directory are needed along with their naming convention. The order of listed elements in the POTCAR must match the order of elements in the POSCAR.

```
active_learning
├── README.md:                  Self-explanatory
├── examples
│   ├── aimd:                   The initial set of AIMD runs for fitting
│   └── active:                 Where input files live for active learning
│       ├── INCAR:              VASP INCAR file for DFT
│       ├── KPOINTS:            VASP KPOINTS file for AIMD and DFT
│       ├── POSCAR:             User generated POSCAR for VASP
│       ├── POTCAR:             User generated POTCAR for VASP
│       ├── active_md.sh:       File linking programs and starting active learning
│       ├── curr.mtp:           File containg template MLIP-2 potential
│       ├── md.in:              LAMMPS input file following desired MD run
│       ├── md_test.in:         LAMMPS input file defining configuations from MD to use in AIMD test
│       └── mlip.ini:           MLIP-2 parameter file
├── python_requirements.txt
└── scripts
    ├── convert
    │   └── poscar2lammps.awk:  Convert POSCAR to LAMMPS readable file
    ├── fit.sh:                 Main script performing active learning
    ├── funcs
    │   ├── convert.py          Convert LAMMPS dump files to POSCARS
    │   ├── gen_parity.sh       Workflow for generating test configurations
    │   └── parity.py           Python parity plots
    └── run_types
        ├── gen_dft.sh:         Generate DFT runs
        └── gen_md.sh:          Generate classical MD runs
```

Outputs are described below for an example run and follow the convention of \<name: \> \<description\>. You need to point to the path where the initial AIMD runs are in the active_md.sh file.

```
run
├── md_dft
│   ├── 0
│   │   ├── dft:      Run DFT for flagged configurations
│   │   ├── md:       Run classical MD with fit potential
│   │   ├── retrain:  Retrained potential
│   │   └── test
│   │       ├── aimd: AIMD from test configurations
│   │       ├── md:   Generate test configurations
│   │       └── ml:   Parity plot of forces and energies
│   └── 1
│       └── md:       Final classical MD (no flagged configurations and can be n > 1)
│       └── test
│           ├── aimd: AIMD from test configurations
│           ├── md:   Generate test configurations
│           └── ml:   Parity plot of forces and energies
└── potential         Folder containing final potential
```
