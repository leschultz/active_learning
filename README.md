# Active learning using VASP, LAMMPS, and MLIP-2

## Dependencies

```
ONEAPI from Intel
mlip-2
VASP.5*
LAMMPS
LAMMPS mlip-2 interface
``` 
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
│   ├── DFT_INCAR:              VASP INCAR file for DFT
│   ├── KPOINTS:                VASP KPOINTS file for AIMD and DFT
│   ├── MD_INCAR:               VASP INCAR file for AIMD
│   ├── POSCAR:                 User generated POSCAR for VASP
│   ├── POTCAR:                 User generated POTCAR for VASP
│   ├── active_md.sh:           File linking programs and starting active learning
│   ├── curr.mtp:               File containg template MLIP-2 potential
│   ├── md.in:                  LAMMPS input file following desired MD run
│   └── mlip.ini:               MLIP-2 parameter file
└── scripts
    ├── convert
    │   └── poscar2lammps.awk:  Convert POSCAR to LAMMPS readable file
    ├── run.sh:                 Main script performing active learning
    └── run_types
        ├── gen_aimd.sh:        Generate AIMD run
        ├── gen_dft.sh:         Generate DFT runs
        └── gen_md.sh:          Generate classical MD runs
```

Outputs are described below for an example run and follow the convention of \<name: \> \<description\>.

```
run
├── aimd:         The AIMD initial training data
├── md_dft:       Beginning of active learning loop
│   ├── 0
│   │   ├── md:   Run classical MD with initial potential
│   │   └── dft:  Run DFT for flagged configurations
│   └── 1
│       └── md:   Final classical MD (no flagged configurations and can be n > 1)
└── potential:    Folder containing final potential
```
