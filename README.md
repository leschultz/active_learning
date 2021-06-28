# Active learning using VASP, LAMMPS, and MLIP-2

The files along with their description follow the convention of <name: > <description> below.

```
vasp_runs
├── README.md
├── examples
│   ├── DFT_INCAR:              VASP INCAR file for DFT
│   ├── KPOINTS:                VASP KPOINTS file for AIMD and DFT
│   ├── MD_INCAR:               VASP INCAR file for AIMD
│   ├── POSCAR:                 User generated POSCAR for VASP
│   ├── POTCAR:                 User generated POTCAR for VASP
│   ├── active_md.sh:           File linking programs and starting active learning
│   ├── curr.mtp:               File containg template MLIP-2 potential
│   ├── md.in:                  LAMMPS input file
│   └── mlip.ini:               MLIP-2 parameter file.
└── scripts
    ├── convert
    │   └── poscar2lammps.awk:  Convert POSCAR to LAMMPS readable file
    ├── run.sh:                 Main script performing active learning
    └── run_types
        ├── gen_aimd.sh:        Generate AIMD run
        ├── gen_dft.sh:         Generate DFT runs
        └── gen_md.sh:          Generate classical MD runs
```
