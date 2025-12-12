# TB-tbG: Tight-binding model solver for twisted bilayer graphene (MPI/Serial version)

## 1. Features

TB-tbG is a Fortran-based program designed for computing electronic properties of twisted bilayer graphene (tbG) using a tight-binding (TB) model. The program supports both MPI-parallel and serial execution.

Key capabilities include:

- **Band structure calculation**: Compute the electronic band structure along high-symmetry k‑paths.
- **Electron Energy Loss Spectroscopy (EELS)**: Calculate the loss function for a given momentum transfer (q) range.
- **MPI‑parallelized computations**: Efficiently distribute k‑point workloads across multiple processes.
- **Flexible input system**: Configure calculations via a simple text file (`input.in`).
- **VASP‑compatible structure input**: Read atomic positions and lattice vectors from a `POSCAR` file.

The code is written in modern Fortran with modular design, separating physical models, I/O, solvers, and parallelization.

## 2. Compilation Instructions

### Prerequisites
- Intel® Fortran Compiler (ifort/ifx) with MPI support (mpiifx)
- Intel® Math Kernel Library (MKL) for linear algebra operations
- Make utility

### Build Steps
1. Clone the repository:
   ```bash
   git clone https://github.com/gawcista/TB-tbG.git
   cd TB-tbG
   ```

2. Compile the program:
   ```bash
   make
   ```
   This will create the executable `tbsolver` in the current directory.

   The Makefile uses the following default settings:
   - Compiler: `mpiifx`
   - Flags: `-O2 -qmkl -module build -Ibuild`
   - Linking: `-qmkl`
   - Object and module files are placed in the `build/` directory.

3. Clean the build:
   ```bash
   make clean
   ```
   Removes the `build/` directory and the executable.

## 3. Quick‑Start Usage

### Preparing Input Files
1. **Structure file**: Provide a VASP‑format `POSCAR` file (default name `cont.vasp`). It must contain the lattice vectors and atomic positions of the twisted bilayer graphene system.

2. **Input configuration**: Create an `input.in` file in the same directory. See Section 4 for a detailed description of all available keywords.

   A minimal example for a band‑structure calculation (`JOB B`):
   ```plaintext
   JOB B
   FPOSCAR cont.vasp
   FEIGEN tb_band.dat
   FKPOINTS KPOINTS.band
   ```

### Running the Program
- **MPI parallel run** (recommended for large calculations):
  ```bash
  mpirun -n 4 ./tbsolver
  ```
  Replace `4` with the desired number of MPI processes.

- **Serial run** (if compiled with MPI but using one process):
  ```bash
  ./tbsolver
  ```

### Output Files
Depending on the job type, the program produces:
- **Band structure**: File specified by `FEIGEN` (default `tb_band.dat`) contains k‑points along the path and corresponding eigenvalues.
- **EELS results**: Loss function and related quantities are written to files with names derived from the `QTAG` keyword.

## 4. Input Keywords (input.in)

All keywords are read by the parser in `src/parser.f90`. Lines starting with `!`, `#`, or `/` are treated as comments and ignored. The order of keywords is arbitrary.

| Keyword | Arguments | Description |
|---------|-----------|-------------|
| **JOB** | `B` or `E` | Type of calculation:<br>`B` – Band structure<br>`E` – Electron Energy Loss Spectroscopy (EELS) |
| **FPOSCAR** | `<filename>` | Path to the VASP‑format POSCAR file (default: `cont.vasp`) |
| **FEIGEN** | `<filename>` | Output file for band eigenvalues (default: `tb_band.dat`). Also enables band‑structure writing. |
| **FKPOINTS** | `<filename>` | VASP‑format file containing k‑point path definitions (default: `KPOINTS`). If present, enables k‑point selection. |
| **IBAND** | `<iband_start> <iband_end>` | Select a specific range of bands to compute (1‑based indices). Enables band selection. |
| **NCACHE** | `<integer>` | Cache size for intermediate arrays (default: 0 – automatic). This number determines how many numbers of k points would be calculated together on each rank. Less number of `ncache` reduce memory cost on each worker rank but brings higher MPI communication cost.|
| **OMEGA** | `<nomega> <omega_i> <omega_f>` | Define the energy (frequency) mesh for the loss function:<br>`nomega` – number of points<br>`omega_i` – starting energy (eV)<br>`omega_f` – ending energy (eV) |
| **QPOINT** | `<nqpts> <qx1 qy1 qz1> <qx2 qy2 qz2>` | Define a q‑point path for EELS:<br>`nqpts` – number of q‑points along the line<br>`(qx1, qy1, qz1)` – start point in reciprocal fractional coordinates<br>`(qx2, qy2, qz2)` – end point |
| **QTAG** | `<character>` | Single‑character label for the q‑point set (default: `D`). `D` means the q vectors are provided in direct (fractional) coordinates; `C` specifies that q vectors are provided in a Cartesian coordinate system |
| **EFERMI** | `<energy>` | Fermi energy in eV (default: 0.0). |
| **TEMPERATURE** | `<T>` | Temperature in Kelvin (default: 300). |
| **DELTA** | `<broadening>` | Broadening parameter (eV) for spectral functions (default: 1e‑3). |

### Example input.in for EELS calculation
```plaintext
JOB E
FPOSCAR cont.vasp
OMEGA 51 0.0 0.5
QPOINT 20 0.0 0.0 0.0 0.1 0.1 0.0
QTAG D
EFERMI -0.05
TEMPERATURE 300
DELTA 0.001
```

## 5. Code Structure

- `src/constants.f90` – Physical constants, global variables, and default settings.
- `src/parser.f90` – Input file parser (keywords listed above).
- `src/tbmodel.f90` – Tight‑binding Hamiltonian construction.
- `src/solver.f90` – Serial solvers for band structure and EELS.
- `src/mpi_solver.f90` – MPI‑parallelized solvers.
- `src/eels.f90` – EELS‑specific routines.
- `src/ioutils.f90` – I/O utilities.
- `src/utils.f90` – General utility functions.
- `src/main.f90` – Main program, initializes MPI and dispatches jobs.
- `src/interface.f90` – Module interfaces.
- `src/types.f90` – Derived data types.

## 6. License and Citation

This software is provided for academic and research use. If you use TB‑tbG in your work, please cite the relevant publications.

---
*For questions or issues, please contact the maintainers or open an issue on the project repository.*
