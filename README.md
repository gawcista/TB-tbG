# TB-tbG: Tight-binding model solver for twisted bilayer graphene (MPI/Serial version)

## 1. Features

TB-tbG is a Fortran-based program designed for computing electronic properties of twisted bilayer graphene (tbG) using a tight-binding (TB) model. The program supports both MPI-parallel and serial execution.

Key capabilities include:

- **Band structure calculation**: Compute the electronic band structure along high-symmetry k‑paths.
- **Density of States (DOS) calculation**: Compute the density of states using Gaussian broadening on k‑point meshes.
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
   - Flags: `-O2 -qmkl -qopenmp -module build -Ibuild -xHost`
   - Linking: `-qmkl -qopenmp`
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
   NCACHE 1
   ```

   A minimal example for a DOS calculation (`JOB D`):
   ```plaintext
   JOB D
   FPOSCAR cont.vasp
   FKPOINTS KPOINTS
   SIGMA 0.1
   NEDOS 1001
   ERANGE -5.0 5.0
   FDOS tb_dos.dat
   DOS_NORMALIZE .true.
   GAUSSIAN_CUTOFF 5.0
   NCACHE 1
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

By default, the program reads `input.in`. A different input file can be passed as the first command-line argument:
```bash
./tbsolver input.band
mpirun -n 4 ./tbsolver input.dos
```

For memory‑limited runs, keep `NCACHE` small. The default `NCACHE 0` uses the lowest-memory automatic mode, processing one k point per group. For pure MPI runs, it is often useful to avoid thread oversubscription:
```bash
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
mpirun -n 4 ./tbsolver
```

### Output Files
Depending on the job type, the program produces:
- **Band structure**: File specified by `FEIGEN` (default `tb_band.dat`) contains k‑points along the path and corresponding eigenvalues.
- **DOS results**: File specified by `FDOS` (default `tb_dos.dat`) contains energy and DOS columns. With `DOS_NORMALIZE .true.`, the DOS is divided by the number of k points and is reported in states/eV per cell.
- **EELS results**: Loss function and related quantities are written to files with names derived from the `QTAG` keyword.

## 4. Input Keywords (input.in)

All keywords are read by the parser in `src/parser.f90`. Lines starting with `!`, `#`, or `/` are treated as comments and ignored. The order of keywords is arbitrary.

| Keyword | Arguments | Description |
|---------|-----------|-------------|
| **JOB** | `B`, `D`, or `E` | Type of calculation:<br>`B` – Band structure<br>`D` – Density of States (DOS)<br>`E` – Electron Energy Loss Spectroscopy (EELS) |
| **FPOSCAR** | `<filename>` | Path to the VASP‑format POSCAR file (default: `cont.vasp`) |
| **FEIGEN** | `<filename>` | Output file for band eigenvalues (default: `tb_band.dat`). Also enables band‑structure writing. |
| **FDOS** | `<filename>` | Output file for DOS data (default: `tb_dos.dat`). |
| **FKPOINTS** | `<filename>` | VASP‑format file containing k‑point definitions (default: `KPOINTS`). Line mode is used for band paths; Gamma or Monkhorst-Pack mesh mode is used for DOS and EELS meshes. |
| **FWANNIER** | `<filename>` | Path to the Wannier90 input file (default: `wannier90`). If present, sets model type to Wannier (`wan`). |
| **IBAND** | `<iband_start> <iband_end>` | Select a specific range of bands to compute (1‑based indices). Enables band selection. |
| **NCACHE** | `<integer>` | Cache size for intermediate arrays (default: 0 – automatic low-memory mode, one k point per group). This number determines how many numbers of k points would be calculated together on each rank. Less number of `ncache` reduce memory cost on each worker rank but brings higher MPI communication cost.|
| **NEDOS** | `<integer>` | Number of DOS energy grid points (default: 1001). Must be greater than 1. |
| **ERANGE** | `<emin> <emax>` | Energy range in eV for DOS (default: `-1.0 1.0`). The upper bound must be greater than the lower bound. |
| **SIGMA** | `<broadening>` | Gaussian broadening for DOS in eV (default: 0.01). Must be positive. |
| **GAUSSIAN_CUTOFF** | `<factor>` | Optional DOS Gaussian cutoff in units of `SIGMA` (default: 0.0, disabled). A positive value limits the energy window for faster DOS accumulation. |
| **DOS_NORMALIZE** | `<logical>` | Whether to divide DOS by the number of k points (default: `.true.`). This gives states/eV per cell for the selected band set. |
| **EELSMODE** | `<integer>` | Select EELS calculation mode: 0 for full calculation, 1 for intraband only, 2 for interband only. |
| **OMEGA** | `<nomega> <omega_i> <omega_f>` | Define the energy (frequency) mesh for the loss function:<br>`nomega` – number of points<br>`omega_i` – starting energy (eV)<br>`omega_f` – ending energy (eV) |
| **QPOINT** | `<nqpts> <qx1 qy1 qz1> <qx2 qy2 qz2>` | Define a q‑point path for EELS:<br>`nqpts` – number of q‑points along the line<br>`(qx1, qy1, qz1)` – start point in reciprocal fractional coordinates<br>`(qx2, qy2, qz2)` – end point |
| **QTAG** | `<character>` | Single‑character label for the q‑point set (default: `D`). `D` means the q vectors are provided in direct (fractional) coordinates; `C` specifies that q vectors are provided in a Cartesian coordinate system |
| **EFERMI** | `<energy>` | Fermi energy in eV (default: 0.0). |
| **TEMPERATURE** | `<T>` | Temperature in Kelvin (default: 300). |
| **DELTA** | `<broadening>` | Broadening parameter (eV) for spectral functions (default: 1e‑3). |
| **CALC_IQR** | `<logical>` | Whether to calculate the imaginary part of the response function (default: `.false.`). |
| **TMPI** | `<logical>` | Enable/disable MPI timer profiling (default: `.true.`). |
| **TDEBUG** | `<logical>` | Enable/disable debug timer profiling (default: `.false.`). |
| **MODEL** | `<string>` | Model type selector: `tbg` for tight-binding graphene (default) or `wan` for Wannier model. |
| **WRITEM** | `<logical>` | Whether to write the IPF matrix to file (default: `.false.`). |
| **FEMATRIX** | `<filename>` | Output file for IPF matrix data (default: `tb_eels.mat`). |
| **TB_EPSILON** | `<value>` | Relative dielectric constant used in EELS Coulomb interaction (default: 4.9). |
| **TB_ONSITE** | `<value>` | Tight-binding onsite energy in eV (default: -0.78). |
| **TB_GAMMA0** | `<value>` | Intralayer nearest-neighbor hopping parameter in eV (default: 2.7). |
| **TB_GAMMA1** | `<value>` | Interlayer coupling parameter in eV (default: 0.48). |

### KPOINTS formats

The `readKPOINTS` routine supports two VASP-like formats:

- **Line mode** for band paths (`JOB B`), where the file provides path endpoints.
- **Automatic mesh mode** for DOS and EELS (`JOB D` and `JOB E`), for example:
  ```plaintext
  Gamma-centered mesh for DOS
  0
  G
  12 12 1
  0.0 0.0 0.0
  ```

Explicit weighted k-point lists are not currently used by the calculation routines.

### Example input.in for band calculation
```plaintext
JOB B
FPOSCAR cont.vasp
FKPOINTS KPOINTS.band
FEIGEN tb_band.dat
IBAND 1 20
NCACHE 1
```

Example `KPOINTS.band`:
```plaintext
High-symmetry path
40
Line-mode
Reciprocal
0.000000 0.000000 0.000000  ! Gamma
0.333333 0.333333 0.000000  ! K
0.333333 0.333333 0.000000  ! K
0.500000 0.000000 0.000000  ! M
0.500000 0.000000 0.000000  ! M
0.000000 0.000000 0.000000  ! Gamma
```

### Example input.in for DOS calculation
```plaintext
JOB D
FPOSCAR cont.vasp
FKPOINTS KPOINTS.dos
FDOS tb_dos.dat
IBAND 1 20
SIGMA 0.05
NEDOS 2001
ERANGE -3.0 3.0
DOS_NORMALIZE .true.
GAUSSIAN_CUTOFF 5.0
NCACHE 1
```

Example `KPOINTS.dos`:
```plaintext
Gamma-centered mesh for DOS
0
G
60 60 1
0.0 0.0 0.0
```

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
