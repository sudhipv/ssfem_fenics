# Scalable Solver for Spectral Stochastic FEM using Domain Decomposition with One/Two-Level Preconditioner

## Overview
This solver targets stochastic steady-state diffusion and static linear elasticity (2D & 3D) problems using domain decomposition methods with scalable one- and two-level preconditioners. The `master` branch keeps a clean copy of the code, while the `cedar`, `graham`, `niagara`, `beluga`, and other branches contain machine-specific adjustments tested on their respective clusters.

## Core Components
1. **FEniCS** – deterministic FE matrix/vector assembly.
2. **PETSc** – sparse storage and linear algebra for stochastic FE matrices.
3. **MPI** – parallel execution of PETSc local routines.

## Tested Software Stack
- **FORTRAN**: GNU GCC
- **MPI**: Open MPI
- **Gmsh**: 3.0.4 (use this exact version)
  - Binary: http://gmsh.info/bin/
  - Source: https://gitlab.onelab.info/gmsh/gmsh/blob/master/README.txt
- **PETSc**: 3.7.5
  - Local install: https://www.mcs.anl.gov/petsc/documentation/installation.html
  - Remote install: see `docs/PETSC_Install.pdf`
- **FEniCS**: 2017.1.0 (use this exact version)
  - Local install: https://fenics.readthedocs.io/en/latest/installation.html
  - Remote install: follow `docs/fenics/fenics_install.pdf` (scripts inside the same folder may be adapted as needed)
  - Troubleshooting: [Compute Canada FEniCS wiki](https://docs.computecanada.ca/wiki/FEniCS)
  - Required Python packages: `meshio`, `lxml`, `mpi4py`
- **ParaView**: 5.4.0
- **Python**: 3.x
- **MATLAB**: 2017

---

## Local Workflow (tested on macOS)
Work from the top-level project directory and use separate terminals for concurrent steps.

### Step 1 – Preprocess (KLE/PCE and mesh generation)
```
$ cd data/klePceData
$ gfortran KLE_PCE_Data_commandLineArg.F90
$ ./a.out          # enter desired nDim (number of random variables)
$ cd ../meshData
```

- **3D domain**
  - Edit `preprocess3D.sh` to tune mesh density (`LC`) and partitions (`NP`).
  - Run `sh preprocess3D.sh` to create `.dat` mesh data.
  - Convert to FEniCS XML:
    ```
    $ cd ../../external/dolfin/src/
    $ python3 gmshData_fenicsXML_3D.py
    ```

- **2D domain**
  - Edit `preprocess2D.sh` (LC, NP).
  - Run `sh preprocess2D.sh` to create `.dat` meshes.
  - Convert to FEniCS XML:
    ```
    $ cd ../../external/dolfin/src/
    $ python3 gmshData_fenicsXML_2D.py
    ```

### Step 2 – Deterministic FEM assembly with FEniCS
```
$ fenicsproject start dolfin
$ cd external/dolfin/src/
```
- **3D domain**
  - `python poisson3D_stochasticDDM_onelevel.py` → one-level DD matrices in `data/Amats/`
  - `python poisson3D_stochasticDDM_twolevel.py` → two-level DD matrices in `data/Amats/`
- **2D domain**
  - `python poisson2D_stochasticDDM_onelevel.py` → one-level DD matrices
  - `python poisson2D_stochasticDDM_twolevel.py` → two-level DD matrices

### Step 3 – SSFEM solve (Fortran/MPI/PETSc)
```
$ cd src/2Dpoisson_*    # or 3Dpoisson_* / 3Delasticity_*
$ cp clusterMakefiles/makefile_mac makefile
$ make all
$ petscexec -np $NP ./a.out   # prints iteration counts and writes to data/vtkOutputs/
```

### Step 4 – Post-process (MATLAB / ParaView)
```
$ cd data/vtkOutputs
```
- For 3D Poisson: `matlabTerminal -r writeVtk`
- For 3D elasticity: `matlabTerminal -r writeVecVtk`
- Visualize: `paraview out_*.vtk` (2D VTK files are created directly).

### Special 2D comparison workflow (all-in-one)
1. Edit `preprocess2D.sh` for partitions and mesh density, then run it.
2. Convert to XML:
   ```
   $ cd ../../external/dolfin/src/
   $ python3 pointElements_xml.py
   ```
3. Start FEniCS (`fenicsproject start dolfin`) and generate DD matrices:
   ```
   $ cd external/dolfin/src/
   $ python stoDDM_poisson.py
   ```
4. Build and run either `src/onelevel` or `src/twolevel`:
   ```
   $ make all
   $ petscexec -np $NP ./a.out
   ```

---

## Remote Cluster Workflow (Cedar, Graham, Niagara, Beluga)

### Step 1 – Clone the repository (private Bitbucket repo)
```
$ git clone https://sudhipv@bitbucket.org/sudhipv/ssfem_fenics.git
$ cd ssfem_fenics
$ git fetch && git checkout <branch>   # defaults to master
```

### Step 2 – Preprocess and process (Gmsh / FEniCS / MPI / PETSc)
Run from `jobscripts/<machine-name>/` (choose `cedar`, `niagara`, `beluga`, or `graham`). Each script has inline comments describing required path edits.

1. `sh preprocess.sh`
   - Prompts for `nDim` (number of random variables) and writes KLE/PCE data to `data/klePceData/`.
   - Control mesh density via `lc` (e.g., `lc=0.1` → `lc=0.09`, keep the first `lc` fixed at `0.1`).
   - Select coarse grid type via `wb` (`1` = wire-basket, `0` = vertex).
   - Set partitions via `NP` (default 32). Mesh/DDM data land in `data/meshData/`.
2. `sh processFEniCS.sh` or `sh processFEniCS_MPI.sh`
   - Adjust job resources inside the script (e.g., `time=0-00:10` → `0-00:15`) and set the desired output directory.
   - Script names invoked by `processFEniCS_MPI.sh` depend on the physics (e.g., `run_fenics_poisson2D_stochasticDDM_parallel_twolevel.sh`, `run_fenics_elasticity3D_stochasticDDM_parallel_twolevel.sh`). Mail notifications trigger only on completion/failure; update addresses in `external/dolfin/src/jobscripts` as needed.
   - If MPI job initialization fails, inspect `/scratch/<user>/ssfem_fenics/data/slurm/` (e.g., `/scratch/sudhipv/sudhipv/ssfem_fenics/data/slurm/`) for details. A successful run prints the random variable settings, BC application, and `==========================Success============================`.

3. `sh processMPI.sh`
   - Controls Fortran compilation and MPI/PETSc execution.
   - Tunable parameters inside the script:
     1. `time=0-00:10` → modify run time (default 10 min).
     2. `tasks-per-node=32` (adjust to match allocation).
     3. `nodes=1` (increase if needed).
     4. `-np` flag (ensure `np = nodes * tasks-per-node`).
     5. `job-name`/`output-file` (e.g., `job-name=d2_lc05_p32`, where `d` = dimensions, `lc` = mesh density, `p` = partitions).
   - Outputs are written to `data/vtkOutputs/`.

### Step 3 – Post-process after remote runs
1. Copy `posprocess3D.sh` (located in `data/vtkOutputs/`) to your local machine, edit the file paths, and run it or follow the manual steps below:
   ```
   $ cd data/vtkOutputs
   $ scp <cluster>:~/scratch/<user>/ssfem_fenics/data/vtkOutputs/*.dat .
   ```
2. For 3D ParaView visualizations also copy:
   ```
   $ scp <cluster>:~/scratch/<user>/ssfem_fenics/data/meshData/points.dat .
   $ scp <cluster>:~/scratch/<user>/ssfem_fenics/data/meshData/tetrahedrons4c.dat .
   $ scp <cluster>:~/scratch/<user>/ssfem_fenics/data/klePceData/pcedata.dat .
   ```
3. Generate VTK:
   - Poisson 3D → `matlabTerminal -r writeVtk`
   - Elasticity 3D → `matlabTerminal -r writeVecVtk`
4. Visualize with `paraview out_*.vtk`.

### Job submission & monitoring
- Submit, monitor, and cancel jobs:
  ```
  $ qsub run_ddm.sh
  $ qstat -u <userID>
  $ qdel <jobID>
  ```
- Inspect outputs once jobs finish:
  ```
  $ more outputfile*
  $ more errorfile*
  ```
- Miscellaneous helpers:
  - `sacct` – view status for recent SLURM jobs.
  - `more slrum-<JobID>.out` – job logs.
  - `git fetch && git checkout graham`
  - `git checkout -f` – discard local edits (useful when pulling latest changes).

---

## Directory Layout
- `src/` – stochastic solvers
  - `3Delasticity_onelevel`, `3Delasticity_twolevel`
  - `3Dpoisson_onelevel`, `3Dpoisson_twolevel`
  - `2Dpoisson_onelevel`, `2Dpoisson_twolevel`
  - Legacy solvers: `onelevel`, `twolevel`
- `jobscripts/`
  - `preprocessKLE.sh`, `preprocessMESH.sh`, `preprocess.sh`
  - `processFEniCS.sh`, `processFEniCS_MPI.sh`
  - `processMPI.sh`
- `external/`
  - `puffin` (MATLAB FEM package)
  - `dolfin` (Python/C++ FEM package)
- `data/`
  - `klePceData`, `mallocData`, `meshData`, `vtkOutputs`
- `docs/`
  - Project documentation (`PETSC_Install.pdf`, `fenics_install.pdf`, etc.)

### Output visualization
- Deterministic case: `paraview out_deterministic.vtk`

### Cleaning
- `sh clean.sh` – remove build and output artifacts.

### Contact
- Sudhi Sharma P V – `sudhisharmapadillath@gmail.com`

## Reference

> **[Scalable Domain Decomposition Methods for Nonlinear and Time-Dependent Stochastic Systems](https://doi.org/10.22215/etd/2023-15817)**

**Authors:** Vasudevan, Padillath and Sharma, Sudhi  
**Institution:** Carleton University (2023)  
**DOI:** [10.22215/etd/2023-15817](https://doi.org/10.22215/etd/2023-15817)

<details>
<summary><b>Click to expand BibTeX citation</b></summary>

```bibtex
@phdthesis{vasudevan2023scalable,
  title={Scalable Domain Decomposition Methods for Nonlinear and Time-Dependent Stochastic Systems},
  author={Vasudevan, Padillath and Sharma, Sudhi},
  year={2023},
  school={Carleton University},
  doi={10.22215/etd/2023-15817}
}
\```
</details>
