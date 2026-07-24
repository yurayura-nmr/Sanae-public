# SANAE
**Simulation Ab Initio for NMR, Anisotropy, and Exchange**

SANAE is a C++ toolkit for biomolecular biophysics simulation and experimental planning.
It covers NMR spectroscopy, fluorescence anisotropy, ITC calorimetry, binding analysis,
kinetics, and more. 
All driven by a very simple text input file.

> **Note**
> I started this project to explore how magnetization changes in NMR experiments. 
> The analytical equations while often elegant rely on approximations and it is tough to keep track of the regime in which what equation remains valid to what extent.
> In SANAE, NMR FID simulations use direct numerical integration of the Bloch-McConnell equations.
> In other words, exchange, relaxation, and precession are handled purely through numerical time-domain integration. No analytical solutions are assumed.

---

## Installation

### Linux / macOS
```bash
cd build
cmake .. && make -j
cd ../bin/
./Sanae <input file>
```

### Windows

SANAE supports two build workflows on Windows: **Visual Studio (MSVC)** and **CMake with MinGW (GCC)**. 
Use the one that fits your environment.

---

#### Option 1: Visual Studio (MSVC)

This step is necessary because my own absolute path for where Eigen is installed will differ from yours.

1. Open the repository folder and double-click the solution file (`.sln`), or open Visual Studio and select **Open a project or solution**.
2. Ensure **Eigen** is available on the include path:
   - Right-click the project → **Properties**
   - Go to **C/C++ → General → Additional Include Directories**
   - Add the path to Eigen (included in Sanae/include)
3. Build the solution (F7 or **Build → Build Solution**).
4. The executable will be placed in the output folder (e.g., `x64\Release\` or `Debug\`). 
5. You can now run Sanae.exe on an input `sanae.in`.

---

#### Option 2: CMake + MinGW (GCC)

This may be simpler if you don't want to install Visual Studio.

1. **Install the build tools** (using [Chocolatey](https://chocolatey.org/) is recommended):
   ```bash
   choco install mingw make cmake
   ```
   Alternatively, install [MinGW-w64](https://www.mingw-w64.org/) and [CMake](https://cmake.org/download/) manually.

2. **Configure and build** using CMake GUI:
   - Open CMake GUI
   - Set source directory to the repository root, build directory to `build/`
   - Click **Configure**, select **"MinGW Makefiles"** as the generator
   - Click **Generate**
   - Close CMake GUI, then run the following in the `build/` directory:
   ```bash
   make -j
   ```
   The executable will be created in the `bin/` subdirectory (or your build output folder).

3. You can now run Sanae.exe on an input `sanae.in`.

---

## Usage

SANAE is driven by a simple input file. 
You can find examples of input files for different applications in Sanae/examples.

```bash
./Sanae                      # reads sanae.in (default) if a sanae.in is found in the current directory
./Sanae my_r2_dispersion.in  # reads a named input file
```

### Example `sanae.in`
```ini
# Example: simple 2-state NMR exchange between states A and B, each having their chemical shift, R2, and kinetic rates between them.
run_nmr2state   = 1
k_ab            = 573.0
k_ba            = 274.0
cs_free         = 0.0
cs_bound        = 1.0
r2_free         = 10.0
r2_bound        = 15.0
mag_free        = 0.55
mag_bound       = 0.34
nucleus         = N
dt              = 0.00001
steps           = 20000
print_every     = 10
output_file     = fid.txt
```

Lines beginning with `#` are comments. Units are the most familiar ones to NMR spectroscopists, e.g. rates in 1/s etc.
Although not recommended, multiple run modes can be active in the same file. 
If any necessary parameters are missing from the file, SANAE will fall back to defaults. 
The `[DEFAULT]` diagnostic printed at runtime shows which values were used.
Check the output to catch this.

### Visualising NMR output
```bash
python3 sanaeplot.py fid.txt
```
To plot a simulated FID and NMR spectrum, this Python script (found in Sanae/py) reads simulation parameters (`dt`, `steps`, `print_every`) directly from the output file header created by SANAE.


---

## Run Modes

### NMR Exchange Simulations


### CPMG Relaxation Dispersion


### Binding and Kinetics


### Relaxation Parameter Estimation


### ITC Calorimetry


### Fluorescence Anisotropy


### Hydrogen Exchange


### Other


---

## Applications

- NMR lineshape simulation across exchange regimes (slow / intermediate / fast)
- CPMG relaxation dispersion prediction and temperature planning
- Activation energy (dH‡) estimation from temperature-dependent kex
- Binding curve and titration simulation (1-, 2-, 3-site; cooperative)
- ITC thermogram simulation for experimental design and model comparison
- Fluorescence anisotropy prediction
- CLEANEX-PM hydrogen exchange rate planning
- Relaxation parameter estimation (τc, J(ω), R2) from molecular weight
- Noise-added FID simulation for SNR planning

---

## Dependencies

- [Eigen](https://eigen.tuxfamily.org) — linear algebra (header-only)
- Python: `numpy`, `scipy`, `matplotlib` (for `sanaeplot.py`)