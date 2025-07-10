# KIM GUI

A graphical user interface for KIM (KiLCA Integral Model) that allows you to:
- Configure KIM namelist parameters through an intuitive interface
- Run KIM simulations with real-time progress monitoring
- Plot and visualize results
- Generate constant plasma profiles

## Features

### Configuration Tab
- **KIM Configuration**: Basic run parameters (input/output paths, debug levels, collision models)
- **KIM Setup**: Physics parameters (magnetic field, mode numbers, frequency settings)
- **KIM Grid**: Grid resolution and spatial parameters
- **Load/Save**: Import and export configuration files in namelist format

### Run & Results Tab
- **Profile Generation**: Create constant plasma profiles for testing
- **Run Controls**: Start/stop KIM simulations with progress tracking
- **Output Monitor**: Real-time display of KIM output messages
- **Results Plotting**: Automatic visualization of electric field components and potential
- **File Management**: Direct access to output directories

## Requirements

- Python 3.6+
- tkinter (usually included with Python)
- numpy
- matplotlib
- f90nml (for namelist file handling)
- h5py (for HDF5 support)

## Installation

1. Ensure all dependencies are installed:
   ```bash
   pip install numpy matplotlib f90nml h5py
   ```

2. Make sure the KIM executable is built:
   ```bash
   cd /path/to/KAMEL
   make KIM
   ```

## Usage

### Launch the GUI

#### From the KIMgui directory:
```bash
python launch_gui.py
```

or

```bash
python kim_gui.py
```

#### From anywhere (after setup):
```bash
python kim_gui_launcher.py
```

or (if you've added the script to your PATH):
```bash
kim-gui
```

#### Global Installation (Optional):
To run the GUI from anywhere, add the KIMgui directory to your PATH or create a symlink:

```bash
# Option 1: Add to PATH (add this to your ~/.bashrc or ~/.zshrc)
export PATH="/path/to/KAMEL/python/KIMgui:$PATH"

# Option 2: Create symlink in a directory already in PATH
ln -s /path/to/KAMEL/python/KIMgui/kim-gui /usr/local/bin/kim-gui
```

### Basic Workflow

1. **Set Run Path**: Choose where to run KIM and store results
2. **Configure Parameters**: Adjust physics and grid parameters in the Configuration tab
3. **Generate Profiles**: Create constant plasma profiles (or use existing ones)
4. **Run Simulation**: Click "Run KIM" to start the simulation
5. **View Results**: Use "Plot Results" to visualize the output fields

### Configuration Parameters

#### KIM Configuration
- **Profile Location**: Directory containing plasma profiles (n.dat, Te.dat, Ti.dat, etc.)
- **Output Path**: Where to save KIM results
- **HDF5 Input/Output**: Enable HDF5 format for input/output files
- **Debug/Status Levels**: Control verbosity of output messages
- **Collision Model**: Choose between 'Krook' and 'FokkerPlanck' models
- **Type of Run**: Select 'electrostatic' or 'electromagnetic' calculation

#### KIM Setup
- **Toroidal Magnetic Field**: Magnetic field strength (typically negative)
- **Major Radius**: Tokamak major radius in cm
- **Mode Numbers**: Poloidal (m) and toroidal (n) mode numbers
- **Frequency**: Mode frequency (0 for eigenmode calculation)
- **Field Type**: Specification of external magnetic perturbation

#### KIM Grid
- **Plasma Radius**: Radial extent of the plasma in cm
- **Grid Dimensions**: Number of points in k-space and configuration space
- **Resonance Parameters**: Control resolution near resonant surfaces
- **Integration Parameters**: Gaussian quadrature settings

## Output Files

KIM generates several output files in the `out/m{m_mode}_n{n_mode}/fields/` directory:
- **phi_sol.dat**: Electric potential solution
- **E_perp.dat**: Perpendicular electric field
- **E_perp_psi.dat**: Perpendicular electric field in flux coordinates
- **E_perp_MA.dat**: Magnetic component of perpendicular electric field

## Troubleshooting

### Common Issues

1. **"KIM executable not found"**
   - Ensure KIM is built: `make KIM` in the KAMEL directory
   - Check that the CODE environment variable is set correctly

2. **"Profiles not found"**
   - Use "Generate Constant Profiles" to create test profiles
   - Ensure profile files exist in the specified profile location

3. **"Cannot import f90nml"**
   - Install the required package: `pip install f90nml`

4. **GUI doesn't start**
   - Ensure tkinter is available (usually included with Python)
   - Try running with: `python -m tkinter` to test tkinter installation

### Tips

- Start with the default configuration and modify parameters gradually
- Use constant profiles for initial testing before using real experimental data
- Check the output text for error messages during runs
- The GUI automatically saves the configuration to the run directory before execution

## Development

The GUI is built using tkinter and integrates with the existing KIMpy interface. Key components:

- **kim_gui.py**: Main GUI implementation
- **launch_gui.py**: Simple launcher script
- **KIMpy integration**: Uses the existing KIMpy class for KIM execution

For modifications or extensions, refer to the KIMpy documentation and the KIM namelist structure.