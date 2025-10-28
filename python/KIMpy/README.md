# Python code for KIM

## kimpy
Handles preparation and runs of KIM.

## kim_data
Handles the output data of KIM.

## KIM WKB
Within the main file wkb.py, this python code calculates the dispersion relation within the KiLCA Integral Model (KIM).
It uses the cxroots module to find the roots to the dispersion equation. For benchmarking,
also the dispersion relation in finite Larmor radius expanded form of Horton (Horton 1999, Rev. of Mod. Phys.)
is used.

KIM comes in 3 different collision models:
- Krook: simple Krook collision model
- collisionless: Krook, but with collision frequencies set to zero
- Fokker_Planck: FP-type collision model with an Ornstein-Uhlenebeck collision
    model including an energy-preserving integral term. This collision model
    results in the susceptibility functions (I00 and I20) imported from susc_funcs.
    The needed Bessel functions are calculated in the Bessel_calculation file.
The modes are implemented with abstract base classes (ABC).
        
Data is handled with dictionaries. Different dicts are: 
- options: contains options for the run including possible collision models
- species_dat: contains all the species specific profiles like density, 
    temperature, as well as derived ones like the thermodynamic forces.
- general_dat: contains data that is not specific to a species, e.g. 
    radial electric field, radial grid,...
- equil_dat: contains data for the magnetic equilibrium of the cylinder that
    is calculated in calc_equilibrium
                
The output data is written to the hdf5 file format using the h5py module. 

CXroots needs: wrapt, numpydoc
