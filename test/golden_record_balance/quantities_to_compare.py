#!/usr/bin/env python3
"""
Quantities to compare in the QL-Balance golden record test.

==============================================================================
USER TODO: Define the quantities to compare!

Specify the HDF5 dataset paths and tolerances for quantities you want to
compare between the golden record and the test output.

Each QuantitySpec has:
    - path: HDF5 path to the dataset (e.g., "/profiles/Te" or "/output/torque")
    - rtol: Relative tolerance (default: 1e-8)
    - atol: Absolute tolerance (default: 1e-15)

To find available paths in your HDF5 file, you can use h5dump or Python:
    ```bash
    h5dump -H your_output.h5
    ```
    or
    ```python
    import h5py
    def print_structure(name, obj):
        print(name)
    with h5py.File("your_output.h5", "r") as f:
        f.visititems(print_structure)
    ```

==============================================================================
"""

from compare_hdf5 import QuantitySpec


# Default tolerances (can be overridden per-quantity)
DEFAULT_RTOL = 1e-8
DEFAULT_ATOL = 1e-15


# List of quantities to compare
# TODO: Add the HDF5 paths you want to compare!
QUANTITIES_TO_COMPARE: list[QuantitySpec] = [
    # Example entries (uncomment and modify as needed):
    #
    # QuantitySpec("/profiles/Te"),  # Electron temperature
    # QuantitySpec("/profiles/Ti"),  # Ion temperature
    # QuantitySpec("/profiles/n"),   # Density
    # QuantitySpec("/output/torque", rtol=1e-6),  # Torque with relaxed tolerance
    # QuantitySpec("/output/D11", rtol=1e-7, atol=1e-12),  # Diffusion coefficient
    #
    # Add your quantities below:
]


# Validate that we have at least one quantity to compare
if not QUANTITIES_TO_COMPARE:
    import warnings
    warnings.warn(
        "No quantities defined in QUANTITIES_TO_COMPARE. "
        "Please edit test/golden_record_balance/quantities_to_compare.py"
    )
