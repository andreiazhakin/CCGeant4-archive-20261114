# Numerical checks

Build the probes against the patched Geant4 installation (CMake, C++17 and multithreading):

```sh
cmake -S tests -B build-tests -DCMAKE_PREFIX_PATH=/path/to/patched-install
cmake --build build-tests --config Release
```

## Deterministic Xe-on-Si check

From the repository root, expand the bundled tables, evaluate them using the Geant4 C++ reader, and compare the result to the independent Python evaluation. On Windows use the `Release` subdirectory and the `.exe` suffix for Visual Studio builds and add the patched installation's DLL directory to PATH.

```sh
python examples/minimal/prepare_tables.py example-input
build-tests/table_kernel_probe example-input/delta/Z054.csv example-input/stopping/Z054.csv 54 xe-probe.csv
python tests/check_xe_kernel.py xe-probe.csv
```

The fixture covers 11 energies, four numbers of bound electrons and four energy thresholds. For each of the 176 combinations it checks the zeroth, first and second moments of the emitted-electron spectrum, the charge-specific stopping, equilibrium stopping and soft coefficient. Units are cm² MeV^k per target electron for the kth moment, and MeV cm² per target atom for stopping and its charge coefficient. Relative tolerance is 3 × 10⁻¹⁰, with an exact check for zero reference values. The probe also exercises delta sampling and checks the kinematic endpoints. The reference grid was evaluated with the converter's Python routines on the original Xe binary packs, independently of the C++ CSV reader. A fresh Windows x64 build gave a maximum relative difference of 5.24 × 10⁻¹⁴.

`validate_table_kernel.py --tables PACK_ROOT --probe xe-probe.csv --out report.json` repeats the comparison directly against local binary packs, including a statistical test of the sampled mean. This optional route needs NumPy and the pack directory containing `delta/` and `stopping/`.

## Event transport and additional probes

The [Cu/Xe/U example](../examples/minimal/README.md) includes expected numerical output and a portable event and primary energy-balance check. The other C++ targets exercise charge-state table coverage, isotope scaling, electron angles, trackless calls, and range cuts. The independent converter checks for zero transition rates and partial spectral moments are in `converters/dirac/`.
