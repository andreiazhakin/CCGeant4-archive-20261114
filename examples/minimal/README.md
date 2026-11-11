# Cu, Xe and U transport through silicon

The application sends fully stripped Cu-63, Xe-124 and U-238 ions at normal incidence through 100 µm of natural silicon. The speed coordinate is (γ−1) × 931.494 MeV = 400 MeV. It uses the Dirac electron spectrum, charge-state transitions, multiple scattering and a 25 µm production range cut. The input seed is 6110101. The default run samples eight events per projectile on two threads. The calculation scores primary-ion kinetic-energy loss and the energy deposited in the silicon separately.

Build against the patched, multithreaded Geant4 11.4.0 installation described in the root README. From the repository root:

```sh
cmake -S examples/minimal/app -B build-example -DCMAKE_PREFIX_PATH=/path/to/patched-install
cmake --build build-example --config Release
python examples/minimal/run_example.py --executable build-example/dirac_subset_example --output example-results
python tests/audit_example_analyze.py example-results
```

On a Visual Studio generator the executable is `build-example/Release/dirac_subset_example.exe`. The runner accepts `--geant4-bin /path/to/patched-install/bin` for DLL lookup. Initialize the standard Geant4 datasets for this installation first.

The runner expands the eight example tables from `data/` into the output directory and writes per-event CSVs, Geant4 logs and a run summary. To use an existing copy of the tables, pass `--tables DIR` with the same directory layout. Python 3.8+ with its standard lzma module is sufficient for the example.

`reference_output.json` gives numerical results for the default fixed-seed run on the reference Windows x64 build. The output columns distinguish energy lost by the primary from energy deposited in the target. For example the eight Xe events give mean values of 172.523 and 152.814 MeV, respectively. Use the event-accounting check above when comparing builds. It checks that every requested event completed, the charge and delta-electron kinematics stay within physical bounds, and the event and primary energy balances are smaller than 1 eV. The deterministic cross-section and stopping-power check is documented in [the tests](../../tests/README.md).
