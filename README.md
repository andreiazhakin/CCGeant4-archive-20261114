# CCGeant4

CCGeant4 models heavy-ion charge exchange and energy-loss fluctuations in
Geant4. It couples discrete capture and loss events from GLOBAL to
charge-dependent stopping normalized to ATIMA and explicit delta-electron
production from frozen-shell Dirac spectra. Applications can select
multiple scattering, native ion single scattering, or tabulated nuclear
elastic scattering with recoil transport. The supplied elastic dataset
covers 57 projectile elements on silicon.

## Complete release

[Version 0.1.0-rc1](https://github.com/andreiazhakin/CCGeant4/releases/tag/v0.1.0-rc1)
is the release candidate prepared for the software article. The
[release guide](docs/release-0.1.0-rc1.md) links the versioned source and all five
runtime-table archives, and gives the installation commands.
The five archives contain 188 files (11.24 GB unpacked, 5.84 GB compressed).
The separate [Legacy stopping table](#legacy-charge-scaled-stopping) is also
included in the checkout.

The [download checksums](#download-checksums) identify the distributed data files.

The separate [article-data archive](https://github.com/andreiazhakin/CCGeant4/releases/tag/article-data-20261109)
contains the recorded simulation events, digitized references, analysis scripts
and figures used in the manuscript. Its [reproduction guide](docs/article-data-20261109.md)
explains how to verify the files and rebuild the plots without rerunning Geant4.

## Build against Geant4 11.4.0

Start from a clean [Geant4 11.4.0 source tree](https://gitlab.cern.ch/geant4/geant4/-/tree/v11.4.0).
Use a C++17 compiler, CMake and a multithreaded Geant4 build. From the CCGeant4 repository root:

```sh
git -C /path/to/geant4-v11.4.0 apply --check /absolute/path/to/CCGeant4/patches/ccgeant4-geant4-11.4.0.patch
git -C /path/to/geant4-v11.4.0 apply /absolute/path/to/CCGeant4/patches/ccgeant4-geant4-11.4.0.patch
cmake -S /path/to/geant4-v11.4.0 -B build-g4 -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DBUILD_SHARED_LIBS=ON -DGEANT4_BUILD_MULTITHREADED=ON -DGEANT4_INSTALL_DATA=ON -DCMAKE_INSTALL_PREFIX=/path/to/patched-install
cmake --build build-g4 --config Release --parallel
cmake --install build-g4 --config Release
```

Initialize the Geant4 environment for that installation. The patch is the
installation method. The same 8 modified and 53 added files can be inspected
individually in [source_changes/](source_changes/).

To verify the source patch locally without building the toolkit, run
`python tools/verify_patch.py --base /path/to/geant4-v11.4.0.zip`. It checks
that the patch applies and the 61 resulting files match the distributed sources in a temporary directory.

## Run the silicon example

The checkout includes compressed input tables for Cu, Xe and U. The runner
expands them automatically. A standard Python installation with
lzma support is sufficient to run it.

```sh
cmake -S examples/minimal/app -B build-example -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/patched-install
cmake --build build-example --config Release
python examples/minimal/run_example.py --executable build-example/dirac_subset_example --output example-results
python tests/audit_example_analyze.py example-results
```

With Visual Studio, pass `build-example/Release/dirac_subset_example.exe`.
The [example](examples/minimal/README.md) gives the beam and target settings,
expected energy loss and energy deposit, and the event-balance check. It runs
with the supplied subset and the standard Geant4 datasets.

## Legacy charge-scaled stopping

In a [modular physics list](docs/physics-configuration.md), select Legacy
before Geant4 initialization.

```cpp
physics->ReplacePhysics(new G4EmStandardPhysics_option4_CurrentCharge(0));
```

This constructor reads its inputs from the environment. [`data/legacy_q2.csv`](data/legacy_q2.csv)
gives the ATIMA equilibrium mass stopping power and the GLOBAL equilibrium
mean square charge for each projectile, target and energy. For an ion of
charge $q$ in units of the elementary charge, Legacy multiplies the ATIMA
value by $q^2/\langle q^2\rangle_{\rm eq}$. The other two inputs are
`global.csv` for charge exchange and `atima_low.csv` for low-energy
equilibrium stopping.

From the repository root, install the [ionization archive](data/README.md#installation)
in a new directory, then set the variables below before launching the application.

```sh
python data/install_runtime.py --archives /path/to/archives --ionization-only --output runtime-tables
```

**PowerShell**

```powershell
$env:G4CC_GLOBAL_TABLE = (Resolve-Path 'runtime-tables/global.csv').Path
$env:G4CC_ATIMA_Q_DEDX_TABLE = (Resolve-Path 'data/legacy_q2.csv').Path
$env:G4CC_ATIMA_EQUILIBRIUM_DEDX_TABLE = (Resolve-Path 'runtime-tables/atima_low.csv').Path
$env:G4CC_BINDING_OVERRIDE_MEV = $null
$env:G4CC_BINDING_OVERRIDE_GEV = '1000'
```

**Bash**

```sh
export G4CC_GLOBAL_TABLE="$PWD/runtime-tables/global.csv"
export G4CC_ATIMA_Q_DEDX_TABLE="$PWD/data/legacy_q2.csv"
export G4CC_ATIMA_EQUILIBRIUM_DEDX_TABLE="$PWD/runtime-tables/atima_low.csv"
unset G4CC_BINDING_OVERRIDE_MEV
export G4CC_BINDING_OVERRIDE_GEV=1000
```

The 1000 GeV override reproduces the screening assumption used in the Legacy
benchmarks. It keeps the charge entering the normalization of the delta-electron
spectrum close to the instantaneous ionic charge over the electron energies
studied. Without an override, Legacy estimates the binding energies of individual
shells. A positive value of `G4CC_BINDING_OVERRIDE_MEV` takes precedence over
`G4CC_BINDING_OVERRIDE_GEV`.

The Weick comparisons used the default values for the following settings.

| Environment variable | Value | Role |
| --- | --- | --- |
| `G4CC_USE_ATIMA_Q_DEDX` | `1` | Use the charge-scaled stopping table. |
| `G4CC_USE_ATIMA_EQUILIBRIUM_DEDX` | `1` | Use tabulated equilibrium stopping at low energies. |
| `G4CC_USE_MOTT_DELTA`, `G4CC_MOTT_MIN_ENERGY_MEV_U` | `1`, `2` | Use the Mott correction for explicit delta electrons from 2 MeV/u. |
| `G4CC_DISABLE_CONTINUOUS_LS_FOR_MOTT_DELTA` | `1` | Suppress the separate continuous Lindhard–Sørensen correction when the Mott correction is active. |
| `G4CC_ENABLE_SMOOTH_GLOBAL_FALLBACK` | `1` | Use smooth transitions at the ends of the table range. |
| `G4CC_GLOBAL_WEIGHT_LOW0_MEV_U`, `G4CC_GLOBAL_WEIGHT_LOW1_MEV_U` | `30`, `40` | Lower transition limits with equilibrium stopping enabled. |
| `G4CC_GLOBAL_WEIGHT_HIGH0_MEV_U`, `G4CC_GLOBAL_WEIGHT_HIGH1_MEV_U` | `1500`, `2000` | Upper transition limits. |
| `G4CC_FLUCTUATION_MODE`, `G4CC_FLUCTUATION_MIX_POWER` | `mix`, `1` | Select Universal or ion fluctuations with the transition weights. |

To try this configuration with the built silicon example, run the executable
with `--legacy 1 --Z 54 --A 124 --N 2 --energy 400 --events 8 --out legacy-xe.csv`.

## Using the model

- [Physics configuration](docs/physics-configuration.md) gives a minimal C++
  registration, data paths, thresholds and the elastic-scattering option.
- [Production tables](data/README.md) identifies the complete ionization and
  silicon elastic banks. The [runtime release](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107)
  provides the matching source and all five table archives.
- [Ionization converter](converters/dirac/README.md) and
  [elastic converter](converters/elastic/README.md) describe rebuilding
  the runtime files from calculated spectra and scattering arrays.
  The elastic tools also include an optional [sparse-table generator](converters/elastic/README.md#generate-a-sparse-coulomb-table-for-one-isotope-pair)
  for a geometric estimate of large-angle Coulomb suppression.
- [Numerical tests](tests/README.md) and [transport benchmarks](validation/README.md)
  provide runnable checks and simulation results for the Weick and Xe–Si comparisons.
- [License and provenance](LICENSING.md) covers the MIT project code, Geant4
  modifications and the numerical inputs. Cite the upstream methods when
  publishing calculations with the corresponding components.


## Download checksums

| Runtime download | SHA-256 |
| --- | --- |
| `CCGeant4-runtime-20261107-ionization.zip` | `0ed618a78ffe98f06fd6125a2692b1ba0b891432534753d86bda04fb0e56ba1d` |
| `CCGeant4-runtime-20261107-elastic-Z26-40.zip` | `38b5f2f859ba06e440ae4a80f6c622823c0d1f52a743d646dfdb4ddbde6b4b9b` |
| `CCGeant4-runtime-20261107-elastic-Z41-55.zip` | `e803dff4502bfa9ccfb5dced17bc881256f83d99cde08a0770ef270cf98aa358` |
| `CCGeant4-runtime-20261107-elastic-Z56-70.zip` | `eb5668bb0a25b4d800067773440145f82b769839d13c88689a39868fa6c9b6a3` |
| `CCGeant4-runtime-20261107-elastic-Z71-82.zip` | `3f919ab7d576508c76e25331dfe41dcc96ba57f8fd76916a319087affb3fa0c3` |
| `legacy_q2.csv` | `17ea430591081e0aade1f0c142937d1567c50c09239db39f5305c6c3521bda8d` |

| Source and article archives | SHA-256 |
| --- | --- |
| `CCGeant4-0.1.0-rc1-source.zip` | `569b07ac723cf043bc03de7fb2d1c44f046928ae34097beb823afbbe163faa12` |
| `CCGeant4-article-data-20261109.zip` | `1049c2c40993b56f68ac5893618574588173987800a036df028fc59075a2b050` |
| `CCGeant4_arxiv_source_20261109.zip` | `0cd9be437f1a15d5ef1fd8b1c1fc8c259e254996171626bcfa2c865f9158f0c0` |

This product includes software developed by Members of the Geant4 Collaboration (http://cern.ch/geant4).
