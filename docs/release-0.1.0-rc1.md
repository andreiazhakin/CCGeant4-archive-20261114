# CCGeant4 0.1.0-rc1

This release candidate packages the Geant4 11.4.0 patch, its 61 added or modified
source files, Dirac and elastic converters, the optional sparse elastic-table
generator, a runnable Cu/Xe/U example, numerical checks and benchmark summaries.
The production tables are distributed in `runtime-20261107`.

## Source and complete transport data

Get the source from the [0.1.0-rc1 release](https://github.com/andreiazhakin/CCGeant4/releases/tag/v0.1.0-rc1).
The release includes `CCGeant4-0.1.0-rc1-source.zip`.
The same files are available through Git:

```sh
git clone --branch v0.1.0-rc1 https://github.com/andreiazhakin/CCGeant4.git
cd CCGeant4
```

Download the five table archives below. They remain in
the versioned runtime release so that existing data references stay valid.

| Download | Contents | Archive size |
| --- | --- | ---: |
| [Ionization and stopping](https://github.com/andreiazhakin/CCGeant4/releases/download/runtime-20261107/CCGeant4-runtime-20261107-ionization.zip) | Full Dirac spectra, ATIMA-normalized stopping, GLOBAL transitions and low-energy ATIMA stopping | 1.230 GB |
| [Elastic Z = 26–40](https://github.com/andreiazhakin/CCGeant4/releases/download/runtime-20261107/CCGeant4-runtime-20261107-elastic-Z26-40.zip) | Silicon tables and the full elastic index | 1.202 GB |
| [Elastic Z = 41–55](https://github.com/andreiazhakin/CCGeant4/releases/download/runtime-20261107/CCGeant4-runtime-20261107-elastic-Z41-55.zip) | Silicon tables | 1.211 GB |
| [Elastic Z = 56–70](https://github.com/andreiazhakin/CCGeant4/releases/download/runtime-20261107/CCGeant4-runtime-20261107-elastic-Z56-70.zip) | Silicon tables | 1.216 GB |
| [Elastic Z = 71–82](https://github.com/andreiazhakin/CCGeant4/releases/download/runtime-20261107/CCGeant4-runtime-20261107-elastic-Z71-82.zip) | Silicon tables | 0.976 GB |

Alternatively, with GitHub CLI:

```sh
gh release download runtime-20261107 --repo andreiazhakin/CCGeant4 --pattern 'CCGeant4-runtime-20261107-*.zip' --dir archives
python data/install_runtime.py --archives archives --output DATA
```

Extract into a new directory. To read the ZIP contents without extracting them, use `--verify-only`. The [README](../README.md#download-checksums) gives the archive checksums.

| Unpacked component | Number of files | Size |
| --- | ---: | ---: |
| Dirac electron spectra | 64 | 1.588 GB |
| ATIMA-normalized stopping tables | 64 | 0.203 GB |
| GLOBAL charge-changing cross sections | 1 | 3.144 GB |
| Low-energy ATIMA stopping | 1 | 0.052 GB |
| Silicon elastic tables and index | 58 | 6.251 GB |
| **Total** | **188** | **11.238 GB** |

Sizes use decimal GB. The archives total 5.835 GB. The dense elastic angular
and energy grids preserve the interference structure during interpolation.
The runtime loads the tables needed by the simulation and shares them between
threads.

Follow the [build instructions](../README.md#build-against-geant4-1140), then set
the data paths as shown in [physics configuration](physics-configuration.md).
The [small example](../examples/minimal/README.md) includes its own Cu/Xe/U
subset for a quick installation check. It does not require the full bank.

## Checks and scope

The release preparation verified all 188 runtime files, application of the
patch to clean Geant4 11.4.0, all 61 resulting source files, the independent
partial-moment tests, the eight bundled example inputs, and loading a sparse Xe–Si table. The [check record](../validation/recorded_checks/release_0.1.0_rc1.json)
identifies these checks. No new transport simulation was needed for this
packaging update.

The [validation directory](../validation/README.md) contains benchmark summaries
and software checks. It is not a complete archive of every calculation in the
evolving manuscript. Upstream atomic solvers and potential databases are
obtained separately. The full offline table-generation workflow is not included;
see the [data provenance](dirac-elastic-data.md).

Independently authored code and documentation use the MIT license. Geant4 code
and adaptations retain the Geant4 license and notices. Calculated tables are
distributed with source credits described in [LICENSING.md](../LICENSING.md).
Use the release tag or full commit identifier when reporting calculations.
