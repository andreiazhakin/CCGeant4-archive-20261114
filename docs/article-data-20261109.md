# Article data and figure reproduction

The [article-data-20261109 release](https://github.com/andreiazhakin/CCGeant4/releases/tag/article-data-20261109)
archives the numerical inputs, recorded simulations and analysis used in the
software manuscript. It accompanies source release
[v0.1.0-rc1](https://github.com/andreiazhakin/CCGeant4/releases/tag/v0.1.0-rc1)
and the complete transport tables in
[runtime-20261107](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107).
The archive manifest records the exact source commit.

Download the data archive from that release. Rebuilding the plots uses the
recorded events and does not require Geant4 or the runtime tables. The
[release index](../validation/article-data-20261109.json) lists both editions.

## Contents

The archive includes:

- GLOBAL charge-evolution and MOCADI energy-loss comparisons, with the saved
  transport results and numerical references.
- Weick experimental digitizations, published-model digitizations, simulation
  events, charge selections, spectra, widths, skewness and source records.
- Xe–Si detector comparisons, including processed experimental spectra,
  numerical peak fits and matched Dirac, Legacy and emopt4 calculations.
- Silicon energy scans, the pure-Universal control, production-cut comparisons
  and the event records used for the 99% exit-angle selection in Section 7.
- Bi loss and deposition at 500 and 1000 MeV/u with standard electron transport and
  MicroElec. Figure 13 and Table 6 compare CCGeant4 Dirac with Geant4 emopt4.
  The archive also contains the Legacy comparison, Xe–Si single scattering and the stopping-power check.
- The Appendix B interpolation check for Cu, Xe and U, including its numerical results
  and reproduction instructions in `workspace/output/dirac_loo_review/`.
- The recorded performance measurements, run configurations and analysis
  scripts. Timing records describe the measured machine and build.

`MANIFEST.json` lists each payload file, its purpose and size. `plot_commands.json` maps manuscript figure labels to scripts,
arguments, inputs and outputs. The `workspace/` tree preserves paths within the
original calculation workspace. `external/` holds the numerical records that
were stored elsewhere, including the companion xenon analysis. Journal PDFs,
full experimental detector DAQ and compiled Geant4 installations are not part
of this archive.

## Verify and reproduce

Extract the archive, then run these commands from its top-level directory:

```sh
python -m pip install -r requirements.txt
python archive_tools.py prepare --output ../CCGeant4-article-replay
cd ../CCGeant4-article-replay
python archive_tools.py plots
```

The original records are kept unchanged. `prepare` creates a separate copy and
resolves the archive paths in that copy. Choose a new destination outside
the extracted archive. The plot commands analyze recorded events or use the
saved numerical summaries; they do not launch transport simulations. Generated
files remain in the replay copy, with command logs and `replay_report.json`.

To reproduce a particular group, use `python archive_tools.py plots --only ID`.
The principal command groups are:

| Command ID | Manuscript comparison |
| --- | --- |
| `charge-evolution` | Bi charge fractions against GLOBAL |
| `legacy-mocadi` | MOCADI comparison and Legacy loss/deposition spectra |
| `weick-variances-spectra` | U and Bi variances and Bi spectra |
| `weick-skewness` | Bi asymmetry |
| `weick-normalized` | Bi variances normalized to collisional straggling |
| `single-scattering` | Native and tabulated Xe–Si single scattering |
| `xenon-detectors` | Xe signals in 300 and 675 µm silicon |
| `section7-angle99` | Energy scans, Universal control and production-cut plots |
| `bismuth-secondary-statistics` | Recompute Bi loss/deposition distributions and Table 6 widths |
| `bismuth-secondary-figure` | Draw Figure 13 for emopt4 and Dirac |

Run `bismuth-secondary-statistics` before `bismuth-secondary-figure` when
recomputing that comparison from events. The complete command map also identifies
optional analyses and recorded table results. Performance timings are read from
the saved runs; the plotting workflow does not measure them again.

## Provenance and citation

Digitization records identify the source publication and figure or table.
Run configurations identify the simulated case and selection. Cite the
[Weick source publications](../validation/README.md#transport-benchmarks) and
the companion xenon experimental work when using those measurements. A
digitization is a numerical reading of a published result, not a new
experimental measurement.

The software license and upstream credits remain those documented in
[LICENSING.md](../LICENSING.md) and the archive's `software-terms/` directory.
For a calculation based on this release, record the source commit, runtime-table
release and article-data release together. This archive preserves the article's
analysis inputs; it does not include upstream atomic solvers or the full offline
table-generation workflow.
