"""Audit A01: explicit GLOBAL zeros, all interpolation paths and equilibrium.

Tests edit copies in --report-dir only; production tables are never rewritten.
The optional baseline source directory must be the unmodified pre-fix snapshot.
"""
import argparse
import copy
import csv
import importlib.util
import json
from pathlib import Path
import sys
import warnings

import numpy as np

from cc_binary import read_pack
from cc_inputs import (read_global, read_atima, equilibrium_from_rates,
                       interpolate_global_rates, TRANSITIONS)
from cc_physics import bracket
from convert import rates_at, normalization_target


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def reject(call, label):
    try:
        call()
    except ValueError:
        return
    raise RuntimeError(f'{label} was not rejected')


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def reference(g, energy):
    """Independent scalar version of the C++ rule, in converter float64 units."""
    e, lo, hi, weights = bracket(g.energies_MeV_u, energy)
    result = np.empty((len(e), len(g.transitions)))
    for i, (low, high, w) in enumerate(zip(lo, hi, weights)):
        for j in range(len(g.transitions)):
            a, b = g.sigma_cm2[0, 0, low, j], g.sigma_cm2[0, 0, high, j]
            result[i, j] = (np.exp(np.log(a)+w*(np.log(b)-np.log(a)))
                            if a > 0 and b > 0 else a+w*(b-a))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--global-csv', type=Path, required=True)
    parser.add_argument('--delta-pack', type=Path, required=True)
    parser.add_argument('--atima-splines', type=Path, required=True)
    parser.add_argument('--baseline-dir', type=Path, required=True)
    parser.add_argument('--report-dir', type=Path, required=True)
    args = parser.parse_args()
    out = args.report_dir
    inputs = out/'inputs'
    inputs.mkdir(parents=True, exist_ok=True)
    with args.global_csv.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        original = [r for r in reader if int(r['Zp']) == 54 and int(r['Zt']) == 14]
    require(len(original) == 4134, 'Expected real Xe/Si: 106 transitions * 39 energies')
    energies = sorted({float(r['E_MeV_u']) for r in original})

    def write(name, rows):
        path = inputs/(name+'.csv')
        with path.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        return path

    paths = {'positive': write('positive', original)}
    for name, zero_energies in [('zero_zero', set(energies)),
                                ('zero_positive', {energies[0]}),
                                ('positive_zero', {energies[1]}),
                                ('both_zero_first_cell', set(energies[:2]))]:
        rows = copy.deepcopy(original)
        for row in rows:
            if (row['Qe_in'] == '0' and row['Qe_out'] == '2'
                    and float(row['E_MeV_u']) in zero_energies):
                row['sigma_cm2'] = '0'
        paths[name] = write(name, rows)
    delta = read_pack(args.delta_pack)[1]
    atima = read_atima(args.atima_splines, [54], [14])
    queries = np.unique(np.r_[energies, np.sqrt(np.array(energies[:-1])*energies[1:]),
                               35., 100., 400., 1500.])
    summaries = []
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with np.errstate(divide='raise', invalid='raise', over='raise'):
            for name, path in paths.items():
                g = read_global(path, 54, 54)
                expected = reference(g, queries)
                a = g.interpolate_rates(queries)[0, 0]
                b = rates_at(g, 0, queries)[0]
                require(np.array_equal(a, expected), f'{name}: GlobalData differs from reference')
                require(np.array_equal(b, expected), f'{name}: rates_at differs from reference')
                require(np.all(np.isfinite(a)) and np.all(a >= 0), f'{name}: invalid interpolation')
                for energy in (energies[0], energies[1], energies[-1], 35., 400.):
                    require(np.array_equal(g.interpolate_rates(energy)[0, 0],
                                           reference(g, [energy])[0]), f'{name}: scalar interpolation')
                p, diagnostic = equilibrium_from_rates(a, g.transitions, return_diagnostics=True)
                normalized, norm_diagnostic = normalization_target(g, 0, 0, atima, delta, queries)
                require(np.array_equal(p, normalized['p_eq']), f'{name}: normalization uses different rates')
                require(np.all(np.isfinite(normalized['B_MeV_cm2'])), f'{name}: invalid normalization')
                residual = float(np.max(np.abs(np.sum(p*normalized['S_full_MeV_cm2'], axis=-1)
                                                /normalized['atima_MeV_cm2']-1)))
                require(residual < 1e-12, f'{name}: equilibrium mean identity')
                summaries.append(dict(name=name, passed=True, energies=len(queries),
                                      compared_rates=int(a.size), mean_identity_error=residual,
                                      stationary_residual=diagnostic['max_stationary_residual']))
            old_inputs = load_module('audit_old_cc_inputs', args.baseline_dir/'cc_inputs.py')
            old_convert = load_module('audit_old_convert', args.baseline_dir/'convert.py')
            g = read_global(paths['positive'], 54, 54)
            require(np.array_equal(g.interpolate_rates(queries),
                                   old_inputs.GlobalData.interpolate_rates(g, queries)),
                    'Positive GlobalData results changed')
            require(np.array_equal(rates_at(g, 0, queries), old_convert.rates_at(g, 0, queries)),
                    'Positive rates_at results changed')
            new_norm, _ = normalization_target(g, 0, 0, atima, delta, queries)
            old_norm, _ = old_convert.normalization_target(g, 0, 0, atima, delta, queries)
            for key in new_norm:
                require(np.array_equal(new_norm[key], old_norm[key]), f'Positive normalization changed: {key}')
            for value in (-1., np.nan, np.inf, -np.inf):
                reject(lambda: interpolate_global_rates([value], [1.], [.5]), f'lower={value}')
                reject(lambda: interpolate_global_rates([1.], [value], [.5]), f'upper={value}')
            for value in (-.01, 1.01, np.nan, np.inf):
                reject(lambda: interpolate_global_rates([0.], [1.], [value]), f'weight={value}')
            for name, value in [('negative', '-1e-20'), ('nan', 'nan'), ('inf', 'inf')]:
                rows = copy.deepcopy(original)
                rows[0]['sigma_cm2'] = value
                path = write(name, rows)
                reject(lambda: read_global(path, 54, 54), name)
            path = write('missing_transition', original[1:])
            reject(lambda: read_global(path, 54, 54), 'missing transition')
            rows = copy.deepcopy(original)
            for row in rows:
                row['sigma_cm2'] = '0'
            all_zero = read_global(write('all_zero_generator', rows), 54, 54)
            require(np.all(all_zero.interpolate_rates(queries) == 0), 'all-zero interpolation')
            reject(all_zero.equilibrium, 'zero generator')
            reject(lambda: normalization_target(all_zero, 0, 0, atima, delta, queries),
                   'zero-generator normalization')
            absorbing = (TRANSITIONS[:, 1] > TRANSITIONS[:, 0]).astype(float)
            p_absorbing = equilibrium_from_rates(absorbing)
            require(np.array_equal(p_absorbing, np.eye(1, 28, 27)[0]),
                    'Unique absorbing-state equilibrium was rejected or changed')
            disconnected = np.ones(len(TRANSITIONS))
            for i, (n, m) in enumerate(TRANSITIONS):
                if (n < 14) != (m < 14):
                    disconnected[i] = 0
            reject(lambda: equilibrium_from_rates(disconnected), 'nonunique stationary distribution')
    result = dict(passed=True, cases=summaries, positive_results_bitwise_unchanged=True,
                  invalid_and_missing_rejected=True, degenerate_generators_rejected=True,
                  unique_absorbing_state_accepted=True,
                  runtime_warnings_or_nan=False, existing_tables_regenerated=False)
    (out/'summary.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
