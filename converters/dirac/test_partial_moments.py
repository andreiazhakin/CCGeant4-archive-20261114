"""Independent quadrature check of affine-energy delta tail moments.

Runs deterministic power-law spectra by default. --pack additionally checks a
real converted delta file. Quadrature evaluates the mapped differential density
directly and does not use the analytic tail/moment code under test.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from cc_binary import read_pack
from cc_physics import delta_moments, delta_tail, integrate_spectrum

LOW_TRANSFER = 0.001


def independent_tmax(e):
    gamma = 1.0 + float(e)/931.494
    return 2.0*0.51099895*(gamma*gamma-1.0)


def synthetic_data():
    energy = np.array([30.0, 100.0, 2000.0])
    transfers = [np.geomspace(LOW_TRANSFER, independent_tmax(e), 73+i*4)
                 for i, e in enumerate(energy)]
    offsets = np.r_[0, np.cumsum([len(t) for t in transfers])]
    density = np.empty((11, int(offsets[-1])))
    for i, t in enumerate(transfers):
        for n in range(11):
            power = [-3., -2., -1., -.5, .25][(i+n) % 5]
            density[n, offsets[i]:offsets[i+1]] = (
                1.e-16*(i+1)*(1+n/20)*(t/LOW_TRANSFER)**power)
    tails = [np.empty_like(density) for _ in range(3)]
    for i, t in enumerate(transfers):
        sl = slice(offsets[i], offsets[i+1])
        for order, values in enumerate(integrate_spectrum(t, density[:, sl])):
            tails[order][:, sl] = values
    return dict(energy_MeV_u=energy, offsets=offsets,
                T_MeV=np.concatenate(transfers), D_cm2_MeV=density,
                tail0_cm2=tails[0], tail1_MeV_cm2=tails[1],
                tail2_MeV2_cm2=tails[2])


def independent_quadrature(data, e, n, cut, zp):
    en = data['energy_MeV_u']
    right = min(max(int(np.searchsorted(en, e, side='right')), 1), len(en)-1)
    left = right-1
    weight = np.log(e/en[left])/np.log(en[right]/en[left])
    n_base = min(n, data['D_cm2_MeV'].shape[0]-1)
    scale = ((zp-n)/(zp-n_base))**2
    result = np.zeros(3)
    quadrature_x, quadrature_w = np.polynomial.legendre.leggauss(8)
    for source, branch_weight in [(left, 1-weight), (right, weight)]:
        sl = slice(data['offsets'][source], data['offsets'][source+1])
        t = data['T_MeV'][sl]
        density = data['D_cm2_MeV'][n_base, sl]
        a = ((independent_tmax(e)-LOW_TRANSFER)
             /(independent_tmax(en[source])-LOW_TRANSFER))
        b = LOW_TRANSFER*(1-a)
        source_cut = (cut-b)/a
        lower = np.maximum(t[:-1], source_cut)
        upper = t[1:]
        keep = upper > lower
        lower, upper = lower[keep], upper[keep]
        if not len(lower):
            continue
        sample = ((upper+lower)[:, None]/2
                  +(upper-lower)[:, None]/2*quadrature_x)
        powers = np.log(density[1:]/density[:-1])/np.log(t[1:]/t[:-1])
        sampled_density = (density[:-1][keep, None]
                           *(sample/t[:-1][keep, None])**powers[keep, None])
        for order in range(3):
            result[order] += branch_weight*np.sum(
                (a*sample+b)**order * sampled_density * quadrature_w
                *(upper-lower)[:, None]/2)
    return scale*result


def expect_value_error(function, name):
    try:
        function()
    except ValueError:
        return
    raise AssertionError('Expected ValueError: '+name)


def check(data, zp, name):
    grid = data['energy_MeV_u']
    ncalc = data['D_cm2_MeV'].shape[0]
    energies = np.unique(np.r_[grid[0], grid[1], grid[len(grid)//2], grid[-1],
                              np.sqrt(grid[0]*grid[1]),
                              np.sqrt(grid[-2]*grid[-1]),
                              np.exp(.37*np.log(grid[-1])+.63*np.log(grid[0]))])
    charge_states = sorted({0, ncalc-1, min(ncalc, 27), 27})
    fractions = [0., .001, .1, .5, .99, 1-1.e-9]
    max_normalized_error = np.zeros(3)
    max_ordinary_relative_error = np.zeros(3)
    max_total_consistency_error = 0.
    cases = 0
    for energy in energies:
        for n in charge_states:
            total = delta_moments(data, energy, zp)[0, n]
            direct_total = independent_quadrature(data, energy, n, LOW_TRANSFER, zp)
            for fraction in fractions:
                cut = LOW_TRANSFER+fraction*(independent_tmax(energy)-LOW_TRANSFER)
                actual = delta_tail(data, energy, n, cut, zp)
                reference = independent_quadrature(data, energy, n, cut, zp)
                difference = np.abs(actual-reference)
                normalized = difference/direct_total
                max_normalized_error = np.maximum(max_normalized_error, normalized)
                if fraction < .999:
                    max_ordinary_relative_error = np.maximum(
                        max_ordinary_relative_error, difference/reference)
                if np.any(difference > 5.e-11*np.abs(reference)+2.e-13*direct_total):
                    raise AssertionError((name, energy, n, fraction, actual, reference))
                if np.any(actual < 0):
                    raise AssertionError('Negative tail moment')
                cases += 1
            tail_total = delta_tail(data, energy, n, LOW_TRANSFER, zp)
            consistency = float(np.max(np.abs(total/tail_total-1)))
            max_total_consistency_error = max(max_total_consistency_error, consistency)
            if consistency > 2.e-12:
                raise AssertionError('Total and cut-tail functions disagree')
            for cut in [independent_tmax(energy), 2*independent_tmax(energy)]:
                if np.any(delta_tail(data, energy, n, cut, zp) != 0):
                    raise AssertionError('Nonzero integral beyond Tmax')
    for energy in [float(grid[0])*.99, float(grid[-1])*1.01, float('nan')]:
        expect_value_error(lambda: delta_tail(data, energy, 0, LOW_TRANSFER, zp),
                           'energy coverage/finite')
    for n in [-1, 28]:
        expect_value_error(lambda: delta_tail(data, float(grid[0]), n, LOW_TRANSFER, zp),
                           'charge-state coverage')
    for cut in [-1., LOW_TRANSFER*.99]:
        expect_value_error(lambda: delta_tail(data, float(grid[0]), 0, cut, zp),
                           'cut below data coverage')
    return dict(dataset=name, Zp=zp, cut_cases=cases,
                max_error_normalized_to_full_moment=max_normalized_error.tolist(),
                max_ordinary_cut_relative_error=max_ordinary_relative_error.tolist(),
                max_total_consistency_relative_error=max_total_consistency_error,
                endpoint_and_domain_checks='passed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pack', type=Path, action='append', default=[])
    args = parser.parse_args()
    reports = [check(synthetic_data(), 29, 'constructed power laws')]
    for path in args.pack:
        header, data = read_pack(path)
        if header['kind'] != 1:
            raise ValueError('Expected a delta pack')
        reports.append(check(data, header['Zp'], str(path)))
    print(json.dumps({'passed': True, 'tests': reports}, indent=2, allow_nan=False))
