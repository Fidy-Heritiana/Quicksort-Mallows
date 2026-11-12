#!/usr/bin/env python3
"""Targeted exact checks of the additional derivations in the detailed manuscript.

Run: python3 verify_details.py
These finite-size checks complement the proofs in the manuscript.
They are not a computer-assisted proof of the asymptotic results.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
import json

from verify_identities import PARAMETERS, MAX_N, quicksort_cost, cost_distributions

ROOT = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()


def verify_parameter(phi):
    pmfs = cost_distributions(phi)
    mu = [sum(w*c for c, w in law.items()) for law in pmfs]
    moments = [[sum(w*(c-mu[n])**r for c, w in law.items())
                for n, law in enumerate(pmfs)] for r in range(5)]
    counts = defaultdict(int)
    p = 1-phi
    for n in range(1, MAX_N+1):
        weights = [p*phi**j/(1-phi**n) for j in range(n)]
        ds = [n-1+mu[j]+mu[n-1-j]-mu[n] for j in range(n)]
        assert sum(w*d for w, d in zip(weights, ds)) == 0
        counts['exact_centring'] += 1
        # The seven terms explicitly displayed in Appendix C.
        source_four = Q(0)
        for j, w in enumerate(weights):
            k, d = n-1-j, ds[j]
            source_four += w*(
                moments[4][j] + 4*moments[3][j]*d
                + 6*moments[2][j]*moments[2][k]
                + 6*moments[2][j]*d*d + 4*moments[3][k]*d
                + 6*moments[2][k]*d*d + d**4)
        # Compare with a direct expectation under the full cost distribution.
        direct_source = moments[4][n]-sum(
            weights[j]*moments[4][n-1-j] for j in range(n))
        assert source_four == direct_source, ('seven terms of T4', phi, n)
        counts['fourth_moment_source'] += 1
        # Intervals for inverse-CDF sampling of the truncated geometric distribution.
        endpoints = [(1-phi**j)/(1-phi**n) for j in range(n+1)]
        assert endpoints[0] == 0 and endpoints[-1] == 1
        for j in range(n):
            assert endpoints[j+1]-endpoints[j] == weights[j]
            assert endpoints[j] < endpoints[j+1]
            counts['sampling_intervals'] += 1
        if n >= 2:
            sn = (1-phi**n)/p
            left = sn*mu[n]-(sn+phi**(n-1))*mu[n-1]
            right = n-2+sn+p*sum(phi**j*mu[j] for j in range(n-1))
            assert left == right, ('first elimination', phi, n)
            assert sn*(mu[n]-2*mu[n-1]+mu[n-2]) == 1+phi**(n-1)
            counts['two_eliminations_for_expectation'] += 1
    assert cost_distributions(1/phi) == pmfs, ('reciprocal-parameter symmetry', phi)
    counts['reciprocal_parameter_distributions'] = len(pmfs)
    if phi == Q(1, 2):
        assert mu[3] == Q(19, 7)
        assert moments[2][3] == Q(10, 49)
        assert moments[3][3] == Q(-30, 343)
        assert moments[4][3] == Q(190, 2401)
        counts['small_size_values'] = 4
    return dict(counts)


def main():
    counts = 0
    for n in range(1, 8):
        for word in permutations(range(n)):
            complement = tuple(n-1-x for x in word)
            assert quicksort_cost(word) == quicksort_cost(complement)
            counts += 1
    result = {
        'status': 'passed',
        'arithmetic': 'exact rational fractions',
        'max_size': MAX_N,
        'permutations_compared_with_their_complements': counts,
        'parameters': {str(phi): verify_parameter(phi) for phi in PARAMETERS},
    }
    (ROOT/'data').mkdir(exist_ok=True)
    (ROOT/'data'/'verification_details.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
