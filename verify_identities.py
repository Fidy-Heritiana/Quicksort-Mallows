#!/usr/bin/env python3
"""Exact checks of factorial, ordinary and centred moments, and their equations.

Python standard library only. Run: python3 verify_identities.py
All calculations use exact rational fractions, with no floating-point tolerance.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
MAX_N, MAX_ORDER = 10, 6
PARAMETERS = (Q(1, 3), Q(1, 2), Q(2, 3))


def falling_factorial(x, r):
    result = 1
    for i in range(r):
        result *= x-i
    return result


def quicksort_cost(word):
    if len(word) < 2:
        return 0
    pivot = word[0]
    return (len(word)-1
            + quicksort_cost([x for x in word[1:] if x < pivot])
            + quicksort_cost([x for x in word[1:] if x > pivot]))


def enumerate_permutations(max_n=7):
    """Integer histogram of (inversions, cost), independent of the distributional recurrence."""
    result, count = {0: {(0, 0): 1}}, 0
    for n in range(1, max_n+1):
        hist = defaultdict(int)
        for word in permutations(range(n)):
            inv = sum(word[i] > word[j] for i in range(n) for j in range(i+1, n))
            hist[inv, quicksort_cost(word)] += 1
            count += 1
        result[n] = hist
    return result, count


def cost_distributions(q):
    result = [{0: Q(1)}]
    for n in range(1, MAX_N+1):
        pmf = defaultdict(Q)
        for j in range(n):
            w = (1-q)*q**j/(1-q**n)
            for a, pa in result[j].items():
                for b, pb in result[n-1-j].items():
                    pmf[n-1+a+b] += w*pa*pb
        assert sum(pmf.values()) == 1
        result.append(dict(pmf))
    return result


def verify_parameter(q, enumeration):
    pmfs = cost_distributions(q)
    f = [[sum(w*falling_factorial(c, r) for c, w in law.items())
          for law in pmfs] for r in range(MAX_ORDER+1)]
    m = [[sum(w*c**r for c, w in law.items()) for law in pmfs]
         for r in range(MAX_ORDER+1)]
    mu = m[1]
    centred = [[sum(w*(c-mu[n])**r for c, w in law.items())
                for n, law in enumerate(pmfs)] for r in range(MAX_ORDER+1)]
    stirling = [[0]*(MAX_ORDER+1) for _ in range(MAX_ORDER+1)]
    stirling[0][0] = 1
    for r in range(1, MAX_ORDER+1):
        for k in range(1, r+1):
            stirling[r][k] = stirling[r-1][k-1]+k*stirling[r-1][k]

    checks = defaultdict(int)
    for n, histogram in enumeration.items():
        total = sum(count*q**inv for (inv, _), count in histogram.items())
        actual = defaultdict(Q)
        for (inv, c), count in histogram.items():
            actual[c] += count*q**inv/total
        assert dict(actual) == pmfs[n], ("enumeration", q, n)
        checks["distributions_compared_with_enumeration"] += 1

    for n in range(MAX_N+1):
        exact_mean = sum((n-k)*(1-q)*(1+q**k)/(1-q**(k+1))
                         for k in range(1, n))
        assert mu[n] == exact_mean
        assert centred[2][n] == f[2][n]+mu[n]-mu[n]**2
        checks["expectation_and_variance"] += 1
        for r in range(MAX_ORDER+1):
            assert m[r][n] == sum(stirling[r][k]*f[k][n] for k in range(r+1))
            assert centred[r][n] == sum(comb(r, k)*(-mu[n])**(r-k)*m[k][n]
                                       for k in range(r+1))
            checks["moment_conversions"] += 1

    # Coefficients of H_r: the binomial sum of z^i F_{r-i}^{(i)}.
    h = [[sum(comb(r, i)*falling_factorial(n, i)*f[r-i][n]
              for i in range(r+1)) for n in range(MAX_N+1)] for r in range(MAX_ORDER+1)]
    d = [[h[r][n]-f[r][n] for n in range(MAX_N+1)] for r in range(MAX_ORDER+1)]
    for r in range(MAX_ORDER+1):
        for n in range(MAX_N+1):
            direct = sum(w*falling_factorial(n+c, r)
                         for c, w in pmfs[n].items())
            assert h[r][n] == direct, ("chain rule", q, n, r)
            checks["chain_rule_derivatives"] += 1

    for r in range(1, MAX_ORDER+1):
        source = [Q(0)]*(MAX_N+1)
        for n in range(1, MAX_N+1):
            mixed = sum(comb(r, k)*h[k][j]*q**(n-1-j)*h[r-k][n-1-j]
                        for k in range(1, r) for j in range(n))
            source[n] = (1-q)*(sum(d[r][j]*q**(n-1-j)+q**j*d[r][j]
                                  for j in range(n))+mixed)
        for n in range(MAX_N+1):
            # Multiply the hierarchy by (1-z)^2 for an exact polynomial coefficient check.
            at = lambda values, k: values[k] if k >= 0 else Q(0)
            lhs = (1-q**n)*(at(f[r], n)-2*at(f[r], n-1)+at(f[r], n-2))
            rhs = at(source, n)-(1+q)*at(source, n-1)+q*at(source, n-2)
            assert lhs == rhs, ("factorial-moment hierarchy", q, n, r)
            checks["hierarchy_coefficients"] += 1

    # Moment recurrences: keep the factorial and centred expansions separate.
    for n in range(1, MAX_N+1):
        for r in range(2, MAX_ORDER+1):
            fact_value, centred_value = Q(0), Q(0)
            for j in range(n):
                k = n-1-j
                w = (1-q)*q**j/(1-q**n)
                toll = n-1+mu[j]+mu[k]-mu[n]
                for a in range(r+1):
                    for b in range(r-a+1):
                        c = r-a-b
                        coefficient = comb(r, a)*comb(r-a, b)
                        fact_value += w*coefficient*f[a][j]*f[b][k]*falling_factorial(n-1, c)
                        centred_value += w*coefficient*centred[a][j]*centred[b][k]*toll**c
            assert fact_value == f[r][n], ("factorial-moment recurrence", q, n, r)
            assert centred_value == centred[r][n], ("centred-moment recurrence", q, n, r)
            checks["recurrences_compared_with_distributions"] += 1

    # Check the summation identity on a separate convolution problem.
    b = [Q(0)]+[Q((-1)**n*(n*n+1), n+1) for n in range(1, MAX_N+1)]
    a = [Q(7, 3)]
    for n in range(1, MAX_N+1):
        a.append(sum((1-q)*q**j/(1-q**n)*a[n-1-j] for j in range(n))+b[n])
        assert a[n] == a[0]+b[n]+(1-q)*sum(b[k]/(1-q**(k+1)) for k in range(1, n))
        if n >= 2:
            assert a[n]-a[n-1] == b[n]-q*(1-q**(n-1))/(1-q**n)*b[n-1]
        checks["exact_summation_and_differences"] += 1
    return dict(checks)


def main():
    enumeration, count = enumerate_permutations()
    results = {str(q): verify_parameter(q, enumeration) for q in PARAMETERS}
    result = {"status": "passed", "arithmetic": "exact rational fractions",
              "max_recurrence_size": MAX_N, "max_moment_order": MAX_ORDER,
              "max_enumeration_size": 7, "permutations_enumerated": count,
              "parameters": results}
    dest = ROOT/"data"/"verification_identities.json"
    dest.parent.mkdir(exist_ok=True)
    dest.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
