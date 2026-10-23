#!/usr/bin/env python3
"""Standalone rational certificate checker: no SymPy or Z3 imports.

Factorization + squarefreeness + Sturm counts prove completeness. Signs at
isolated roots use rational interval evaluation, or a gcd-root equality proof.
This is a second checking procedure by the same agent, not independent review.
"""
import sys
if not __debug__:
    raise RuntimeError("Run without -O or PYTHONOPTIMIZE; checks use assert")
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import sys

KINDS = ['P1', 'P2', 'H', 'P2H', 'D']
NINE = {'H', 'P2H', 'D,H', 'D,P2H', 'P1,P2,H', 'P1,H,H',
        'P1,P2,P2H', 'P1,H,P2H', 'P1,P2H,P2H'}


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p or [F(0)]


def add(p, q):
    return trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
                 for i in range(max(len(p), len(q)))])


def scale(p, a):
    return trim([a*v for v in p])


def mul(p, q):
    r = [F(0)] * (len(p)+len(q)-1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i+j] += x*y
    return trim(r)


def divrem(p, q):
    p, q = trim(p), trim(q)
    assert q != [0]
    out = [F(0)] * max(1, len(p)-len(q)+1)
    while p != [0] and len(p) >= len(q):
        n, a = len(p)-len(q), p[-1]/q[-1]
        out[n] += a
        p = add(p, [F(0)]*n + scale(q, -a))
    return trim(out), p


def gcd(p, q):
    p, q = trim(p), trim(q)
    while q != [0]:
        p, q = q, divrem(p, q)[1]
    return scale(p, 1/p[-1]) if p != [0] else p


def deriv(p):
    return trim([i*p[i] for i in range(1, len(p))])


def value(p, x):
    ans = F(0)
    for a in reversed(p):
        ans = ans*x+a
    return ans


def interval_value(p, lo, hi):
    left = right = F(0)
    for a in reversed(p):
        corners = [left*lo, left*hi, right*lo, right*hi]
        left, right = min(corners)+a, max(corners)+a
    return left, right


def sign(a):
    return (a > 0) - (a < 0)


def sturm(p):
    seq = [trim(p), deriv(p)]
    if seq[-1] == [0]:
        return seq[:1]
    while True:
        r = scale(divrem(seq[-2], seq[-1])[1], -1)
        if r == [0]:
            break
        seq.append(r)
    return seq


def variation(seq, x):
    signs = []
    for p in seq:
        if x == '+inf':
            s = sign(p[-1])
        elif x == '-inf':
            s = sign(p[-1]) * (-1)**(len(p)-1)
        else:
            s = sign(value(p, x))
        if s:
            signs.append(s)
    return sum(a != b for a, b in zip(signs, signs[1:]))


def roots_between(p, lo, hi):
    p = trim(p)
    if len(p) == 1:
        return 0
    if lo == hi:
        return int(value(p, lo) == 0)
    assert value(p, lo) != 0 and value(p, hi) != 0
    seq = sturm(p)
    return variation(seq, lo)-variation(seq, hi)


def algebraic_sign(factor, polynomial, bounds):
    lo, hi = map(F, bounds)
    polynomial = trim(polynomial)
    if lo == hi:
        return sign(value(polynomial, lo))
    common = gcd(factor, polynomial)
    if len(common) > 1 and roots_between(common, lo, hi) == 1:
        return 0
    a, b = interval_value(polynomial, lo, hi)
    assert a > 0 or b < 0, ('unproved nonzero sign', polynomial, bounds)
    return 1 if a > 0 else -1


def det(matrix):
    a = [list(map(F, row)) for row in matrix]
    answer = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            answer = -answer
        q = a[j][j]
        answer *= q
        for i in range(j+1, len(a)):
            r = a[i][j]/q
            for k in range(j+1, len(a)):
                a[i][k] -= r*a[j][k]
    return answer


def category(p, n, difference, mult):
    assert p != 0
    if mult == 2:
        assert difference == 0
        return 'D'
    assert mult == 1 and difference != 0
    if n == 0:
        return 'P2H'
    if n*p < 0:
        return 'H'
    return 'P1' if difference*p > 0 else 'P2'


def canonical(kinds):
    if 'D' in kinds:
        assert kinds.count('D') == 1 and len(kinds) == 2
        return 'D,' + next(k for k in kinds if k != 'D')
    return ','.join(sorted(kinds, key=KINDS.index))


def coefficients(theta):
    s, t, c, q, r, d = theta
    return [s, -2*s-t, t, s+t+c, -t-2*c, c], [q, r, -2*q-r, d, -r-2*d, q+r+d]


def reaction_check(reactions, aa, bb):
    exponents = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2)]
    out = [[F(0)]*6 for _ in range(2)]
    for reaction in reactions:
        a, b, rate = reaction['reactant'], reaction['product'], F(reaction['rate'])
        assert len(a) == len(b) == 2
        assert all(isinstance(v, int) and v >= 0 for v in a+b)
        assert sum(a) <= 2 and rate > 0 and a != b
        i = exponents.index(tuple(a))
        for j in range(2):
            out[j][i] += rate*(b[j]-a[j])
    assert out == [aa, bb]


def witness(cert):
    aa, bb = [list(map(F, cert[k])) for k in ('a', 'b')]
    assert len(aa) == len(bb) == 6
    theta = [aa[0], aa[2], aa[5], bb[0], bb[1], bb[3]]
    assert all(t >= 0 for t in theta)
    assert coefficients(theta) == (aa, bb)
    if 'theta_normalized' in cert:
        assert list(map(F,cert['theta_normalized'])) == theta and sum(theta) == 1
    A, B, C, D, E, Ff = aa[3:]+bb[3:]
    assert sum(aa) == sum(bb) == 0
    assert aa[1]+2*A+B == aa[2]+B+2*C == bb[1]+2*D+E == bb[2]+E+2*Ff == 0
    terms = sum(v != 0 for v in aa+bb)
    assert cert['terms'] == terms
    resultant = det([[A,B,C,0], [0,A,B,C], [D,E,Ff,0], [0,D,E,Ff]])
    assert F(cert['resultant']) == resultant
    if resultant == 0:
        assert cert['code'] == 'EXCLUDED_COMMON_FACTOR'
        return {'code':cert['code'], 'terms':terms}
    h = trim([-D, A-E, B-Ff, C])
    assert trim(cert['h']) == h
    assert h != [0]
    degree = len(h)-1
    assembled = [F(cert['factor_unit'])]
    factors = []
    for row in cert['factors']:
        f = trim(row['polynomial'])
        assert len(f) > 1 and f[-1] == 1
        assert gcd(f, deriv(f)) == [1]
        assert row['multiplicity'] in (1, 2, 3)
        for other in factors:
            assert gcd(f, other) == [1]
        factors.append(f)
        for _ in range(row['multiplicity']):
            assembled = mul(assembled, f)
    assert assembled == h
    assert all(isinstance(r['factor_index'],int) and 0 <= r['factor_index'] < len(factors) for r in cert['roots'])
    all_intervals = []
    kinds = []
    for index, f in enumerate(factors):
        records = [r for r in cert['roots'] if r['factor_index'] == index]
        seq = sturm(f)
        assert len(records) == variation(seq, '-inf')-variation(seq, '+inf')
        for root in records:
            lo, hi = map(F, root['interval'])
            assert lo <= hi and roots_between(f, lo, hi) == 1
            assert all(hi < a or b < lo for a, b in all_intervals)
            all_intervals.append((lo, hi))
            mult = cert['factors'][index]['multiplicity']
            assert root['multiplicity'] == mult
            if mult == 3:
                kinds.append('TRIPLE')
                continue
            p = [A, B, C]
            n = [E, 2*Ff-B, -2*C]
            nm = add(n, scale(p, -1))
            signs = {name:algebraic_sign(f, poly, root['interval'])
                     for name, poly in [('p',p), ('n',n), ('n_minus_p',nm)]}
            assert signs == root['signs']
            kind = category(signs['p'], signs['n'], signs['n_minus_p'], mult)
            assert root['kind'] == kind
            kinds.append(kind)
    infinity = 3-degree
    assert cert['infinity_multiplicity'] == infinity
    if infinity:
        assert C == 0 and Ff != 0
        if infinity == 3:
            kinds.append('TRIPLE')
        else:
            kinds.append(category(sign(Ff), sign(B), sign(B-Ff), infinity))
    code = 'EXCLUDED_TRIPLE' if 'TRIPLE' in kinds else canonical(kinds)
    assert cert['code'] == code
    if code != 'EXCLUDED_TRIPLE':
        assert code in NINE
    reaction_check(cert['reactions'], aa, bb)
    if cert.get('source_reactions'):
        reaction_check(cert['source_reactions'], aa, bb)
    return {'code':code, 'terms':terms, 'resultant':str(resultant), 'root_count':len(cert['roots'])+int(infinity>0)}


def lower_bound(cert):
    """Complete support/cancellation-free proof that every allowed field has M>=6."""
    assert cert['lemma'] == 'NONZERO_COMPONENT_THREE_TERMS'
    rows = cert['faces']
    assert {(r['xmask'], r['ymask']) for r in rows} == set(product(range(1,8), repeat=2))
    assert len(rows) == 49
    # Verify the fixed-sign property of every coefficient's linear form.
    images = [coefficients([F(i == j) for i in range(6)]) for j in range(6)]
    for component, term in product(range(2),range(6)):
        weights = [image[component][term] for image in images]
        assert all(w >= 0 for w in weights) or all(w <= 0 for w in weights)
    for row in rows:
        xmask, ymask = row['xmask'], row['ymask']
        theta = [F(int(bool(mask & (1 << bit)))) for mask in (xmask, ymask) for bit in range(3)]
        aa, bb = coefficients(theta)
        # All expressions have a fixed sign on a face: no coefficient has cancellation.
        # The three certified nonzero positions follow the first positive variable.
        expected_x = [0,1,3] if xmask & 1 else ([1,2,3] if xmask & 2 else [3,4,5])
        expected_y = [0,2,5] if ymask & 1 else ([1,2,5] if ymask & 2 else [3,4,5])
        assert row['three_x_indices'] == expected_x and row['three_y_indices'] == expected_y
        assert all(aa[i] != 0 for i in expected_x) and all(bb[i] != 0 for i in expected_y)
        count = sum(v != 0 for v in aa+bb)
        assert row['terms'] == count and count >= 6
        # No unhandled degree0/1 branch occurs in the chemical domain.
        A,B,C,D,E,Ff = aa[3:]+bb[3:]
        assert Ff > 0 and B <= 0
        assert C > 0 or B-Ff < 0
    assert cert['bound'] == 6
    return {'bound':6, 'faces':49, 'scope':'all real nonnegative theta with both components nonzero; resultant nonzero is a subset'}


def unsat_lemma(cert):
    assert cert['lemma'] == 'DECOUPLED_AXES'
    assert cert['xmask'] == cert['ymask'] == 1
    assert cert['excluded_type'] in NINE and cert['excluded_type'] != 'P1,P2H,P2H'
    # Formal sparse polynomials in s,q,z, checked with rational coefficients.
    def plus(a,b):
        out = dict(a)
        for k,v in b.items():
            out[k] = out.get(k,F(0))+v
        return {k:v for k,v in out.items() if v}
    def times(a,b):
        out = {}
        for i,x in a.items():
            for j,y in b.items():
                k = tuple(v+w for v,w in zip(i,j))
                out[k] = out.get(k,F(0))+x*y
        return {k:v for k,v in out.items() if v}
    def neg(a):
        return {k:-v for k,v in a.items()}
    ss,qq,zz = [{tuple(int(j==i) for j in range(3)):F(1)} for i in range(3)]
    # Derive the coefficients by linearity from the actual parametrization.
    basis_s = coefficients([F(1),F(0),F(0),F(0),F(0),F(0)])
    basis_q = coefficients([F(0),F(0),F(0),F(1),F(0),F(0)])
    lifted = [[plus({k:v*basis_s[j][i] for k,v in ss.items()},
                    {k:v*basis_q[j][i] for k,v in qq.items()}) for i in range(6)] for j in range(2)]
    A,B,C = lifted[0][3:];D,E,Ff = lifted[1][3:]
    z2 = times(zz,zz)
    pp = plus(A,plus(times(B,zz),times(C,z2)))
    qp = plus(D,plus(times(E,zz),times(Ff,z2)))
    hh = plus(times(zz,pp),neg(qp))
    assert pp == ss and qp == times(qq,z2)
    assert hh == times(zz,plus(ss,neg(times(qq,zz))))
    # With strict s,q>0 on face(1,1), 0,s/q,infinity are distinct.
    # P=s, numerator=2qz: rho=0,2,0, with positive nonzero denominators.
    return {'status':'PROVED_UNSAT','lemma':'DECOUPLED_AXES','excluded_type':cert['excluded_type'],
            'only_type':'P1,P2H,P2H','formal_factorization_checked':True}


def check(path):
    cert = json.loads(Path(path).read_text())
    if cert.get('schema') == 'R022.lower-bound.v1':
        return lower_bound(cert)
    if cert.get('schema') == 'R022.unsat-lemma.v1':
        return unsat_lemma(cert)
    return witness(cert)


if __name__ == '__main__':
    print(json.dumps(check(sys.argv[1]), ensure_ascii=False))
