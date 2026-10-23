#!/usr/bin/env python3
"""Portable short-lemma checker for all-real-coefficient face exclusions.

Sparse rational polynomial identities + fixed-sign positive sums. No Z3/SymPy.
The real-root arguments for the small templates are in paper/paper-v1.tex, Section 5.
"""
from pathlib import Path
from collections import Counter
from itertools import permutations
import json
import sys

if not __debug__:
    raise RuntimeError("Run without -O or PYTHONOPTIMIZE; checks use assert")
import exact_check as q
NAMES=['s','t','c','q','r','d','z'];ZERO=(0,)*7


class P:
    def __init__(self,value=0):
        if isinstance(value,P):self.d=dict(value.d)
        elif isinstance(value,dict):self.d={k:q.F(v) for k,v in value.items() if v}
        else:self.d={ZERO:q.F(value)} if value else {}
    @staticmethod
    def variable(i):
        powers=list(ZERO);powers[i]=1;return P({tuple(powers):1})
    def __add__(self,other):
        other=P(other);out=dict(self.d)
        for k,v in other.d.items():out[k]=out.get(k,q.F(0))+v
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.d.items()})
    def __sub__(self,other):return self+-P(other)
    def __rsub__(self,other):return P(other)+-self
    def __mul__(self,other):
        other=P(other);out={}
        for a,x in self.d.items():
            for b,y in other.d.items():
                k=tuple(i+j for i,j in zip(a,b));out[k]=out.get(k,q.F(0))+x*y
        return P(out)
    __rmul__=__mul__
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        out=P(1)
        for _ in range(n):out=out*self
        return out
    def __eq__(self,other):return self.d==P(other).d
    def positive(self):
        # All active theta variables are strictly positive. z is unrestricted.
        return bool(self.d) and all(v>0 and k[-1]==0 for k,v in self.d.items())
    def record(self):return [{'powers':list(k),'coefficient':str(v)} for k,v in sorted(self.d.items())]


def determinant(matrix):
    total=P(0)
    for permutation in permutations(range(4)):
        term=P((-1)**sum(permutation[i]>permutation[j] for i in range(4) for j in range(i+1,4)))
        for i,j in enumerate(permutation):term=term*matrix[i][j]
        total=total+term
    return total


def data(case,swap=False):
    assert 1<=case['xmask']<=7 and 1<=case['ymask']<=7
    theta=[P.variable(j) if mask&(1<<i) else P(0) for mask,offset in [(case['xmask'],0),(case['ymask'],3)] for i,j in [(i,offset+i) for i in range(3)]]
    if swap:theta=theta[3:]+theta[:3]
    aa,bb=q.coefficients(theta)
    A,B,C,D,E,F=aa[3:]+bb[3:]
    assert (A+B+C)==theta[0] and (D+E+F)==theta[3]
    assert A.positive() and F.positive()
    return theta,(A,B,C,D,E,F)


def common_factor(case):
    theta,co=data(case);s,t,c,qq,r,d=theta;A,B,C,D,E,F=co
    if not(s==0 and qq==0):return None
    assert A+B+C==0 and D+E+F==0
    resultant=determinant([[A,B,C,P(0)],[P(0),A,B,C],[D,E,F,P(0)],[P(0),D,E,F]])
    assert resultant==0
    return {'only_type':None,'lemma':'COMMON_POINT_U_EQUALS_V','resultant_identity_zero':True}


def forced_axes(case):
    theta,co=data(case);s,t,c,qq,r,d=theta;A,B,C,D,E,F=co;required=[]
    if d==0:
        assert (A-E).positive()  # h'(0)=A+r>0, so root0 is simple
        assert E==-r
        required.append('P2H' if r==0 else 'H')
    if c==0:
        assert (F-B).positive()  # degree2, simple infinity
        assert B==-t
        required.append('P2H' if t==0 else 'H')
    return required


def small_face(case):
    common=common_factor(case)
    if common:return common
    x,y=case['xmask'],case['ymask']
    templates={(1,1):'AXIS_AXIS',(1,4):'AXIS_DIFFERENCE_SQUARE',(1,2):'AXIS_CROSS'}
    swap=False
    if (x,y) not in templates:
        assert (y,x) in templates,('not a small-face template',x,y)
        x,y=y,x;swap=True
    theta,co=data(case,swap);A,B,C,D,E,F=co;z=P.variable(6)
    pp=A+B*z+C*z*z;qp=D+E*z+F*z*z
    hh=z*pp-qp;nn=E+(2*F-B)*z-2*C*z*z
    assert C==0 and A.positive()
    template=templates[x,y]
    if template=='AXIS_AXIS':
        a,b=A,F
        assert a.positive() and b.positive()
        assert pp==a and qp==b*z*z and hh==z*(a-b*z) and nn==2*b*z
        code='P1,P2H,P2H'
    elif template=='AXIS_DIFFERENCE_SQUARE':
        a,b=A,D
        assert a.positive() and b.positive()
        assert pp==a and qp==b*(1-z)**2
        assert hh==-b*z*z+(a+2*b)*z-b
        delta=(a+2*b)**2-4*b*b
        assert delta==a*a+4*a*b and (delta-a*a).positive()
        assert nn==-2*b+2*b*z
        # Δ>a²>0; finite rho=1±sqrt(Δ)/a are H/P1, infinity=P2H.
        code='P1,H,P2H'
    else:
        a,b=A,F
        assert a.positive() and b.positive()
        assert pp==a and qp==b*z*(z-1) and hh==z*(a+b-b*z)
        assert nn==-b+2*b*z
        assert (a+b).positive() and (2*a+b).positive()
        code='P1,H,P2H'
    resultant=determinant([[A,B,C,P(0)],[P(0),A,B,C],[D,E,F,P(0)],[P(0),D,E,F]])
    assert resultant==a*a*b*b and resultant.positive()
    return {'only_type':code,'lemma':template,'species_swap':swap,'resultant_positive':True,
            'scope':'all strictly positive active real coefficients on this face; no normalization needed'}


def check_unsat(cert):
    assert cert['schema']=='R022.face-unsat.v1'
    case=cert['case'];assert case['type'] in q.NINE
    proof=common_factor(case)
    if proof:return {'status':'PROVED_UNSAT',**proof}
    required=forced_axes(case)
    target=Counter(case['type'].split(','));needed=Counter(required)
    if any(target[k]<n for k,n in needed.items()):
        return {'status':'PROVED_UNSAT','lemma':'FORCED_SIMPLE_AXIS_KINDS','required_simple_kinds':required}
    try:proof=small_face(case)
    except AssertionError:raise NotImplementedError('no portable proof template for this face/type')
    assert proof['only_type']!=case['type'],'not an impossible type on this face'
    return {'status':'PROVED_UNSAT',**proof}


def lower_seven(cert):
    assert cert['schema']=='R022.lower-seven-eight.v1'
    base=cert['base_lower_six'];q.lower_bound(base)
    low=[r for r in base['faces'] if r['terms']<8]
    assert len(low)==8
    assert {(r['xmask'],r['ymask']) for r in cert['low_faces']}=={(r['xmask'],r['ymask']) for r in low}
    proofs=[]
    for row in low:
        proof=small_face(row)
        assert proof['only_type'] in [None,'P1,H,P2H','P1,P2H,P2H']
        proofs.append({'xmask':row['xmask'],'ymask':row['ymask'],**proof})
    types=sorted(q.NINE-{'P1,H,P2H','P1,P2H,P2H'})
    assert sorted(cert['types'])==types and cert['bound']==8
    return {'status':'PASS','bound':8,'types':types,'all_sub8_faces':proofs,
            'scope':'all real coefficients satisfying setup; proof is not a rational sampling argument'}


if __name__=='__main__':
    cert=json.loads(Path(sys.argv[1]).read_text())
    print(json.dumps(lower_seven(cert) if cert['schema']=='R022.lower-seven-eight.v1' else check_unsat(cert),ensure_ascii=False))
