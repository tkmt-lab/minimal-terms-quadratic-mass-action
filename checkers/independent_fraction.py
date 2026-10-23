#!/usr/bin/env python3
"""Independent verifier: standard library, rational elimination and quadratic fields.

Does NOT import sympy_certificate or read its output. Root enumeration uses the
rational root theorem + the quadratic formula, not numerical approximations.
Generic lower-bound identities use Laurent polynomial arithmetic with sqrt(D)^2
= a^2+4ab; positivity is justified by the accompanying written proof.
"""
import sys
if not __debug__:
    raise RuntimeError("Run without -O or PYTHONOPTIMIZE; checks use assert")
from fractions import Fraction as F
from itertools import product, permutations
from math import comb, gcd, isqrt, lcm
from pathlib import Path
import copy
import json
import resource
import time

ROOT=Path(__file__).resolve().parents[1]
MON=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
FREE=(0,2,5,6,7,9)
NAMES=('s','t','c','q','r','d')
EXPECTED={1:['H'],4:['P2H'],5:['D','H'],8:['D','P2H'],
          11:['P1','P2','H'],12:['P1','H','H'],14:['P1','P2','P2H'],
          15:['P1','H','P2H'],16:['P1','P2H','P2H']}

class Q:
    """a+b sqrt(d), d>=0; every comparison uses rational squares."""
    def __init__(self,a=0,b=0,d=0):
        self.a,self.b,self.d=F(a),F(b),F(d)
        assert self.d>=0
        n,m=isqrt(self.d.numerator),isqrt(self.d.denominator)
        if n*n==self.d.numerator and m*m==self.d.denominator:
            self.a+=self.b*F(n,m);self.b=F(0)
        if not self.b:self.d=F(0)
    def __add__(self,o):
        o=o if isinstance(o,Q) else Q(o)
        assert not self.d or not o.d or self.d==o.d
        return Q(self.a+o.a,self.b+o.b,self.d or o.d)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b,self.d)
    def __sub__(self,o):return self+-o if isinstance(o,Q) else self+Q(-F(o))
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        o=o if isinstance(o,Q) else Q(o)
        assert not self.d or not o.d or self.d==o.d
        d=self.d or o.d
        return Q(self.a*o.a+self.b*o.b*d,self.a*o.b+self.b*o.a,d)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=o if isinstance(o,Q) else Q(o)
        norm=o.a*o.a-o.b*o.b*o.d
        assert norm
        return self*Q(o.a/norm,-o.b/norm,o.d)
    def __rtruediv__(self,o):return Q(o)/self
    def sign(self):
        a,b=self.a,self.b
        if not b:return (a>0)-(a<0)
        if not a:return (b>0)-(b<0)
        if (a>0)==(b>0):return 1 if a>0 else -1
        delta=a*a-b*b*self.d
        return ((delta>0)-(delta<0))*(1 if a>0 else -1)
    def __eq__(self,o):
        if not isinstance(o,(Q,F,int)):return False
        return (self-o).sign()==0
    def __str__(self):
        return str(self.a) if not self.b else f'({self.a})+({self.b})*sqrt({self.d})'

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p
def evaluate(p,z):
    out=0
    for a in reversed(p):out=out*z+a
    return out
def deriv(p):return [i*p[i] for i in range(1,len(p))] or [F(0)]
def divide_linear(p,r):
    q=[F(0)]*(len(p)-1);q[-1]=p[-1]
    for i in range(len(q)-2,-1,-1):q[i]=p[i+1]+r*q[i+1]
    assert p[0]+r*q[0]==0
    return trim(q)
def divisors(n):
    n=abs(n)
    return {d for i in range(1,isqrt(n)+1) if n%i==0 for d in (i,n//i)}
def real_roots(p):
    """Returns distinct roots with multiplicities, and nonreal degree."""
    p=trim(list(map(F,p))); roots=[];nonreal=0
    if len(p)==4:
        scale=lcm(*(a.denominator for a in p));ip=[int(a*scale) for a in p]
        cand=[F(0)] if not ip[0] else sorted({F(k*n,d) for n in divisors(ip[0])
                                                           for d in divisors(ip[-1]) for k in (-1,1)})
        rr=next((r for r in cand if evaluate(p,r)==0),None)
        assert rr is not None,'This small certificate checker requires a rational cubic root'
        roots.append(Q(rr));p=divide_linear(p,rr)
    if len(p)==3:
        c,b,a=p;disc=b*b-4*a*c
        if disc<0:nonreal=2
        else:roots += [Q(-b/(2*a),sgn/(2*a),disc) for sgn in (-1,1)]
    elif len(p)==2:roots.append(Q(-p[0]/p[1]))
    else:assert len(p)==1
    grouped=[]
    for rr in roots:
        old=next((v for v in grouped if v[0]==rr),None)
        if old is None:grouped.append([rr,1])
        else:old[1]+=1
    return grouped,nonreal

def solve(matrix,rhs):
    a=[list(map(F,row))+[F(v)] for row,v in zip(matrix,rhs)]
    n=len(rhs)
    for j in range(n):
        k=next(k for k in range(j,n) if a[k][j])
        a[j],a[k]=a[k],a[j]
        pivot=a[j][j];a[j]=[z/pivot for z in a[j]]
        for k in range(n):
            if k!=j:
                f=a[k][j];a[k]=[z-f*w for z,w in zip(a[k],a[j])]
    return [row[-1] for row in a]

def derive_basis():
    # Six equilibrium/Jacobian equations from monomial exponents at (1,1).
    rows=[]
    for component in (0,1):
        for op in (None,0,1):
            rows.append([F(1 if op is None else MON[j%6][op]) if j//6==component else F(0)
                         for j in range(12)])
    dependent=[j for j in range(12) if j not in FREE]
    basis=[]
    for selected in FREE:
        values=solve([[row[j] for j in dependent] for row in rows],[-row[selected] for row in rows])
        vector=[F(0)]*12;vector[selected]=F(1)
        for j,k in zip(dependent,values):vector[j]=k
        assert all(sum(a*b for a,b in zip(row,vector))==0 for row in rows)
        basis.append(vector)
    return basis

def from_parameters(values,basis):
    assert len(values)==6
    return [sum(F(v)*base[j] for v,base in zip(values,basis)) for j in range(12)]

def from_homogeneous(values):
    # Direct binomial translation P(x-1,y-1), Q(x-1,y-1).
    out=[]
    for offset in (0,3):
        c={ij:F(0) for ij in MON}
        for value,(i,j) in zip(map(F,values[offset:offset+3]),MON[3:]):
            for ii in range(i+1):
                for jj in range(j+1):
                    c[ii,jj]+=value*comb(i,ii)*comb(j,jj)*(-1)**(i+j-ii-jj)
        out.extend(c[ij] for ij in MON)
    return out

def determinant(m):
    result=F(0)
    for p in permutations(range(len(m))):
        sign=(-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
        value=F(sign)
        for i,j in enumerate(p):value*=m[i][j]
        result+=value
    return result

def label(rho,mult):
    if mult==2:
        assert rho==1;return 'D'
    assert mult==1 and rho!=1
    if rho.sign()==0:return 'P2H'
    if rho.sign()<0:return 'H'
    return 'P1' if (rho-1).sign()>0 else 'P2'

def verify_coefficients(coeff,expected_type=None,expected_M=None):
    assert len(coeff)==12
    assert all(coeff[j]>=0 for j in FREE),'kinetic sign'
    for offset in (0,6):
        assert sum(coeff[offset:offset+6])==0,'equilibrium'
        for direction in (0,1):
            assert sum(coeff[offset+j]*MON[j][direction] for j in range(6))==0,'Jacobian'
    A,B,C=coeff[3:6];D,E,F_=coeff[9:12]
    R=determinant([[A,B,C,0],[0,A,B,C],[D,E,F_,0],[0,D,E,F_]])
    assert R!=0,'homogeneous common factor'
    h=trim([-D,A-E,B-F_,C]);degree=len(h)-1
    assert degree in (2,3)
    roots,nonreal=real_roots(h)
    assert sum(mult for _,mult in roots)+nonreal==degree,'exhaustiveness'
    lines=[]
    for rr,mult in roots:
        assert evaluate(h,rr)==0
        dh=h
        for _ in range(mult):
            assert evaluate(dh,rr)==0;dh=deriv(dh)
        assert evaluate(dh,rr)!=0,'multiplicity'
        p=evaluate([A,B,C],rr);qq=evaluate([D,E,F_],rr)
        assert p.sign()>0 and qq==rr*p
        # Quotient derivative, not the producer's h'/p as primary calculation.
        rho=(evaluate([E,2*F_],rr)*p-qq*evaluate([B,2*C],rr))/(p*p)
        assert rho==1-evaluate(deriv(h),rr)/p
        lines.append({'root':str(rr),'multiplicity':mult,'rho':str(rho),
                      'label':label(rho,mult),'radial':str(p)})
    if degree==2:
        assert C==0 and F_>0
        rho=Q(B/F_)
        lines.append({'root':'infinity','multiplicity':1,'rho':str(rho),
                      'label':label(rho,1),'radial':str(F_)})
    found=[k for k,ls in EXPECTED.items() if sorted(ls)==sorted(a['label'] for a in lines)]
    assert len(found)==1,'unsupported type'
    if expected_type is not None:assert found[0]==expected_type,'type'
    M=sum(bool(a) for a in coeff)
    if expected_M is not None:assert M==expected_M,'term count'
    return {'type':found[0],'M':M,'coefficients':list(map(str,coeff)),
            'R':str(R),'h_ascending':list(map(str,h)),'lines':lines,
            'nonreal_degree':nonreal}

class L:
    """Small Laurent ring Q[a,+/-1,b,+/-1,delta]/(delta^2-a^2-4ab)."""
    def __init__(self,v=0,terms=None):
        self.t={k:F(x) for k,x in (terms if terms is not None else {(0,0,0):v}).items() if x}
    def __add__(self,o):
        o=o if isinstance(o,L) else L(o);t=self.t.copy()
        for k,v in o.t.items():t[k]=t.get(k,0)+v
        return L(terms=t)
    __radd__=__add__
    def __neg__(self):return L(terms={k:-v for k,v in self.t.items()})
    def __sub__(self,o):return self+(-o if isinstance(o,L) else -F(o))
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        o=o if isinstance(o,L) else L(o);t={}
        for (i,j,k),v in self.t.items():
            for (ii,jj,kk),w in o.t.items():
                exponent=(i+ii,j+jj,k+kk)
                if exponent[2]==2:
                    for a,b,co in [(2,0,1),(1,1,4)]:
                        e=(exponent[0]+a,exponent[1]+b,0);t[e]=t.get(e,0)+co*v*w
                else:t[exponent]=t.get(exponent,0)+v*w
        return L(terms=t)
    __rmul__=__mul__
    def __eq__(self,o):return not (self-o).t

def check_generic_families():
    a=L(terms={(1,0,0):1});b=L(terms={(0,1,0):1})
    ia=L(terms={(-1,0,0):1});ib=L(terms={(0,-1,0):1})
    delta=L(terms={(0,0,1):1})
    assert a*ia==1 and b*ib==1 and delta*delta-a*a==4*a*b
    families=[('I',[0,0,b],[L(0),a*ib],[L(0),L(2)]),
              ('II',[b,-2*b,b],[(a+2*b-delta)*ib*F(1,2),(a+2*b+delta)*ib*F(1,2)],
               [1-delta*ia,1+delta*ia]),
              ('III',[0,-b,b],[L(0),1+a*ib],[-b*ia,2+b*ia])]
    records=[]
    for name,qq,roots,rhos in families:
        for z,rho in zip(roots,rhos):
            assert a*z-evaluate(qq,z)==0
            assert evaluate(deriv(qq),z)*ia==rho
        # Factoring h exhausts its two simple finite roots, for every a,b>0.
        assert -b*(roots[0]+roots[1]) == -(a-qq[1])
        assert -b*roots[0]*roots[1] == -qq[0]
        # Homogeneous resultant for P=a*u^2 is a^2*b^2>0.
        generic_matrix=[[a,0,0,0],[0,a,0,0],[qq[0],qq[1],b,0],[0,qq[0],qq[1],b]]
        det=L(0)
        for p in permutations(range(4)):
            term=L((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4)))
            for i,j in enumerate(p):term=term*generic_matrix[i][j]
            det=det+term
        assert det==a*a*b*b
        records.append({'family':name,'identities':'PASS','type':16 if name=='I' else 15})
    return records

def lower_bound(basis):
    for j in range(12):
        signs={1 if b[j]>0 else -1 for b in basis if b[j]}
        assert len(signs)==1,'linear coefficient cancellation would invalidate support enumeration'
    # Canonical sparse homogeneous forms from the independently solved basis.
    assert basis[0][3:6]==[1,0,0] and basis[1][3:6]==[1,-1,0]
    assert basis[2][3:6]==[1,-2,1]
    assert basis[3][9:12]==[0,0,1] and basis[4][9:12]==[0,-1,1]
    assert basis[5][9:12]==[1,-2,1]
    cells=[];small=[]
    for bits in product((0,1),repeat=6):
        # This evaluation counts an entire sign cell because cancellation was excluded above.
        coeff=from_parameters(bits,basis)
        mx=sum(bool(a) for a in coeff[:6]);my=sum(bool(a) for a in coeff[6:])
        active=[NAMES[i] for i,b in enumerate(bits) if b]
        row={'active':active,'Mx':mx,'My':my,'M':mx+my}
        if not mx or not my:row['case']='zero-component'
        elif mx+my>=8:row['case']='M-at-least-8'
        else:
            key=tuple(active)
            cases={('s','q'):'I',('s','d'):'II',('c','q'):'II-swapped',
                   ('s','r'):'III',('t','q'):'III-swapped'}
            if key in cases:row['case']=cases[key]
            else:
                assert key in [('c','d'),('c','r'),('t','d')]
                # Every active generator vanishes on u=v, identically in its positive weight.
                assert all(sum(basis[i][3:6])==0 and sum(basis[i][9:12])==0
                           for i,b in enumerate(bits) if b)
                row['case']='common-factor'
            small.append(row)
        cells.append(row)
    assert len(cells)==64 and len(small)==8
    assert sorted(a['M'] for a in small)==[6,6,6,6,7,7,7,7]
    assert sum(a['case']=='common-factor' for a in small)==3
    assert all(a['Mx']>=3 and a['My']>=3 for a in cells if a['Mx'] and a['My'])
    return {'derived_basis':[[str(x) for x in row] for row in basis],
            'cells':cells,'small_cells':small,'generic_families':check_generic_families(),
            'minima':{str(k):(6 if k in (15,16) else 8) for k in EXPECTED}}

def mutation_tests(basis,witnesses):
    base=from_parameters(witnesses[-1]['parameters'],basis)
    tests=[]
    for name,coeff,ty,m in [
        ('negative-kinetic-sign',from_parameters([-1,0,0,1,0,0],basis),16,6),
        ('common-factor',from_parameters([0,0,1,0,0,1],basis),16,6),
        ('wrong-type',base,15,6),('wrong-count',base,16,5),
        ('broken-equilibrium',[base[0]+1]+base[1:],16,6)]:
        try:verify_coefficients(coeff,ty,m)
        except AssertionError:tests.append(name)
        else:raise AssertionError(f'mutation survived: {name}')
    # Primitive checks protect both sign boundaries and exact multiplicity.
    assert Q(1,-1,5).sign()<0 and Q(1,1,5).sign()>0
    assert Q(1,-1,1).sign()==0
    rr,n=real_roots([1,-2,1]);assert n==0 and len(rr)==1 and rr[0][1]==2
    assert Q(1,1,5)*Q(1,-1,5)==-4
    return {'rejected':tests,'arithmetic_boundary_checks':'PASS'}

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(120,125))
    start=time.monotonic();data=json.loads((ROOT/'certificates/independent-witnesses.json').read_text());basis=derive_basis()
    out={'method':'stdlib rational elimination, rational-root theorem, quadratic fields, Laurent identities'}
    for key in ('witnesses','paper_controls'):
        out[key]=[]
        for e in data[key]:
            coeff=(from_parameters(e['parameters'],basis) if 'parameters' in e else from_homogeneous(e['ABCDEF']))
            assert coeff==from_parameters([coeff[j] for j in FREE],basis)
            out[key].append({'id':e['id'],**verify_coefficients(coeff,e['type'],e['M'])})
    out['excluded_controls']=[]
    for e in data['excluded_controls']:
        coeff=from_homogeneous(e['ABCDEF'])
        try:verify_coefficients(coeff)
        except AssertionError as err:
            assert str(err)=='homogeneous common factor'
            out['excluded_controls'].append({'id':e['id'],'excluded':True})
        else:raise AssertionError('excluded control accepted')
    out['lower_bound']=lower_bound(basis)
    out['mutation_tests']=mutation_tests(basis,data['witnesses'])
    out['elapsed_seconds']=time.monotonic()-start
    out['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    out['status']='PASS'
    (ROOT/'work').mkdir(exist_ok=True)
    (ROOT/'work/independent-fraction-check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','elapsed_seconds','peak_rss_bytes')}))

if __name__=='__main__':main()
