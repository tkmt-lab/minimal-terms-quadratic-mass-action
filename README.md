# Minimal term counts in quadratic mass-action systems

Nine configurations of invariant-line types have minimum component-monomial
counts **8, 8, 8, 8, 8, 8, 8, 6, 6** in the domain defined below. This package
contains the [article, draft v1](paper/paper-v1.pdf), its [TeX source](paper/paper-v1.tex),
exact witnesses, lower-bound certificates and two implementations for checking them.
It is intended for readers who want to inspect the proof or reproduce the exact checks.

The question comes from [Yeung and Erban, arXiv:2610.06130v1](https://arxiv.org/abs/2610.06130v1).
The geometric framework builds on [Lyagina (1951)](https://www.mathnet.ru/eng/rm6840)
and [Newton (1978)](https://doi.org/10.1137/1020007). We acknowledge those foundations;
the cited authors did not participate in or endorse this work.

## Results

Here `P2H` is one boundary line type. A comma separates invariant lines, and
the order of the lines does not matter. The best previous count in this table
means the best among Yeung and Erban's Table 1 and their Section 4 six-term
example, both in the pinned v1. It is not a claim about every earlier publication.

| Line-type multiset | Best count in that comparison | Minimum | Paper witness `(s,t,c,q,r,d)` |
|---|---:|---:|---|
| H | 9 | 8 | (0,0,3,1,0,2) |
| P2H | 8 | 8 | (1,0,0,2,0,1) |
| D,H | 10 | 8 | (0,0,9,4,0,2) |
| D,P2H | 8 | 8 | (2,0,0,3,0,1) |
| P1,P2,H | 9 | 8 | (0,0,56,25,0,12) |
| P1,H,H | 9 | 8 | (1,1,0,0,0,2) |
| P1,P2,P2H | 8 | 8 | (31,0,0,33,0,64) |
| P1,H,P2H | 6 | 6 | (1,0,0,0,0,1) |
| P1,P2H,P2H | 8 | 6 | (1,0,0,1,0,0) |

Five upper bounds decrease; four previously attained bounds are proved optimal
in the stated domain. Configurations are calculated from coefficients. The
printed labels in rows 14–16 and for the Section 4 example do not all match
those calculations; the article explains the comparison without imposing a
new historical numbering convention.

## Run the checks

Use Python 3.11 or later on macOS or Linux. Only the standard library is needed;
the checks use no network, solver or numerical root approximation. Run from
the package root:

```sh
python3 -B checkers/verify.py --output work/verification.json
python3 -B tools/regenerate.py --check
```

Success means exit code 0 and a JSON summary with `"status": "PASS"`.
The detailed output records nine paper witnesses checked by both implementations,
nine separately constructed witnesses plus eleven source and boundary controls,
two excluded controls, the lower-bound calculations, and seventeen rejected
corruptions. Optimized Python (`-O` or `PYTHONOPTIMIZE`) is refused because the
inherited checkers use assertions. The check process has a 120-second CPU limit.

The checks recompute the following:

- Chemical coefficient signs, equilibrium at `(1,1)`, zero Jacobian, term
  count, nonzero homogeneous resultant and a positive-rate reaction realization.
- All real invariant directions, including infinity; polynomial factorization,
  root multiplicities and exact signs. The research checker uses Sturm sequences
  and rational isolating intervals. The second implementation uses rational
  elimination, the rational-root theorem and exact quadratic-field arithmetic.
- All 49 supports with both components nonzero, independently all 64 supports
  including zero components, and the eight supports with fewer than eight terms.
  Symbolic identities reduce the admissible small supports to three families.
  The inequalities for arbitrary positive parameters are explained in the proof.
- Agreement of the paper's parameter table with the certificates, the result
  table in `claims.json`, file hashes, deterministic regeneration and rejection
  of corrupted witnesses and lower-bound inputs.

The second implementation is designed for the supplied examples: it requires
a rational root when the slope polynomial is cubic. It is not a general-purpose
classifier for arbitrary cubic coefficients. Neither implementation is formally
verified. Passing these programs is distinct from a proof-assistant theorem or
human peer review.

## Domain and definitions

There are two species, positive rate constants, reactant molecularity at most
two, and no bound on product molecularity. At a selected positive equilibrium,
the entire Jacobian is zero. Positive concentration scaling places the
equilibrium at `(1,1)`. In monomial order `(1,x,y,x^2,xy,y^2)`, the component
coefficient vectors are

```text
theta = (s,t,c,q,r,d) >= 0
a = (s, -2s-t, t, s+t+c, -t-2c, c)
b = (q, r, -2q-r, d, -r-2d, q+r+d)
M = number of nonzero entries of a + number of nonzero entries of b
```

The same monomial in both components counts twice. Counts refer to concentration
coordinates, before translating the equilibrium. With `u=x-1, v=y-1`, write
`P=A*u^2+B*u*v+C*v^2` and `Q=D*u^2+E*u*v+F*v^2`. We require the homogeneous
resultant `Res(P,Q)` to be nonzero, excluding zero components and common factors
over the complex projective line. The minima range over all real nonnegative
parameters satisfying this condition, not only rational witnesses or a finite grid.

For a finite invariant slope, `h(z)=z*P(1,z)-Q(1,z)=0` and
`rho=1-h'(z)/P(1,z)`. A simple root has type `P1` when `rho>1`, `P2` when
`0<rho<1`, `H` when `rho<0`, and `P2H` when `rho=0`. A double root has type `D`
and `rho=1`; this does not include triple roots. Infinity is checked in its
own chart. These are multisets of line types, not a claim that nine classes
are distinct under every notion of topological equivalence.

The lower bound uses the absence of cancellation in the signed coefficient
sums. Every nonzero component has at least three terms. Of the eight supports
with `M<8`, three have zero resultant; the other five reduce, by exchanging
species, to `P=a*u^2` with `Q=b*v^2`, `b*(u-v)^2` or `b*v*(v-u)`, where `a,b>0`.
Only the last two configurations in the table occur. Section 5 of the article
gives the complete argument, including the parameter inequalities.

## Scope and verification history

The result does not minimize the number of reactions. It does not cover zero
resultant, boundary equilibria, different molecularity constraints, or term
counts minimized under arbitrary affine changes of concentration coordinates.

The bounded literature comparison did not locate the same table of minima.
Full-text comparison with the classification domains of Lyagina and Newton
and a broader priority search remain incomplete. No priority claim is made.

A separate verification determined the same nine minima while withholding
the research-side minima, proofs, implementations and certificates. It shared
the problem specification, the six-parameter cone, term-counting rule,
line-type tests and published examples. Its two implementations were produced
within one verification effort; agreement is not independent human review.
This package includes its standard-library implementation. The earlier SymPy
run is part of the recorded research history, not an additional runtime
dependency or a third implementation rerun by this package.

See [claims.json](claims.json) for the claims and exclusions,
[VERIFICATION-RECORD.md](VERIFICATION-RECORD.md) for provenance and assurance,
and [PROVENANCE.json](PROVENANCE.json) for input hashes and adaptations.

## Regenerate and inspect

```sh
python3 -B tools/regenerate.py --output work/regenerated
python3 -B checkers/verify.py --witness-dir work/regenerated --output work/regenerated-check.json
```

Regeneration applies the positive time scales in `certificates/index.json` to
the normalized source certificates. It recreates the paper's integer-coefficient
witnesses, scaled reaction rates and resultants exactly. The lower-bound
certificate is supplied as a small input; `verify.py` derives and checks its
support counts and family identities. Detailed verdicts are regenerated by
the checker, not accepted as evidence from stored output. The larger exploratory
search and its unresolved cases are not required for these minima and are not
redistributed.

To rebuild the article, use an existing LaTeX installation with the packages
listed in its preamble. For example:

```sh
mkdir -p work/paper
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=work/paper paper/paper-v1.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=work/paper paper/paper-v1.tex
```

A rebuild may have a different PDF hash because of engine and metadata changes.
The included PDF is the unchanged draft v1. Building it is optional and is not
part of the standard-library certificate checks.

## Author, citation and license

Tatsuhiko Takimoto. Package version **1.0.0**; article **draft v1**, 24 October 2026.
Use [CITATION.cff](CITATION.cff) to cite the exact version. AI tools assisted the
research, software and writing; the author is responsible for the work.

Code is under Apache-2.0. Original prose and numerical data are under CC BY 4.0.
See [LICENSE](LICENSE), [Apache-2.0.txt](Apache-2.0.txt) and [CC-BY-4.0.txt](CC-BY-4.0.txt).
Cited papers are not redistributed, and these licenses do not relicense them.
