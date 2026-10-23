# Verification record

Package 1.0.0, prepared on 24 October 2026. This record separates the historical
research checks from what the portable package executes.

## Fixed inputs

The article is copied byte for byte from research commit
`0859caa140e062225a6f52d3550f04fa0c95c8a0`, integrated in input commit
`bc151e1d1dba992467996658c6b7cf7f5cee947f`. The included PDF has seven pages and SHA-256
`b9f687157ce100973944c52ef3df7069abf73e934baa17df02b1b2b6e5b30f97`.
Its author/date line still says draft v1. The package version is 1.0.0.

`PROVENANCE.json` identifies the source files, their hashes and packaging
adaptations. `manifest.json` pins every distributed file except itself.
Normalized research certificates, the lower-bound certificate and the separate
verification's input examples are unchanged. The nine paper certificates apply
positive time scales to the normalized witnesses. Coefficients, the slope
polynomial and rates scale linearly; the homogeneous resultant scales to the
fourth power. Root locations, monic factors and line types stay unchanged.

## Evidence carried by the package

The research checker reconstructs the cone, equilibrium and zero Jacobian;
checks the coefficient count and resultant; factors the invariant-slope
polynomial; verifies multiplicities, complete real-root counts and isolating
interval signs; treats the vertical direction; and reconstructs the ODE from
positive-rate reactions. Its lower-bound checker enumerates the 49 nonzero
component supports and verifies sparse-polynomial identities for the eight
supports below eight terms. The passage from those identities to inequalities
for arbitrary positive real parameters is the displayed proof in Section 5.

The separate standard-library implementation solves the six linear equilibrium
and Jacobian constraints by rational elimination, classifies the supplied
examples using exact quadratic fields, enumerates all 64 zero-or-positive
supports, and checks three generic families in a Laurent polynomial ring.
Its nine own witnesses, eleven published-example or boundary controls, two
resultant-zero exclusions and five corrupted inputs are rerun here. The
portable entry point additionally checks all nine paper witnesses through that
implementation, cross-checks the small supports, and exercises twelve
corruptions of research certificates. Optimized Python is explicitly refused.

The proof is over all real nonnegative parameters in the defined domain. It
does not depend on treating a finite set of samples as exhaustive. The broader
exploratory feasibility search included unresolved larger-support cases; those
cases are not premises of the minima and are omitted from this package.

## Independence and review limits

The separate verification originally had the specification, not the research
minima, implementations, certificates or proof. Shared preparation included the
six-parameter cone, fixed-sign coefficient counting, invariant-line criteria,
and published examples. It was not blind to those examples. Its SymPy and
standard-library implementations agreed on the nine minima and their controls.
This package reruns the standard-library implementation, not the earlier SymPy
run. Both derivations were available during packaging; adaptation and cross-checks
here are not a new independent derivation.

The mathematical functions of the three inherited checker modules are retained.
Packaging changes add optimization guards, relocate input/import/output paths,
and convert Linux resource units in the optional standalone second checker.
The package entry point, manifest check and regeneration tool are new integration
code. That integration code does not serve as a third mathematical implementation.

Draft v1 had two AI fresh-reader review rounds; the recorded final round accepted
the revised draft. This is not human peer review. Package assembly did not launch
another reviewer. No proof-assistant formalization or verified compiler is claimed.

## Literature boundary

The comparison is with Yeung and Erban, arXiv:2610.06130v1, Table 1 and Section 4,
using explicit coefficients instead of inconsistent printed labels. The five
decreased bounds include the six-term realization of `P1,P2H,P2H`; the known
Section 4 six-term example already realizes `P1,H,P2H`. The article also records
an admissible equality case for the discrepancy between Lemma 6(ii) and its
proof equation (3.43). No strict version of that inequality is used here.

The bounded comparison also examined the related problems cited in the article.
It did not establish worldwide priority. Lyagina's and Newton's full historical
classification domains have not been compared in full text with the domain here.
The line-type multisets, nonzero resultant, unrestricted products and component
term-count definition must stay attached to the result.

## Reproduction and integrity

`python3 -B checkers/verify.py --output work/verification.json` recomputes the
claims described above and records the interpreter, platform, elapsed time and
peak process memory. `python3 -B tools/regenerate.py --check` compares the nine
regenerated certificates byte for byte. Neither command downloads inputs.
The CI workflow runs both commands on Python 3.11 and 3.14 on Linux. A workflow
definition is not evidence that hosted CI has run.

The package is prepared locally. Publishing, repository creation, DOI allocation
and archival deposit are separate actions. Citation and deposit metadata are
provided without inventing identifiers or release URLs.
