# Entry point for a computational reader

Read README.md, claims.json and VERIFICATION-RECORD.md. The article is
paper/paper-v1.tex; paper/paper-v1.pdf is its unchanged draft v1.

Run `python3 -B checkers/verify.py` and `python3 -B tools/regenerate.py --check`.
Python 3.11+, standard library, macOS or Linux. No network is needed. Never run
the checkers under optimized Python: their assertions are part of the checks.

Keep the evidence types separate:

- `exact-computation`: exact rational, polynomial and quadratic-field checks;
  the programs have not been formally verified.
- `independent-implementation`: agreement with a separately developed
  verification that shared the specification and mathematical preparation.
  Packaging adaptations were made after that verification, with both sets of
  results available. Do not call package assembly a new blind verification.
- `prose`: the proof over arbitrary real parameters and its inequalities.
  Symbolic identities and support enumeration support it, but do not turn it
  into a proof-assistant theorem.

Do not claim a minimum number of reactions, priority, human peer review,
historical-domain equivalence, or a complete general classifier. Preserve the
explicit nonzero-resultant domain and credit to Yeung–Erban, Lyagina and Newton.
Numeric source row labels locate examples; line-type multisets identify the results.

The package contains no publication identifier or remote configuration. A local
commit, metadata file or passing CI is not evidence of publication. External
publication and communication require the author's instruction.

After an intentional file change, inspect it before running
`python3 -B tools/update_manifest.py`, then rerun the checks. Do not refresh hashes
to conceal an unexplained mismatch. Keep generated output in ignored `work/`.
