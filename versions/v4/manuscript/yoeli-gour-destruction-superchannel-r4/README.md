# Fixed points of superchannels for dephasing-covariant channels

Extending the two-sided destruction construction of Saxena, Chitambar and
Gour (2019), this note answers the fixed-superchannel question in Yoeli and Gour,
*Instability as a Quantum Resource*, arXiv:2602.18323v1, p. 5,
Open directions. For complete dephasing the answer is positive in dimensions
1 and 2 and negative in every dimension at least 3. For complete
depolarization it is negative in every dimension at least 2.

## Mathematical contents

- Theorem 1: a completely positive Choi action fixing all dephasing-covariant
  channels in dimension at least 3 is the identity.
- Corollary 3: an explicit physical projection works for any two-sector
  pinching, including qubit complete dephasing.
- Corollary 4: the corresponding rigidity holds for all uniformity-preserving
  channels in every dimension at least 2 and for conditional nonuniformity.
  With an unrestricted side system, the dephasing classification is still
  determined by the dimension of the dephased factor.
- The two-sided destruction construction fixes the smaller two-sided-invariant
  class. The classical-channel case is credited to Saxena, Chitambar and Gour
  (2019 preprint; 2020 publication); the uniform-replacer example is credited
  to Yoeli and Gour (2026).
- Corollary 5: the natural linear projection onto covariant maps fails
  complete positivity, with an explicit negative Choi principal determinant.
- Corollary 6: an affine density-matrix set need not admit a self-adjoint
  resource-destroying map, giving a counterexample to Theorem 9.2 as stated
  in arXiv:2402.05474v2, printed p. 442.

Each result has a written proof in `note.tex`. The source coverage table
distinguishes answered questions, consequences not separately posed, cited
earlier results, the counterexample and items not covered.

## Formal and finite coverage

The pinned Lean source is
https://github.com/the-omega-institute/trureturing/blob/a252cc46f0e145418a978217b8abd41a73fd820d/D5/S3/QuantumChannels/YoeliGourDestructionSuperchannelRefutation.lean

Its `result` declaration negates the universal exact-fixed-CPTP-set assertion
using qutrit complete dephasing. Its `fixes_dio_identity` helper establishes
the qutrit rigidity statement. The dimension-uniform extensions and the
remaining corollaries have the written proofs specified above.

`verify.py` supplies exact finite checks over the Gaussian rationals:
the qutrit spanning/overlap family, the complete qubit fixed-space mask,
the seven-vector qubit uniformity family, two-sector pinching coefficient
masks for selected sector ranks in dimensions 2 through 4, the two negative
determinants, a full-rank replacer identity and the qubit group-covariance
distinction. Their exact scope and expected results are in `verification.json`.
Finite checks certify these witnesses rather than all parameter values.

## Build and check

Compile with Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a standard LaTeX installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded; there are no external figures or data files.
With Python 3.10 or later, run the finite checks:

```sh
python3 -m pip install -r requirements.txt
python3 verify.py > actual-verification.json
diff -u verification.json actual-verification.json
```

`metadata.json` records the scientific scope, sources, creators and license.
`source-manifest.json` gives the exact sizes and SHA-256 values of the
package members other than the manifest itself.

## Creators and license

1. Haobo Ma — ChronoAI Pte Ltd, Singapore.
2. Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

CC BY 4.0 (Creative Commons Attribution 4.0 International).
