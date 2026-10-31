# Independent audit and presentation recommendations

Audit conclusion: previously reported errors, necessary Schmidt ranks, and constructed bond dimensions passed independent validation; the numerical tables require no changes.

The audit covers 15 velocity arrays, direct SVD at 90 scale cuts, and full-grid reconstruction from 45 MPS files reread from disk. It does not call the previous tensor-permutation, rank-selection, or reconstruction helpers. Matrices are independently assembled from physical-coordinate coarse/fine blocks. Direct SVD also checks the large cuts previously treated through Gram eigenvalues.

Every velocity-array SHA-256 matches its export record. Regenerated 15- and 27-qubit fields match pointwise. The maximum normalized difference between old and new spectra is 9.597866e-14; the maximum entropy difference is 5.595524e-14 bits; the maximum difference between reconstruction-error measurements is 3.469447e-17.

The original findCritChi3D.m was also executed on a 64^3 grid. Its spectrum difference is 3.662037e-15, and maximum single-cut ranks under the original joint-error convention at the three tolerances are [475.0, 504.0, 511.0].

## Rank bounds depend on tensorization

For an MPS with physical dimensions $d_1,\ldots,d_L$, the cut-rank bound is $r_\ell\leq\min(\prod_{j\leq\ell}d_j,\prod_{j>\ell}d_j)$. This experiment uses $N=2^L$ with one $d=8$ site per scale, giving $\chi_{\max}=8^{\lfloor L/2\rfloor}$. For even $L$, this equals $N^{3/2}$; for odd $L$, it equals $N^{3/2}/\sqrt{8}$.

With $3L$ binary sites, the bound becomes $2^{\lfloor3L/2\rfloor}$. Treating the $N\times N\times N$ array as three $N$-dimensional sites gives maximum rank $N$, but the central core contains $Nr_1r_2$ scalars. A smaller numerical chi therefore need not mean fewer parameters. Chi values from different tensorizations are not directly comparable.

The audit also examines both cuts in the conventional three-site tensorization with physical dimension N. For the 512³ field at 1% componentwise error, necessary ranks are 510-511 against a maximum of 512. Absolute chi changes with tensorization, while nearly full ranks persist in this additional comparison; all orderings have not been exhausted.

| n | Grid | Octal bound | Necessary rank at 1% componentwise error | Verified construction | Construction/bound |
|---:|---:|---:|---:|---:|---:|
| 15 | 32³ | 64 | 64 | 64 | 100.0000% |
| 18 | 64³ | 512 | 475 | 475 | 92.7734% |
| 21 | 128³ | 512 | 512 | 512 | 100.0000% |
| 24 | 256³ | 4096 | 3427 | 3427 | 83.6670% |
| 27 | 512³ | 4096 | 3994 | 4027 | 98.3154% |

## Why the ranks are large

The direct numerical cause is a long Schmidt-spectrum tail at multiple scale cuts. A 1% relative L2 error requires retaining at least 99.99% of the squared norm, not 99%.

| Maximum bond cap for the 512³ grid | Unavoidable relative-error lower bound in at least one component |
|---:|---:|
| 256 | 49.161732% |
| 512 | 40.170358% |
| 1024 | 29.585927% |
| 2048 | 16.813996% |
| 3072 | 7.926782% |
| 3994 | 0.999218% |

These necessary bounds follow from optimal low-rank matrix approximation at original-field cuts, independently of TT-SVD tuning. Long tails and high necessary ranks are numerical properties of the current scale ordering.

Phase structure, derivatives, products, and division by density in nonlinear velocity decoding do not guarantee low MPS rank. These are plausible mechanisms; their individual contributions have not been established through ablation. Physical-space smoothness and coarse/fine-scale separability are different conditions.

## Error conventions

The main table requires each velocity component to satisfy its own relative L2 tolerance. The original findCritChi3D.m sums discarded squared norms over the three components, divides by their combined squared norm, and tests each cut separately. The main-table condition is stricter. This is not a numerical error, but comparisons must use matching definitions. Also, the largest necessary single-cut rank bounds the full-MPS bond dimension from below; it need not equal the smallest uniform cap satisfying the full reconstruction error.

At 1% error, the original joint convention gives maximum single-cut ranks 64, 475, 512, 3423, and 3990. The componentwise convention gives lower bounds 64, 475, 512, 3427, and 3994. The conclusion does not depend on this small difference. Both 1% criteria refer to L2 error, so the original code receives a squared-error budget of 1e-4.

## Supported conclusions

For these realizations, manuscript parameters, and NCS octal scale ordering, the velocity fields show little low-rank MPS compression. Constructed gauge-reduced NVPS is 97.03%-100% of the dense count; stored dense cores require 1.216-2.025 times the dense storage. Necessary-rank storage bounds also exceed dense arrays, excluding conservative uniform-cap selection as the sole cause of the storage result.

This provides numerical evidence for assessing the inference that these velocity fields naturally admit low-bond MPS representations. It does not directly test spinor-state MPS compression or compare complete quantum and classical evolution with matched outputs. Classical computers also store compact generating circuits; high MPS rank alone establishes neither quantum advantage nor improved physical flow quality. Not every bit ordering or tensor-network structure was tested.

## Presentation recommendations

The candidate diagnostic figure is figures/mps_velocity_20261029/mps_recheck_diagnostics.png/pdf.

1. Plot chi/chi_max, showing necessary bounds and verified constructions, to distinguish the chain-parity staircase from physical scaling.
2. Include the unavoidable-error versus bond-cap curve for the 512³ field to explain why small caps fail. Plot measured MPS reconstruction errors as separate symbols; do not label the lower-bound curve as measured MPS error.
3. Show gauge-reduced NVPS and actual core storage separately, normalized by dense velocity storage. Do not use the circuit-description gap between different output tasks as evidence of quantum advantage.
4. Place Schmidt spectra, component entropies, the two additional tolerances, the original joint convention, and complete validation records in supplementary material.
5. Label every result as one realization and disclose grid-dependent settings in the caption. Do not infer a new asymptotic exponent from five sizes.

Data files are upper_bounds_and_error_conventions.csv, unavoidable_error_curves.csv, independent_validation.json, and all recomputed spectra under direct_spectra/.

Source: [Gourianov et al., Eq. (5) and Methods](https://arxiv.org/html/2106.05782). Original public MATLAB files are preserved in original/.
