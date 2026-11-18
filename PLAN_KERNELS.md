# Precompiled-bookkeeping ("perf plan") kernels for the IMSRG(3) commutators

Opt-in with the environment variable `IMSRG_PERF_PLAN=1`. Without it, every term runs its existing code, unchanged.

## Headline numbers

i7-1355U laptop, 8–10 threads, OpenBLAS, 14N reference. "Existing code" is the same build without `IMSRG_PERF_PLAN`.

| | existing code | perf plan | speed-up |
|---|---|---|---|
| δ_C(14O → 14N), full IMSRG(3), emax 3, E3max 5 (887 commutators) | 2 h 03 m | 2 m 54 s | 42× |
| δ_C, full IMSRG(3), emax 2, E3max 6 (4990 commutators) | 80 min | 4.9 min | 16× |
| δ_C, full IMSRG(3), emax 2, E3max 4 (2054 commutators) | 13 min | 1.6 min | 8× |
| δ_C, full IMSRG(3), emax 3, E3max 9 (965 commutators) | not run (about 9 h, estimated from 965 × 34 s) | 17.3 min, peak 3.2 GB | ≈ 30× |
| one IMSRG(3) commutator, scalar, emax 3, E3max 9 | 34.4 s | 1.18 s | 29× |
| one IMSRG(3) commutator, isospin-changing Y, emax 3, E3max 9 | 66.5 s | 2.58 s | 26× |
| one IMSRG(3) commutator, scalar, emax 4, E3max 8 | out of memory (7 GB) | 5.8 s, peak 5.5 GB | runs |

All δ_C values agree with the existing code to all 16 digits, or to rounding (2e-12 relative at emax 2, E3max 4). The per-term breakdown is under **Speed** below.

## What changes

The maths of each term is unchanged: each plan kernel produces the same matrix elements as the term it replaces, to rounding (≤ 1e-12 relative in all tests). Only the bookkeeping changes:

- **Pair-antisymmetric product basis.** Per 3B channel, P = |(ab)J_ab c; J⟩ with a ≤ b, ordered by the third orbit. Stored kets map to P through two sparse matrices built once per model space:
  - R, the recoupling from `GetKetIndex_withRecoupling`;
  - A, the antisymmetriser: the direct, P(ik) and P(jk) terms with their 6j symbols and phases.
- **Unpack once per commutator.** For each 3B operator, W = R X Rᵀ and the stored-ket block X are computed once and shared by all terms of one `CommutatorScalarScalar` call, without checks: `perfplan_begin_commutator()` and `perfplan_end_commutator()` bracket the IMSRG(3) terms there. The half-unpacked R X is computed by the few terms that read it (comm331, comm233_pp_hh, comm332_ppph) while they run, so it is not resident during the memory-heavy terms (about 5 % of a commutator).
- **Reuse across commutators.** In a BCH series every nested commutator has the same Ω. When an operator is unpacked again (a new commutator, or a direct call of a single term), each channel's stored matrix is compared bitwise with the cached one, and only the channels that changed are recomputed. W depends on nothing else, so results are unchanged.
- **E3max < 3 emax.** The pph states of comm333_pph and the column pairs of comm233_ph carry no 3-body cut, but the back-transformation only reads entries whose two P states both exist. These two kernels, and their channel-changing (`_tz`) versions, compute and store only those entries, as blocks sorted by pair energy. The `_tz` versions also skip, in their 3B Pandya fills and GEMMs, the columns and intermediate states that cannot contribute under the cut (a product state P(g, a) needs e(g) + e_a <= E3max). At emax 4, E3max 8 that is about a quarter of the work and of the memory. For E3max = 3 emax nothing changes.
- **emax_3body < emax.** Only orbits with e ≤ emax_3body occur in 3B states, so everything indexed by 3B states is sized by that space, not by emax:
  - the pair list puts the pairs of the 3B space first (g < ng3), and the per-channel product-state and pph tables are ng3 × norb, read through `Plan::pidx`;
  - the dense 6j tables use the j and J ranges of the 3B space, except for the arguments that really run over all orbits (the summed orbit of comm223, the row orbit of comm232), with bounds-checked reads (`six6`);
  - the cross-coupled rows and summed pairs of comm233_ph, comm332_pphh and comm332_ppph, and the pph pair classes, only include pairs of the 3B space; 2B outputs of the [3B, 3B] terms skip kets outside it (they are exactly zero).

  For an isospin-changing Y (emax_3body 4, E3max 6) one commutator takes 2.0 s at emax 4, 3.8 s at emax 6 and 6.4 s at emax 8, at 3.0-4.5 GB.

  For emax_3body = emax all of this is the identity, and the results are bitwise unchanged. With emax_3body = 4, E3max 6, one scalar commutator takes 1.4 s at emax 4, 1.9 s at emax 6, 2.7 s at emax 8 and 7.3–12.7 s at emax 10 (on different days; 2.8–7.4 s of it IMSRG(2)-level terms); before, emax 6 ran out of memory.
- **Fills are table lookups and gathers into W**, blocked by the single orbits. Within a block the 6j weights depend only on the pair J values, so they are small tables rather than per-element `GetSixJ` calls.
- **GEMMs** for the contractions.
- **Back-transformation** into a dense W_Z per channel, then Z = Aᵀ W_Z A with column-contiguous sparse products.
- **Bulk 3B access.** `ThreeBodyME::GetChannelMatrix_pn` / `AddChannelMatrix_pn` read and write a whole channel block with one `ch_start` lookup, instead of one per element.
- **One-body embeddings** (`comm133`) are sparse matrices with a cached pattern, multiplied column by column.
- **Small products as loops.** In comm223 the per-block products (about 40 × 10 × 10) are column axpys rather than BLAS calls.
- **Channels are scheduled largest first** (dynamic scheduling).
- **Memory.** After the terms with large intermediates the freed memory is given back to the operating system (glibc `malloc_trim`), and comm333_pph's channel-changing kernel forms Z̄ in place, so later terms start from a lower resident size.
- **Single-precision W (optional).** Compiled with `-DIMSRG_PLAN_W_FLOAT`, W (and the off-diagonal W of channel-changing operators) is stored as `float`; everything else stays double. This halves the largest plan intermediate. Commutators then agree with the double build to about 1e-9 in their norms, and δ_C at emax 3, E3max 5 to 1.4e-8 relative. The default build is unchanged.
- **IMSRG(2), channel-changing Y.** For a scalar X and a J = 0 Y that changes isospin and/or parity, two IMSRG(2) terms use block products (in `Commutator.cc`):
  - comm122ss: `comm122ss_blocks` instead of the element-by-element `comm122ss_slower`, 4–9× faster at emax 8–10;
  - comm222_ph: `comm222_phss_blocks` instead of the general tensor route (`comm222_phst`). It forms the Pandya-transformed Z̄ blocks between cross-coupled channels that differ in Tz and/or parity with two GEMMs per block pair, as the scalar comm222_phss does within one channel; 2.4–3× faster at emax 10.
  - The isospin-changing IMSRG(2) commutator at emax 10 goes from 30.5 s to 10.0 s, and at emax 8 from 4.8 s to 1.8 s (same process, interleaved, 10 threads).

## How the kernels work (for readers new to the code)

Each source file has its own header with the term's formula and the steps in this notation; this section is the common part.

**Stored kets and the product basis.** A 3-body channel has fixed $(J, \pi, T_z)$. `ThreeBodyStorage_pn` stores matrix elements between antisymmetric kets $|k\rangle = |(ab)J_{ab}\,c;\,J\rangle$ with $a \le b \le c$. Many commutator terms are easier to write in the *pair-antisymmetric product basis*
$$
|p\rangle = |(ab)J_{ab}\,c;\,J\rangle, \qquad a \le b,\ c \text{ arbitrary},
$$
where only the pair $(ab)$ is antisymmetric and the third orbit $c$ is free. Two sparse matrices connect the two bases. $R$ expresses a product state in stored kets, $|p\rangle = \sum_k R_{pk}\,|k\rangle$ (the coefficients of `GetKetIndex_withRecoupling`), so the matrix of an operator $X$ in the product basis is
$$
W_X = R\,X\,R^{T}, \qquad (W_X)_{pq} = \langle p|X|q\rangle = \sum_{kl} R_{pk}\, X_{kl}\, R_{ql}.
$$
Going back, a result $W_Z$ built in the product basis is antisymmetrised by
$$
Z = A^{T}\, W_Z\, A, \qquad Z_{kl} = \sum_{pq} A_{pk}\,(W_Z)_{pq}\,A_{ql},
$$
where the columns of $A$ hold the direct, $P_{ik}$ and $P_{jk}$ terms of the $J$-coupled permutation operator $P(ij/k)^{J_1,J}$ (with the 6j symbols and phases of the original back-transformations), so that $A^T \cdots A$ implements $P(ij/k)\,P(lm/n)$.

**Unpacking.** For each operator in a commutator the kernels need, per channel, the stored block $X$, the half-unpacked $R\,X$ and $W_X$. These are computed once per `Commutator(X, Y)` and shared by all terms. Across commutators a channel is recomputed only if its stored block changed (bitwise check), which helps nested commutators $[\Omega,[\Omega,H]]$, where $\Omega$ is the same every time.

**Example: comm333_pph.** The term
$$
Z_{ijklmn} = -\tfrac12\, P(ij/k)\,P(lm/n) \sum_{abc} \left(n_a n_b \bar n_c + \bar n_a \bar n_b n_c\right) \left(X_{ijc,abn}\, Y_{abk,lmc} - Y_{ijc,abn}\, X_{abk,lmc}\right)
$$
becomes a matrix product after a particle–hole (Pandya) recoupling. In each particle–hole channel $J_{ph}$, with "external" states $e = (g_{12}, n)$ (a pair $g_{12} = (ij)J_{12}$ and a single orbit $n$) and internal states $\alpha = (g_{ab}, c)$ with a nonzero occupation factor:
$$
\bar X_{e\alpha} = \sum_{J'} (2J'+1) \begin{Bmatrix} j_n & J_{ph} & J_{12} \\ j_c & J' & J_{ab} \end{Bmatrix} (W_X^{J'})_{P(g_{12},c),\,P(g_{ab},n)},
\qquad
\bar Z = h_Y\, \bar X\, \bar Y^{T} - h_X h_Y\, (\ldots)^{T},
$$
with $\bar Y$ built the same way and weighted by the occupations, and $h = \pm 1$ for (anti-)Hermitian operators. Then, back in each 3-body channel,
$$
(W_Z)_{P(g_{12}, I_3),\,P(g_{45}, I_6)} = \sum_{J_{ph}} (2J_{ph}+1) \begin{Bmatrix} j_3 & J & J_{12} \\ j_6 & J_{ph} & J_{45} \end{Bmatrix} \bar Z^{J_{ph}}_{(g_{12}, I_6),\,(g_{45}, I_3)},
\qquad Z \mathrel{+}= A^T W_Z A .
$$
The other terms follow the same pattern: unpack, gather into dense matrices, multiply, transform back.

**Why it is faster.** The original code evaluates each recoupled 3-body matrix element with `GetME_pn` (hash lookups and recoupling sums per element). Here all bookkeeping that depends only on the model space ($R$, $A$, index tables, the 6j values that occur) is built once. Inside a block of fixed single orbits the 6j weights depend only on the pair $J$ values, so they come from small tables. The contractions are BLAS GEMMs.

**E3max < 3 emax.** A product state $|p\rangle$ exists only if $e_a + e_b + e_c \le E_{3\max}$ ($e = 2n + l$). The external states of comm333_pph and the column pairs of comm233_ph carry no such cut, but $\bar Z_{e_1 e_2}$ is only read if both product states it connects exist, e.g. $e(g_{12}) + e_{I_3} \le E_{3\max}$ and $e(g_{45}) + e_{I_6} \le E_{3\max}$. Sorting the states by pair energy makes the needed part of every block a contiguous corner, which is all that is computed and stored.

## Coverage

| term | scalar X, Y | scalar X, J = 0 Y that changes isospin and/or parity (e.g. the Fermi operator) |
|---|---|---|
| comm333_pph_hhpss | `_plan` | `_plan_tz` |
| comm233_phss | `_plan` | `_plan_tz` |
| comm223ss | `_plan` | `_plan_tz` |
| comm232ss | `_plan` | `_plan_tz` |
| comm332_pphhss | `_plan` | `_plan_tz` |
| comm332_ppph_hhhpss | `_plan` | `_plan_tz` |
| comm233_pp_hhss | `_plan` | `_plan_tz` |
| comm133ss | `_plan` | `_plan_tz` |
| comm331ss | `_plan` | `_plan_tz` |
| comm333_ppp_hhhss | `_plan` | existing code (already GEMM-based; a plan version measured no faster) |
| comm231ss | `_plan` | `_plan_tz` |
| comm132ss | `_plan` | `_plan_tz` |
| comm330ss | existing code | existing code |

The dispatch sits at the top of each term, after the existing check for reduced inputs. The plan kernels handle:
- the dE3max and occnat3 cuts: each unpacked operator and each 3B output is masked, as described below;
- the special IMSRG(3) flags of the scalar comm223 (`perturbative_triples`, `pert_trip_novvv`, `imsrg3_no_qqq`, `imsrg3_only_vvv`) and comm232 (`imsrg3_no_qqq`, `imsrg3_valence_2b`).

They fall back to the existing code for tensor operators (J-rank > 0), for `emax_3body < emax` in comm231 (its plan kernel loops over the full 2B space, and the existing code is faster there), and, in the channel-changing kernels, for the special flags (which the existing code applies only to the scalar terms).

**dE3max and occnat3 cuts.** A 3B state $|abc\rangle$ takes part only if $\sum_{x=a,b,c} |e_x - e_F| \le \texttt{dE3max}$ and $\prod_{x} n^{\rm nat}_x (1 - n^{\rm nat}_x) \ge \texttt{occnat3cut}$. The cut depends only on the three orbits, so it is the same for a product state and a stored ket and commutes with $R$ and $A$. The plan zeroes the rows and columns of failing kets in each unpacked operator, which carries over to $R X$ and $W$, and masks each 3B output block.

**Existing terms that apply the cuts incorrectly.** Under a non-default dE3max or occnat3 cut, four existing terms differ from the plan kernels because the existing code applies the cut incorrectly: `comm233_pp_hhss` (and `_tz`), `comm231ss`, `comm331ss`, and `comm232ss_slow` (the channel-changing path of comm232ss). The fixes are a separate bug-fix PR. Until it is merged, the ctest skips these four comparisons in its dE3max case.

## Validation

- **Reference tests.** `UnitTest::TestCommutators(X, Y)` with `IMSRG_PERF_PLAN=1` (N14, emax 1 and 2) compares each production term with `ReferenceImplementations`. It passes for:
  - X anti-Hermitian, Y Hermitian;
  - X and Y both anti-Hermitian;
  - an isospin-changing Y, both Hermitian and anti-Hermitian;
  - `TestCommutators_IsospinChanging()`.
- **Term by term**, plan vs existing code, with random 1+2+3B operators: ≤ 7e-13 relative at emax 1–3.
- **E3max < 3 emax.** All 12 scalar and 8 isospin-changing plan kernels, and the full commutator, vs existing code at (emax, E3max) = (2, 4), (2, 5), (3, 5), (3, 6), (3, 7): ≤ 5e-16 relative. The reference tests also pass at emax 1, E3max 2.
- **E3max-blocked `_tz` kernels** (comm233_ph, comm333_pph): the blocked path is taken at (emax, E3max) = (2, 4), (3, 5), (3, 6), (3, 7) for (T, P) = (1, 0), (0, 1), (1, 1); against the existing `_tz` terms <= 3.5e-16 relative.
- **emax_3body < emax.** Every kernel and the full commutator vs existing code at emax 3, emax_3body 2, for (T, P) = (0, 0), (1, 0), (0, 1): ≤ 7e-16 relative (in `test/UnitTestPerfPlan.py`). The changes for this case leave results bitwise identical for emax_3body = emax (checked on six model spaces, scalar, isovector and parity-odd Y, with and without an E3max cut).
- **Full `Commutator(X, Y)`**, with and without `IMSRG_PERF_PLAN`, scalar and isospin-changing Y: |dZ|/|Z| ≈ 3e-16 at emax 3, E3max 9.
- **Cuts.** All 23 kernels and the full commutator (scalar and isovector Y) vs existing code with the cut fixes applied, at dE3max = 2, dE3max = 3.5, occnat3cut = 3e-4 (with fractional natural occupations) and dE3max = 4 with occnat3cut = 1e-4, emax 2: |dZ| ≤ 1.5e-12 for |Z| ~ 10³.
- **Parity-changing Y** (T = 0 and T = 1, emax 1 and 2): every `_plan_tz` kernel vs its existing term, ≤ 4.3e-13.
  - The existing comm223ss aborts for such a Y (`std::out_of_range` from `map::at`), and so does the existing full `Commutator(X, Y)`. It reads 2B blocks as `GetMatrix(ch_bra, ch_ket)`, but only blocks with ch_bra ≤ ch_ket are stored.
  - With that lookup fixed (not part of this change), the existing comm223ss passes `UnitTest::Test_comm223ss` against `ReferenceImplementations`, and comm223ss_plan_tz agrees with it to ≤ 3e-13.
  - With comm223ss switched off on both sides, the full commutator agrees: |dZ| ≤ 8.5e-13 for |Z| ~ 2×10³.
- **comm122ss_blocks** vs `comm122ss_slower`: (T, P) = (1, 0), (0, 1), (1, 1), Hermitian and anti-Hermitian Y, emax 2–10: ≤ 2e-16 relative.
- **comm222_phss_blocks** vs the existing tensor route, the same cases, emax 2–10: ≤ 5.8e-16 relative; vs `ReferenceImplementations::comm222_phss` at emax 2: ≤ 2.7e-16. δ_C(¹⁴O → ¹⁴N) in IMSRG(2) at emax 4 is identical with the plan on and off.
- **comm331 and comm132 for emax_3body < emax** (plan kernels, scalar and `_tz`) vs existing code at (emax, emax_3body, E3max) = (4, 2, 6), (4, 3, 7), (5, 3, 8): ≤ 4e-15 relative.
- **Hunter-gatherer Magnus.** δ_C at emax 2, E3max 4, full IMSRG(3), with `SetHunterGatherer(true)`: identical to all 16 digits with the plan on and off (10.5 s instead of 141 s).
- **Flags.** comm223 and comm232 with `imsrg3_no_qqq`, `imsrg3_only_vvv`, `imsrg3_valence_2b`, and `perturbative_triples` with and without `pert_trip_novvv` (E_mp2 agrees to all printed digits).
- **End to end:** δ_C(14O → 14N), full IMSRG(3), 14N reference, Magnus flow plus transformation of the Fermi operator:
  - emax 3, E3max 9: agrees with the earlier plan run to 8e-13 relative (rounding-level reordering in comm223), same number of commutators;
  - emax 3, E3max 5: identical to all 16 digits to the run with the existing kernels (2 m 54 s instead of about 2 h).
  - emax 2, E3max 6: identical to all 16 digits, 4990 commutators (295 s instead of 4814 s).

## Speed (i7-1355U laptop, 8 threads, emax 3, E3max 9, 14N)

One full `Commutator(X, Y)` with random X (anti-Hermitian) and Y (Hermitian), 1+2+3B, plan on and off in the same process, best of 2 (`ClaudAssist/imsrg_design/termtable.py`). The one-time plan setup is excluded. Unpacking a 3B operator is counted in whichever term needs it first.

| term | scalar Y: existing | scalar Y: plan | speed-up | isovector Y: existing | isovector Y: plan | speed-up |
|---|---|---|---|---|---|---|
| comm233_phss | 13.26 s | 0.318 s | 42× | 25.02 s | 0.638 s | 39× |
| comm333_pph_hhpss | 11.25 s | 0.395 s | 28× | 15.45 s | 0.615 s | 25× |
| comm223ss | 5.18 s | 0.075 s | 69× | 12.40 s | 0.098 s | 126× |
| comm332_pphhss | 2.02 s | 0.074 s | 27× | 2.31 s | 0.134 s | 17× |
| comm232ss | 1.28 s | 0.103 s | 12× | 9.11 s | 0.133 s | 68× |
| comm332_ppph_hhhpss | 0.45 s | 0.033 s | 14× | 0.63 s | 0.026 s | 24× |
| comm233_pp_hhss | 0.32 s | 0.049 s | 6× | 0.45 s | 0.095 s | 5× |
| comm331ss | 0.16 s | 0.039 s | 4× | 0.22 s | 0.103 s | 2× |
| comm133ss | 0.13 s | 0.017 s | 7× | 0.20 s | 0.134 s | 1.5× |
| comm333_ppp_hhhss | 0.12 s | 0.026 s | 5× | 0.11 s | 0.072 s | existing code in both |
| comm231ss | 0.08 s | 0.015 s | 6× | 0.05 s | 0.014 s | 3× |
| comm132ss | 0.06 s | 0.002 s | 38× | 0.08 s | 0.001 s | 70× |
| comm330ss | 0.02 s | 0.015 s | existing code in both | 0 | 0 | |
| **whole commutator (wall)** | **34.4 s** | **1.18 s** | **29×** | **66.5 s** | **2.58 s** | **26×** |

Relative difference of the results: 2.8e-16 (scalar) and 3.8e-16 (isovector). The rows marked "existing code in both" don't use a plan kernel. Their times differ only because of the unpack accounting above.

Same measurement at smaller sizes, scalar Y:

| | existing code | perf plan | speed-up |
|---|---|---|---|
| emax 3, E3max 6 | 9.2 s | 0.42 s | 22× |
| emax 2, E3max 6 | 0.79 s | 0.04 s | 20× |

The full δ_C job at emax 3, E3max 9 (HF, Magnus flow to η < 1e-5, 965 commutators, ds_max 1, transformation of H and the Fermi operator) takes 17.3 min and peaks at 3.2 GB.

### Where the time goes in a δ_C calculation

δ_C(14O → 14N), full IMSRG(3), emax 3, E3max 5, step by step (`ClaudAssist/overnight/deltaC_steps.py`): 14N reference, HF, Coulomb, Magnus flow to η < 1e-5, ds_max 1, BCH threshold 1e-6.
- Existing-code column: the same job run earlier with `deltaC_run.py`, which records only the flow and the total.
- Plan column: 10 threads.

| step | existing code | perf plan | speed-up |
|---|---|---|---|
| model space, read interaction, Coulomb and Fermi operators, HF, normal ordering, 3B allocation | < 0.1 s | 0.07 s | |
| IMSRG flow: solve for Ω (860 commutators) | 7082 s | 147 s | 48× |
| transform H and the Fermi operator with Ω (27 commutators) | 269 s | 26 s | 10× |
| re-normal-order with respect to 16O, 2 × 2 diagonalisations, δ_C | < 0.1 s | 0.08 s | |
| **total** | **7351 s (2 h 03 m)** | **173.5 s** | **42×** |

- Almost all of the time is in the IMSRG(3) commutators, which is what the plan kernels replace.
- With the plan, the flow is spread over the terms. Over the whole run, from the profiler:

| term | time |
|---|---|
| comm233_phss | 51 s |
| comm333_pph_hhpss | 41 s |
| comm332_pphhss | 20 s |
| comm232ss | 19 s |
| comm223ss | 8.6 s |
| all 2B and 1B terms together | 4 s |

- The Fermi-operator transformation (12 isospin-changing commutators, 24 s) has the smaller speed-up (10×), because its commutators are more expensive per call.
- δ_C = 4.081948125356227e-02, identical to the earlier plan run, which agreed with the existing code to all 16 digits.

Most of the remaining time is in the comm333_pph and comm233_ph contractions (GEMM, close to the BLAS rate) and in their back-transformation gathers.

## E3max < 3 emax and larger model spaces (same laptop, 8 threads, scalar IMSRG(3) commutator, random X, Y)

| | existing code | perf plan |
|---|---|---|
| emax 3, E3max 6 | 9.2 s | 0.42 s |
| emax 4, E3max 8 | out of memory after 96 s (7 GB WSL limit) | 5.8 s, peak 5.5 GB |

Before the E3max-restricted blocks, comm233_ph alone ran out of memory here (7 GB), and comm333_pph took 10.4 s.

**Memory.** In units of $U = \sum_{\text{3B channels}} n_P^2$ doubles (0.3 GB at emax 3, E3max 9; 1.0 GB at emax 4, E3max 8; 3.2 GB at emax 4, E3max 10), one commutator holds about:
- $W$ for both operators: 2U (1U with `-DIMSRG_PLAN_W_FLOAT`);
- during comm333_pph, its stored $\bar Z$ corners: about 1U; during comm233_ph, its $Z3_{ph}$ blocks: about 0.7U;
- the plan tables (index tables, 6j tables, the comm233_ph column tables): 1–1.5 GB at emax 4, independent of the operators;
- plus per-thread work matrices.

Measured peak of one commutator at emax 4, E3max 8: 5.5 GB for a scalar Y and 7.1 GB for an isospin-changing Y (4.6 and 5.9 GB with single-precision W).

**Parallelism.** The plan kernels parallelize over channels (one channel per thread, largest first). With the largest channel at about 9 % of the work, one node can speed a term up by at most about 10× however many cores it has. On many-core nodes, letting the largest channels use several threads (nested or BLAS-level threading) would be the next step.

## What limits larger model spaces

With the plan, what limits a calculation depends on the size of the 3B space relative to emax. $N$ below is the number of single-particle states (140, 336, 660, 1144, 1820 at emax 4, 6, 8, 10, 12); the exponents are fits to the J-coupled channel sizes, not the m-scheme counts ($N^6$, $N^9$).

| regime | what dominates | numbers | possible next step |
|---|---|---|---|
| large emax, small 3B space (e.g. emax 10–12, emax_3body 4, E3max 6–8) | the IMSRG(2)-level scalar terms, which are existing code: comm222_pp_hh_221ss, comm222_phss, comm122ss | emax 10, emax_3body 4, E3max 6: 12.7 s per commutator, 7.4 s of it IMSRG(2)-level; the 3B terms (5.3 s) grow only slowly with emax. The pp-ladder work $\sum_{\rm ch} n_{\rm ch}^3$ grows as $N^{4.3}$ (×7.5 from emax 10 to 12), 2B storage as $N^{3}$ (0.29 → 1.18 GB per scalar operator) | faster Pandya fills and inverse transformations in comm222_phss (element-wise today); a Magnus flow holds ~10 operators, so at emax 12 the 2B storage alone is 10–15 GB |
| larger 3B space (emax_3body ≥ 5, E3max ≥ 10) | memory of $W$, then the comm333_pph and comm233_ph GEMMs | $W$ per operator, $(e_{3b}, E_{3\max})$: (4, 12) 4.8 GB, (5, 12) 28 GB, (6, 12) 66 GB, (6, 16) 317 GB; a commutator holds about 5 $W$-sized blocks. For E3max = 3 emax, $\sum n_P^2 \propto N^{5}$ and the GEMM work $\sum n_P^3 \propto N^{7.3}$ | single-precision $W$ (halves it); building $W$ per channel when needed; distributing channels over nodes (MPI) |
| many cores | the largest channels | the largest pph channel is ~9 % of the comm333_pph GEMM work and the largest 3B back-transformation channel ~7 %, which caps those stages at ~11× | nested or BLAS-level threading for the largest channels |
| long flows with split Magnus | memory of the stored $\Omega_i$, each a full 1+2+3B operator | grows with the number of splits | hunter-gatherer (keeps two Ω; truncated BCH product, δ_C changed by ~0.5 % at emax 2, E3max 4), or imsrg++'s scratch directory for Ω |

Two machine notes. They aren't part of this change, but they matter on similar laptops:
- **OpenBLAS core detection.** OpenBLAS can misdetect hybrid Intel CPUs under WSL as `Prescott` (SSE3 kernels). `OPENBLAS_CORETYPE=Haswell` made `dgemm` 3.6× faster there.
- **Build flags.** `-march=native` and `-DARMA_NO_DEBUG` gave no measurable change at 8 threads.

## Building on this

### Testing a change

`test/UnitTestPerfPlan.py` (registered with ctest, about 40 s at emax 2) compares every `_plan` and `_plan_tz` kernel with its existing term:
- for a scalar Y and for J = 0 Y with (T, P) = (1, 0), (0, 1), (1, 1);
- for E3max < 3 emax and for dE3max = 2.

It then compares the full `Commutator(X, Y)` with the plan on and off, and runs the `ReferenceImplementations` tests with the plan on. To check a new kernel, add its term to `PLAN_TERMS` or `TZ_TERMS`.

For timings, `test/BenchmarkPerfPlan.py` times one full commutator with the existing code and with the plan, in the same process (best of `--reps`), term by term, and prints a markdown table (options: `--emax`, `--e3max`, `--emax3b`, `--isovector`, `--ref`, `--plan-only`). Timings from separate processes on a laptop can vary by 15–30 %.

### The shared pieces (`PerfPlan.inc`, namespace `perfplan`)

| piece | what it is |
|---|---|
| `get(X, Z)` → `Plan&` | Everything that depends only on the model space and the occupations: the product basis (`Pst`, `PIdx`), `R` and `A` per channel, the pph channels of comm333_pph, the 6j tables. Built on first use; rebuilt when the model space or the occupations change. |
| `unpack(op, P)`, `unpack_k`, `unpack_rx` | Per 3B channel: W = R X Rᵀ (type `wmat`: double, or float with `-DIMSRG_PLAN_W_FLOAT`) and the stored block X, for an operator whose 3B blocks are diagonal in the channel, cached (see "Reuse across commutators" above); `unpack_rx` returns R X, computed from the cached X when called. |
| `unpack_off(op, P)` → `WOffCache` (`PerfPlan333tz.inc`) | The same for an operator with off-diagonal 3B blocks (cb ≠ ck). It loops over whatever channel pairs the operator stores, and `idx[cb·nch + ck]` gives the block. |
| `AtWA(Ab, WZ, Ak)` | Z(cb, ck) = Aᵀ W_Z A, with column-contiguous sparse products. |
| `cut3(ms)`, `mask3`, `pass3` | The dE3max / occnat3 cut on 3B states, applied to unpacked inputs and to outputs. |
| `ThreeBodyME::GetChannelMatrix_pn` / `AddChannelMatrix_pn` | Read or add a whole 3B channel block at once. |
| `perfplan_begin_commutator()` / `perfplan_end_commutator()` (`Commutator.cc`) | Bracket one `CommutatorScalarScalar` call, so that X and Y are unpacked once for all its terms. |

### Adding a kernel

1. Start from the closest `PerfPlan*.inc`. The `_tz` files are the templates for operators with off-diagonal channel blocks. Give it the same kind of header: the formula as written in `IMSRG3Commutators.cc`, then the steps in the notation above.
2. First line of the kernel: a guard that calls the existing term for every case the kernel doesn't cover. The plan never has to be complete.
3. Unpack with `get` / `unpack` / `unpack_off`. Loop over channels with `#pragma omp parallel for schedule(dynamic, 1)`, largest first. Mask with `mask3` and write with `AddChannelMatrix_pn`.
4. Declare it in `IMSRG3Commutators.hh`, `#include` it at the end of `IMSRG3Commutators.cc`, and add a dispatch block at the top of the existing term, after its reduced-input check (copy one of the existing blocks). Bind it in `pyIMSRG.cc`.
5. Add it to `test/UnitTestPerfPlan.py`.

### What is not covered, and where to start

- **Tensor operators (J-rank > 0).** These still use the existing tensor terms. The product basis, R and A are per channel and stay as they are, and `unpack_off` already handles any pair of stored channel blocks. What changes is the angular-momentum algebra in the fills and back-transformations: an extra recoupling with the operator rank, and the reduced-matrix-element convention. I'd expect the `_tz` kernels to be the right starting point, since they already work with bra and ket channels that differ.
- **comm330ss** (0-body output, cheap), and comm231 for `emax_3body < emax`: existing code.
- **Special flags** (`imsrg3_no_qqq` etc.) in the channel-changing kernels: existing code, as in the existing channel-changing terms.
- **Further IMSRG(2) speed-ups for a channel-changing Y (optional).** This case is covered: comm122 and comm222_ph have block kernels, and the remaining large term, comm222_pp_hh_221ss, is existing code that is already GEMM-based (about 3 s of the 10 s at emax 10). If more speed is needed here, the next step would be in `comm222_phss_blocks`, whose time is split about evenly between the element-wise Pandya fills of Ȳ and the element-wise inverse transformation.
- **Memory.** W is cached for both operators (5.2 GB per operator at emax 4, E3max 12 in double precision). Beyond single precision, the next steps would be distributing the channels over several nodes (MPI), or building W per channel when needed (this needs the kernels reorganised by channel, since most of them read W of many channels at once).
- **Many-core nodes.** Parallelism is over channels, and the largest pph channel is about 9 % of the comm333_pph GEMM work, which caps that stage at about 11×. Threading inside the largest GEMMs (nested OpenMP or threaded BLAS for those calls) is the next thing to try. On the laptop it was slower, so it isn't included.
- **Making it the default.** After more production use, `IMSRG_PERF_PLAN` could become opt-out. The existing terms would stay as the reference implementation that the test compares against.

## Files

- `src/Commutator/PerfPlan*.inc`: one file per kernel, included at the end of `IMSRG3Commutators.cc`.
- `src/Commutator/IMSRG3Commutators.cc/.hh`: dispatch at the top of each term, and declarations.
- `src/Commutator/Commutator.cc`: `perfplan_begin_commutator()` / `perfplan_end_commutator()` around the IMSRG(3) terms, `comm122ss_blocks` and `comm222_phss_blocks`.
- `src/ThreeBodyME/*`: bulk channel accessors.
- `src/pyIMSRG.cc`: Python bindings for the `*_plan` / `*_plan_tz` functions.
- `test/UnitTestPerfPlan.py`: the test above (registered in `CMakeLists.txt` and `src/CMakeLists.txt`).
- `test/BenchmarkPerfPlan.py`: the benchmark above (not part of ctest).
