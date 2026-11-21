#!/usr/bin/env python3
# Benchmark of the perf-plan kernels (IMSRG_PERF_PLAN, see PLAN_KERNELS.md) on this machine: one full IMSRG(3)
# Commutator(X, Y) with random 1+2+3B operators, existing code vs plan, in the same process (best of --reps), with the
# time of every IMSRG(3) term. Prints a markdown table. Not part of ctest.
#
#   python3 test/BenchmarkPerfPlan.py                          # emax 3, E3max 9, scalar Y
#   python3 test/BenchmarkPerfPlan.py --emax 4 --e3max 8       # larger space (the existing code needs > 7 GB here)
#   python3 test/BenchmarkPerfPlan.py --isovector              # isospin-changing Y (as for the Fermi operator)
#   python3 test/BenchmarkPerfPlan.py --emax 8 --emax3b 4 --e3max 6 --plan-only
#
# Threads: set OMP_NUM_THREADS; with OpenBLAS, OPENBLAS_NUM_THREADS=1 (BLAS is called inside OpenMP loops).

import argparse, os, time
import pyIMSRG
from pyIMSRG import ModelSpace, UnitTest, Operator, Commutator, IMSRGProfiler

ap = argparse.ArgumentParser()
ap.add_argument('--emax', type=int, default=3)
ap.add_argument('--e3max', type=int, default=-1, help='3B cut e1+e2+e3 (default 3 emax_3body)')
ap.add_argument('--emax3b', type=int, default=-1, help='emax of the 3B space (default emax)')
ap.add_argument('--isovector', action='store_true', help='Y with isospin rank 1 (scalar X)')
ap.add_argument('--ref', default='O16', help='reference (default O16)')
ap.add_argument('--valence', default='', help='valence space (default: same as --ref)')
ap.add_argument('--reps', type=int, default=2)
ap.add_argument('--plan-only', action='store_true', help='skip the existing code (e.g. when it does not fit in memory)')
a = ap.parse_args()

val = a.valence or a.ref
ms = ModelSpace(a.emax, a.emax3b, a.ref, val) if a.emax3b >= 0 else ModelSpace(a.emax, a.ref, val)
if a.e3max >= 0:
    ms.SetE3max(a.e3max)
ms.PreCalculateSixJ()
Commutator.SetUseIMSRG3(True)
ut = UnitTest(ms)
ut.SetRandomSeed(5)
TY = 1 if a.isovector else 0
X = ut.RandomOp(ms, 0, 0, 0, 3, -1)
Y = ut.RandomOp(ms, 0, TY, 0, 3, +1)
if Y.IsReduced():
    Y.MakeNotReduced()

TERMS = ['comm233_phss', 'comm333_pph_hhpss', 'comm223ss', 'comm332_pphhss', 'comm232ss', 'comm332_ppph_hhhpss',
         'comm233_pp_hhss', 'comm331ss', 'comm133ss', 'comm333_ppp_hhhss', 'comm231ss', 'comm132ss', 'comm330ss']
prof = IMSRGProfiler()


def run(plan):
    if plan:
        os.environ['IMSRG_PERF_PLAN'] = '1'
    else:
        os.environ.pop('IMSRG_PERF_PLAN', None)
    before = {t: prof.GetTimer(t) for t in TERMS}
    t0 = time.perf_counter()
    Z = Commutator.Commutator(X, Y)
    wall = time.perf_counter() - t0
    return wall, {t: prof.GetTimer(t) - before[t] for t in TERMS}, Z


modes = [True] if a.plan_only else [True, False]
for plan in modes:   # warm-up, not timed (one-time plan set-up; creates the profiler timers)
    run(plan)
best = {}
for r in range(a.reps):
    for plan in modes:
        w, d, Z = run(plan)
        if plan not in best or w < best[plan][0]:
            best[plan] = (w, d, Z)

threads = os.environ.get('OMP_NUM_THREADS', '?')
print(f'\nemax {a.emax}, emax_3body {a.emax3b if a.emax3b >= 0 else a.emax}, E3max {ms.GetE3max()}, {"isospin-changing" if TY else "scalar"} Y, '
      f'OMP_NUM_THREADS={threads}, best of {a.reps}')
if a.plan_only:
    wp, dp, Zp = best[True]
    print('\n| term | perf plan |\n|---|---|')
    for t in sorted(TERMS, key=lambda t: -dp[t]):
        print(f'| {t} | {dp[t]:.3f} s |')
    print(f'| **whole commutator (wall)** | **{wp:.2f} s** |')
else:
    wp, dp, Zp = best[True]
    wo, do, Zo = best[False]
    print(f'relative difference of the results: {(Zp - Zo).Norm() / Zo.Norm():.1e}')
    print('\n| term | existing code | perf plan | speed-up |\n|---|---|---|---|')
    for t in sorted(TERMS, key=lambda t: -do[t]):
        sp = f'{do[t] / dp[t]:.1f}x' if dp[t] > 1e-4 else ''
        print(f'| {t} | {do[t]:.3f} s | {dp[t]:.3f} s | {sp} |')
    print(f'| **whole commutator (wall)** | **{wo:.2f} s** | **{wp:.2f} s** | **{wo / wp:.1f}x** |')
