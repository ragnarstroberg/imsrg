#!/usr/bin/env python3
# Perf-plan kernels (IMSRG_PERF_PLAN=1, see PLAN_KERNELS.md) against the existing code they replace.
#
# Each *_plan / *_plan_tz kernel is compared with its existing term on random 1+2+3B operators, for
#   - scalar Y, and J = 0 Y that changes isospin, parity or both;
#   - E3max < 3 emax, emax_3body < emax, and a dE3max cut;
# then the IMSRG(2) block kernels for a channel-changing Y (comm122ss, comm222_phss), the full Commutator(X, Y) with
# the plan on and off, and the ReferenceImplementations tests with the plan on.
# Runs in about a minute. To check a new kernel, add it to PLAN_TERMS / TZ_TERMS.

import os
import pyIMSRG
from pyIMSRG import ModelSpace, UnitTest, Operator, Commutator

PLAN_TERMS = ['comm333_pph_hhpss', 'comm233_phss', 'comm223ss', 'comm332_pphhss', 'comm232ss', 'comm332_ppph_hhhpss',
              'comm233_pp_hhss', 'comm133ss', 'comm331ss', 'comm132ss', 'comm231ss', 'comm333_ppp_hhhss']
TZ_TERMS = [t for t in PLAN_TERMS if t != 'comm333_ppp_hhhss']   # comm333_ppp keeps its existing channel-changing code
TOL = 1e-10   # on |Z_plan - Z_existing| / |Z_existing|

Commutator.SetUseIMSRG3(True)
failures = []


def plan(on):
    if on:
        os.environ['IMSRG_PERF_PLAN'] = '1'
    else:
        os.environ.pop('IMSRG_PERF_PLAN', None)


def make_ms(emax, e3max=None, de3max=None, emax3b=None):
    ms = ModelSpace(emax, emax3b, 'N14', 'p-shell') if emax3b is not None else ModelSpace(emax, 'N14', 'p-shell')
    if e3max is not None:
        ms.SetE3max(e3max)
    if de3max is not None:
        ms.SetdE3max(de3max)
    ms.PreCalculateSixJ()
    return ms


def random_ops(ms, ty, py, seed=5):
    ut = UnitTest(ms)
    ut.SetRandomSeed(seed)
    X = ut.RandomOp(ms, 0, 0, 0, 3, -1)    # scalar, anti-Hermitian
    Y = ut.RandomOp(ms, 0, ty, py, 3, +1)  # J = 0, Hermitian
    if Y.IsReduced():
        Y.MakeNotReduced()
    return ut, X, Y


def new_Z(ms, ty, py):
    Z = Operator(ms, 0, ty, py, 3)
    Z.ThreeBody.SetMode('pn')
    Z.SetHermitian()   # [anti-Hermitian, Hermitian] is Hermitian
    if Z.IsReduced():
        Z.MakeNotReduced()
    return Z


def check(label, Z_plan, Z_ref):
    rel = (Z_plan - Z_ref).Norm() / max(Z_ref.Norm(), 1e-300)
    ok = rel <= TOL
    print(f'{"ok  " if ok else "FAIL"} {label}: |Z| {Z_ref.Norm():.4e}  rel. diff {rel:.1e}', flush=True)
    if not ok:
        failures.append(label)


# Under a dE3max / occnat3 cut, the existing comm233_pp_hhss, comm231ss, comm331ss and the channel-changing path of
# comm232ss apply the cut incorrectly (fixed in a separate bug-fix PR); the plan kernels apply it as the other terms do.
CUT_BUGGY_EXISTING = {'comm233_pp_hhss', 'comm231ss', 'comm331ss', 'comm232ss'}


def kernels(tag, ms, ty, py, cut=False):
    ut, X, Y = random_ops(ms, ty, py)
    suffix = '_plan' if (ty, py) == (0, 0) else '_plan_tz'
    for term in (PLAN_TERMS if suffix == '_plan' else TZ_TERMS):
        if py == 1 and term == 'comm223ss':
            continue   # the existing comm223ss aborts for a parity-changing Y (separate fix); nothing to compare to
        if cut and term in CUT_BUGGY_EXISTING:
            continue   # nothing correct to compare to until the cut fixes are merged; remove this skip then
        plan(False)
        Z0 = new_Z(ms, ty, py)
        getattr(Commutator, term)(X, Y, Z0)
        Z1 = new_Z(ms, ty, py)
        getattr(Commutator, term + suffix)(X, Y, Z1)
        check(f'{tag} T={ty} P={py} {term}{suffix}', Z1, Z0)


# term by term
ms = make_ms(2)
for ty, py in ((0, 0), (1, 0), (0, 1), (1, 1)):
    kernels('emax 2', ms, ty, py)
for ty, py in ((0, 0), (1, 0)):
    kernels('emax 2 E3max 4', make_ms(2, e3max=4), ty, py)
    kernels('emax 2 dE3max 2', make_ms(2, de3max=2.0), ty, py, cut=True)
for ty, py in ((0, 0), (1, 0), (0, 1)):
    kernels('emax 3 emax_3body 2', make_ms(3, emax3b=2), ty, py)

# IMSRG(2) terms with block kernels for a channel-changing Y (comm122ss_blocks, comm222_phss_blocks), plan on vs off
ms3 = make_ms(3)
for ty, py in ((1, 0), (0, 1), (1, 1)):
    ut, X, Y = random_ops(ms3, ty, py)
    for term in ('comm122ss', 'comm222_phss'):
        Zs = []
        for on in (False, True):
            plan(on)
            Zs.append(new_Z(ms3, ty, py))
            getattr(Commutator, term)(X, Y, Zs[-1])
        check(f'emax 3 T={ty} P={py} {term} blocks', Zs[1], Zs[0])
plan(False)

# full commutator, plan on vs off (comm223ss aborts in the existing code for parity-changing Y, so P = 0 only)
for tag, msf in (('emax 2', ms), ('emax 3 emax_3body 2', make_ms(3, emax3b=2))):
    for ty in (0, 1):
        ut, X, Y = random_ops(msf, ty, 0)
        plan(False)
        Z0 = Commutator.Commutator(X, Y)
        plan(True)
        Z1 = Commutator.Commutator(X, Y)
        check(f'{tag} T={ty} full Commutator', Z1, Z0)

# ReferenceImplementations tests, with the plan on
plan(True)
ms1 = make_ms(1)
ut, X, Y = random_ops(ms1, 0, 0, seed=11)
ok = ut.TestCommutators(X, Y)
ok &= ut.TestCommutators_IsospinChanging()
print(f'{"ok  " if ok else "FAIL"} emax 1 TestCommutators + TestCommutators_IsospinChanging with the plan on', flush=True)
if not ok:
    failures.append('reference tests')
for ty in (0, 1):   # comm223 with a parity-changing Y: comm223ss_plan_tz vs ReferenceImplementations::comm223ss
    ut, X, Y = random_ops(ms1, ty, 1)
    ok = ut.Test_comm223ss(X, Y)
    print(f'{"ok  " if ok else "FAIL"} emax 1 T={ty} P=1 Test_comm223ss with the plan on', flush=True)
    if not ok:
        failures.append(f'Test_comm223ss T={ty} P=1')
plan(False)

passed = not failures
print('passed? ', passed, '' if passed else failures)
exit(not passed)
