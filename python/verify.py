"""Compare C with supplied Python formulas and independent SciPy/QUADPACK.

Usage: python python/verify.py --exe build/bspline_demo --report docs/precision_report.json
"""
import argparse
import json
import platform
from pathlib import Path
import subprocess
import sys
import warnings
import numpy as np
import scipy
from scipy.integrate import quad, IntegrationWarning
from scipy.interpolate import BSpline
from scipy.optimize import brentq
import reference as ref


def independent_arc(P, U):
    curve = BSpline(U, P, 4)
    derivative = curve.derivative()
    breaks = np.unique(U)
    def speed(u):
        return np.linalg.norm(derivative(u))
    def integrate(a, b):
        if a == b:
            return 0.
        with warnings.catch_warnings():
            warnings.simplefilter('error', IntegrationWarning)
            value, _ = quad(speed, a, b, epsabs=2e-11, epsrel=2e-13, limit=200)
        return value
    prefix = np.r_[0., np.cumsum([integrate(a, b) for a, b in zip(breaks[:-1], breaks[1:])])]
    def cumulative(u):
        i = np.clip(np.searchsorted(breaks, u, side='right') - 1, 0, len(breaks) - 2)
        return prefix[i] + integrate(breaks[i], u)
    return curve, derivative, cumulative, float(prefix[-1])


def check(exe, case):
    args = case['args']
    continuous = case.get('continuous', False)
    command = [exe, '--json'] + (['--continuous'] if continuous else []) + [str(x) for x in args]
    process = subprocess.run(command, text=True, capture_output=True, check=True)
    actual = json.loads(process.stdout)
    ax, ay, ox, oy, radius, branch, d_init, emin, emax = args
    scale = max(1., radius / 150.)
    d, e, P, U, geo, _ = ref.fit_bspline_by_bisection(
        (ax, ay), (ox, oy), radius, emin, emax, branch=branch, d_init=d_init,
        continuous=continuous)
    metrics = {}
    def test(name, error, tolerance):
        error = float(error)
        if not np.isfinite(error) or error > tolerance:
            raise AssertionError(f'{case["name"]}: {name}={error:.6e} > {tolerance:.6e}')
        metrics[name] = dict(error=error, tolerance=tolerance)
    test('d', abs(actual['d'] - d), 1e-12 * scale)
    test('control_points', np.max(np.abs(np.asarray(actual['controls']) - P)), 1e-9 * scale)
    test('knots', np.max(np.abs(np.asarray(actual['knots']) - U)), 5e-13)
    test('sampled_max_error', abs(actual['max_error'] - e), 1e-9 * scale)
    assert emin <= actual['max_error'] <= emax
    samples = np.asarray(actual['samples'])
    us = samples[:, 0]
    python_points = np.array([ref.de_boor(u, P, 4, U) for u in us])
    Q, Up = ref.compute_derivative_curve(P, 4, U)
    python_derivatives = np.array([ref.de_boor(u, Q, 3, Up) for u in us])
    test('positions_vs_python', np.max(np.abs(samples[:, 1:3] - python_points)), 2e-9 * scale)
    test('derivatives_vs_python', np.max(np.abs(samples[:, 3:5] - python_derivatives)), 2e-6 * scale)
    independent, derivative, cumulative, total = independent_arc(P, U)
    test('positions_vs_scipy', np.max(np.abs(samples[:, 1:3] - independent(us))), 2e-9 * scale)
    test('derivatives_vs_scipy', np.max(np.abs(samples[:, 3:5] - derivative(us))), 2e-6 * scale)
    independent_cumulative = np.array([cumulative(u) for u in us])
    test('cumulative_vs_quadpack', np.max(np.abs(samples[:, 5] - independent_cumulative)), 1e-8 * scale)
    test('total_vs_quadpack', abs(actual['length'] - total), 1e-8 * scale)
    # Original GL has an unconditional n_max exit; record its accuracy separately.
    old_length, old_info = ref.arc_length_bspline_gauss(P, 4, U, err=1e-6)
    original_length_error = abs(old_length - total)
    u_old, point_old, length_old, _ = ref.find_u_by_arc_length(.5, 1.2, P, 4, U)
    inv = actual['inverse']
    test('default_inverse_point_vs_python', np.max(np.abs(np.asarray(inv['point']) - point_old)), 3e-6 * scale)
    test('default_inverse_length_vs_python', abs(inv['length'] - length_old), 3e-6 * scale)
    test('default_inverse_true_residual', abs(cumulative(inv['u']) - cumulative(.5) - 1.2), 1.01e-6)
    max_residual = 0.
    max_u_error = 0.
    for u0, target, u, x, y, length in actual['inverses']:
        base = cumulative(u0)
        residual = abs(cumulative(u) - base - target)
        expected_u = brentq(lambda v: cumulative(v) - base - target, u0, 1., xtol=1e-14)
        max_residual = max(max_residual, residual)
        max_u_error = max(max_u_error, abs(u - expected_u))
    test('inverse_true_residual', max_residual, 1.1e-8 * scale)
    test('inverse_parameter_vs_brent', max_u_error, 1e-9)
    dense = np.linspace(0, 1, 200001)
    # Matches the user's dense speed trapezoid, vectorized for runtime only.
    dense_speed = np.linalg.norm(ref.evaluate_many(dense, Q, 3, Up), axis=1)
    trapezoid = getattr(np, 'trapezoid', None) or np.trapz
    dense_length = float(trapezoid(dense_speed, dense))
    dense_error = abs(dense_length - total)
    test('dense_trapezoid_reference', dense_error, 5e-4 * scale)
    dense_points = ref.evaluate_many(dense, P, 4, U)
    relative = dense_points - geo['A']
    v = geo['AC']
    line_projection = np.outer(relative @ v / np.dot(v, v), v)
    line_errors = np.linalg.norm(relative - line_projection, axis=1)
    circle_errors = np.abs(np.linalg.norm(dense_points - geo['O'], axis=1) - radius)
    dense_geometric_error = float(np.max(np.minimum(line_errors, circle_errors)))
    test('dense_geometric_sample_max', dense_geometric_error, emax)
    print(f'{case["name"]}: PASS, controls={len(P)}, C length error={abs(actual["length"]-total):.3e}, '
          f'original GL error={original_length_error:.3e}, inverse residual={max_residual:.3e}', flush=True)
    return dict(name=case['name'], args=args, continuous=continuous, passed=True,
                d=d, controls=len(P), max_error=e, c_length=actual['length'],
                independent_length=total, original_python_length=old_length,
                original_python_length_error=original_length_error,
                original_gl_difference_estimate=sum(old_info['local_err']),
                original_gl_max_order=max(old_info['n_used']),
                dense_trapezoid_error=dense_error, dense_geometric_max_error=dense_geometric_error,
                dense_geometric_within_interval=bool(emin<=dense_geometric_error<=emax), metrics=metrics)


def cases():
    yield dict(name='original-input', args=[500., -150., 0., 0., 150., -1, 20., .04, .05])
    yield dict(name='opposite-branch', args=[500., -150., 0., 0., 150., 1, 20., .04, .05])
    yield dict(name='translated', args=[623., -227., 123., -77., 150., -1, 20., .04, .05])
    angle = .73
    ax = 500*np.cos(angle) + 150*np.sin(angle)
    ay = 500*np.sin(angle) - 150*np.cos(angle)
    yield dict(name='rotated', args=[float(ax), float(ay), 0., 0., 150., -1, 20., .04, .05])
    for scale in [.1, 10.]:
        yield dict(name=f'scale-{scale}', args=[500*scale, -150*scale, 0., 0., 150*scale,
                                                -1, 20*scale, .04*scale, .05*scale])
    for branch in [-1, 1]:
        yield dict(name=f'continuous-{branch}', continuous=True,
                   args=[500., -150., 0., 0., 150., branch, 20., .04, .05])
    rng = np.random.default_rng(20261029)
    for i in range(4):
        radius = float(rng.uniform(50, 300))
        angle = float(rng.uniform(-np.pi, np.pi))
        distance = float(radius * rng.uniform(1.5, 5))
        center = rng.uniform(-200, 200, size=2)
        a = center + distance * np.array([np.cos(angle), np.sin(angle)])
        scale = radius / 150
        yield dict(name=f'seeded-{i}', args=[float(a[0]), float(a[1]), float(center[0]),
                   float(center[1]), radius, -1 if i % 2 else 1, 20*scale, .04*scale, .05*scale])


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--exe', default='build/bspline_demo')
    parser.add_argument('--report', default='docs/precision_report.json')
    options = parser.parse_args()
    results = [check(str(Path(options.exe).resolve()), case) for case in cases()]
    report = dict(passed=True, case_count=len(results), environment=dict(
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        platform=platform.system()+' '+platform.machine()),
        notes=['All tolerances are absolute and expressed in the input coordinate units.',
               'GL successive-order differences are estimates, not rigorous error bounds.',
               'A sampled maximum is not a certified continuous-curve maximum.'], cases=results)
    Path(options.report).parent.mkdir(parents=True, exist_ok=True)
    Path(options.report).write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(f'All {len(results)} cases passed; report: {options.report}')


if __name__ == '__main__':
    main()
