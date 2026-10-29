"""Render a C result without putting plotting dependencies into the C library."""
import argparse
import json
import subprocess
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(__doc__)
parser.add_argument('--exe', default='build/bspline_demo')
parser.add_argument('--output', default='trajectory.png')
parser.add_argument('--continuous', action='store_true')
args = parser.parse_args()
command = [args.exe, '--json'] + (['--continuous'] if args.continuous else [])
data = json.loads(subprocess.check_output(command, text=True))
controls, samples = np.asarray(data['controls']), np.asarray(data['samples'])
theta = np.linspace(0, 2*np.pi, 500)
fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
ax.plot(150*np.cos(theta), 150*np.sin(theta), 'k--', label='Theoretical circle')
ax.plot([500, -500], [-150, -150], 'k-', label='Theoretical line')
ax.scatter(controls[:, 0], controls[:, 1], s=5, color='red', label='Control points')
ax.plot(samples[:, 1], samples[:, 2], color='blue', label='Quartic B-spline')
ax.set(aspect='equal', xlabel='x', ylabel='y', title=f'd={data["d"]:.6g}, sampled error={data["max_error"]:.6g}')
ax.grid(alpha=.3)
ax.legend()
fig.savefig(args.output, dpi=160)
print(args.output)
