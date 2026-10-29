"""Generate the explanatory figures from the actual formulas (PNG and SVG).

python tools/generate_doc_figures.py
Requires requirements-docs.txt and a Chinese font (Noto Sans CJK or Microsoft YaHei).
"""
from pathlib import Path
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, Arc
from scipy.interpolate import BSpline
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
import reference as ref

OUT = ROOT / 'docs' / 'images'
OUT.mkdir(parents=True, exist_ok=True)
FONT_FILES = [Path('C:/Windows/Fonts/msyh.ttc'), Path('/mnt/c/Windows/Fonts/msyh.ttc'),
              Path('/System/Library/Fonts/PingFang.ttc')]
font_path = next((str(p) for p in FONT_FILES if p.exists()), None)
if font_path is None:
    for family in ['Noto Sans CJK SC', 'Noto Sans CJK JP', 'Microsoft YaHei']:
        try:
            font_path = font_manager.findfont(family, fallback_to_default=False)
            break
        except ValueError:
            pass
if font_path is None:
    raise RuntimeError('Install a Chinese font, for example fonts-noto-cjk.')
font_manager.fontManager.addfont(font_path)
plt.rcParams.update({'font.family': font_manager.FontProperties(fname=font_path).get_name(),
                     'font.size': 12, 'axes.unicode_minus': False,
                     'svg.fonttype': 'path', 'savefig.facecolor': 'white'})
BLUE, ORANGE, GREEN, DARK, GRAY = '#1764a0', '#b45309', '#127047', '#17212b', '#718096'


def save(fig, name):
    for extension in ['png', 'svg']:
        path = OUT / f'{name}.{extension}'
        fig.savefig(path, dpi=170, bbox_inches='tight', pad_inches=.25)
        if extension == 'svg':
            # Matplotlib emits trailing spaces inside multiline path data.
            # Keep newline separators while making generated text Git-clean.
            path.write_text('\n'.join(line.rstrip() for line in path.read_text(encoding='utf-8').splitlines())
                            + '\n', encoding='utf-8')
    plt.close(fig)


def arrow(ax, a, b, color=BLUE, lw=2, **kwargs):
    ax.annotate('', xy=b, xytext=a,
                arrowprops=dict(arrowstyle='->', color=color, lw=lw, **kwargs))


def mark(ax, point, text, offset=(8, 8), color=DARK):
    ax.scatter(*point, s=32, color=color, zorder=6)
    ax.annotate(text, point, xytext=offset, textcoords='offset points', color=DARK,
                ha='left' if offset[0] >= 0 else 'right', va='bottom' if offset[1] >= 0 else 'top')


def clean(ax):
    ax.set_aspect('equal')
    ax.axis('off')


geo = ref.build_trajectory_geometry((500., -150.), (0., 0.), 150., -1)
A, O, B, C = (geo[k] for k in ['A', 'O', 'B', 'C'])
D = float(np.linalg.norm(A-O))
theta = float(np.arccos(150 / D))
ct, st = float(np.cos(theta)), float(np.sin(theta))

# 1. Actual example geometry; no circular direction is implied by this diagram.
fig, ax = plt.subplots(figsize=(10.4, 4.9))
ax.add_patch(Circle(O, 150, fill=False, color=BLUE, lw=2))
ax.plot([C[0], A[0]], [C[1], A[1]], color=DARK, lw=2)
ax.plot([O[0], A[0]], [O[1], A[1]], color=GRAY, ls='--', lw=1.5)
ax.plot([O[0], B[0]], [O[1], B[1]], color=ORANGE, lw=2)
ax.plot([0, 17, 17], [-133, -133, -150], color=GRAY)
ax.add_patch(Arc(O, 105, 105, theta1=-90, theta2=np.degrees(np.arctan2(-150, 500)),
                 color=ORANGE, lw=1.4))
ax.text(47, -68, r'$\theta$', color=ORANGE, fontsize=16)
mark(ax, A, 'A：起点\n(500, −150)', (-8, -15))
mark(ax, O, 'O：圆心\n(0, 0)', (8, 9))
mark(ax, B, 'B：切点 (0, −150)\n绕一整圈后仍回到这里', (13, -16))
mark(ax, C, 'C_end：终点\n(−500, −150)', (0, -15))
ax.text(325, -30, f'D = |OA| ≈ {D:.3f}', rotation=-16.7, ha='center', color=GRAY)
ax.text(-16, -75, 'OB 的长度\nR = 150', ha='right', color=ORANGE)
arrow(ax, (380, -150), (235, -150), DARK)
arrow(ax, (-180, -150), (-340, -150), DARK)
ax.text(270, -125, 'AB = 500', ha='center')
ax.text(-270, -125, 'BC_end = 500', ha='center')
ax.text(0, 195, 'A → B → 绕圆一整圈回到 B → C_end', ha='center', fontsize=15)
ax.set(xlim=(-555, 555), ylim=(-265, 225))
clean(ax)
save(fig, '01-geometry')

# 2. Basis coordinates, two tangent branches, and cos/sin projections.
fig, ax = plt.subplots(figsize=(7, 6.4))
ax.add_patch(Circle((0, 0), 1, fill=False, lw=1.5, color=GRAY))
arrow(ax, (-.2, 0), (1.4, 0), GRAY, 1.3)
arrow(ax, (0, -1.2), (0, 1.3), GRAY, 1.3)
ax.text(1.05, .08, 'e：朝向 A 的方向', ha='center')
ax.text(-.04, 1.31, 'n：e 逆时针转 90°', ha='center')
for sign, label, color in [(-1, 'branch = −1', BLUE), (1, 'branch = +1', GREEN)]:
    pos = (ct, sign*st)
    arrow(ax, (0, 0), pos, color, 2)
    ax.plot([ct, ct], [0, sign*st], ls='--', color=color, lw=1.3)
    ax.plot([0, ct], [sign*st, sign*st], ls=':', color=GRAY)
    ax.scatter(*pos, color=color, s=34)
    ax.text(ct+.1, sign*st, label+'\n切点的单位方向', va='center', color=color)
ax.add_patch(Arc((0, 0), .75, .75, theta1=-np.degrees(theta), theta2=0, color=BLUE))
ax.text(.38, -.25, 'θ ≈ 73.30°', color=BLUE)
ax.plot([0, ct], [0, 0], color=ORANGE, lw=4)
ax.annotate('ct = cos θ\n≈ 0.2873', (.5*ct, 0), xytext=(-.55, .25),
            ha='center', color=ORANGE, arrowprops=dict(arrowstyle='->', color=ORANGE))
ax.text(-.14, -.53, 'st = sin θ\n≈ 0.9578\n下方取负号', ha='right', color=BLUE)
ax.text(.4, -1.32, '这是局部 e、n 坐标；不是世界 x、y 坐标。', ha='center')
ax.set(xlim=(-1.15, 1.75), ylim=(-1.5, 1.5))
clean(ax)
save(fig, '02-tangent-basis')

# 3. Deliberately coarse controls make the smoothing visible.
P, h, _ = ref.generate_trajectory_points(A, O, 150., 80., -1)
U = ref.generate_knot_vector(P, 4)
curve = BSpline(U, P, 4)
values = curve(np.linspace(0, 1, 3000))
fig, ax = plt.subplots(figsize=(9.7, 5.5))
angle = np.linspace(0, 2*np.pi, 800)
ax.plot(150*np.cos(angle), 150*np.sin(angle), color=GRAY, ls='--', label='理论路径：圆与切线')
ax.plot([-500, 500], [-150, -150], color=GRAY, ls='--')
ax.plot(P[:, 0], P[:, 1], color=ORANGE, lw=.8, alpha=.6, label='相邻控制点连线')
ax.scatter(P[:, 0], P[:, 1], s=27, facecolors='white', edgecolors=ORANGE, zorder=5,
           label='控制点：来自理论路径')
ax.plot(values[:, 0], values[:, 1], color=BLUE, lw=2.4, label='实际生成的四次 B 样条')
index = int(np.argmax(P[:, 1]))
ax.annotate('曲线通常不穿过内部控制点', P[index], xytext=(300, 78),
            arrowprops=dict(arrowstyle='->', color=DARK), ha='center')
ax.text(0, -257, f'教学用粗采样：d = 80，实际间距 ≈ {h:.2f}，控制点 {len(P)} 个\n'
        '正式默认结果是 347 个控制点；这里故意放大偏差以看清差别。', ha='center')
ax.set(xlim=(-550, 550), ylim=(-300, 220))
ax.legend(loc='upper center', ncol=2, bbox_to_anchor=(.5, 1.06), frameon=False)
clean(ax)
save(fig, '03-controls-and-curve')

# 4. Six equally spaced-by-chord controls give exactly one internal knot.
small = np.array([[0, 0], [1, 0], [1.6, .8], [2.6, .8], [3.2, 0], [4.2, 0]])
small_U = ref.generate_knot_vector(small, 4)
small_curve = BSpline(small_U, small, 4)
small_t = np.r_[0, np.cumsum(np.linalg.norm(np.diff(small, axis=0), axis=1)) / 5]
fig, axes = plt.subplots(3, 1, figsize=(8.8, 6.9), gridspec_kw={'height_ratios':[2, 1, 1]})
ax = axes[0]
ax.plot(small[:, 0], small[:, 1], color=ORANGE, ls=':', lw=1)
ax.scatter(small[:, 0], small[:, 1], color=ORANGE, s=32)
vv = small_curve(np.linspace(0, 1, 500))
ax.plot(vv[:, 0], vv[:, 1], color=BLUE, lw=2)
for i, p in enumerate(small):
    ax.annotate(f'P{i}', p, xytext=(0, 9), textcoords='offset points', ha='center')
ax.text(2.1, 1.32, 'P：二维空间里的控制点（六点教学例子）', ha='center', fontsize=14)
ax.set(xlim=(-.2, 4.5), ylim=(-.2, 1.6))
clean(ax)
for ax in axes[1:]:
    ax.plot([0, 1], [0, 0], color=GRAY)
    ax.set(xlim=(-.1, 1.12), ylim=(-.7, .8))
    ax.axis('off')
ax = axes[1]
for i, t in enumerate(small_t):
    ax.scatter(t, 0, color=ORANGE, s=28)
    ax.text(t, -.18, f't{i}\n{t:.1f}', ha='center', va='top')
ax.text(.5, .48, 't_i：由相邻控制点弦长得到的参数标签', ha='center')
ax = axes[2]
for x in [0, 1]:
    for j in range(5):
        ax.scatter(x, j*.09, s=18, color=BLUE)
    ax.text(x, -.18, f'{x}\n重复 5 次', ha='center', va='top')
ax.scatter(.5, 0, color=BLUE, s=28)
ax.text(.5, -.18, '0.5\n唯一内部节点', ha='center', va='top')
ax.text(.5, .58, 'U：一维参数轴上的分界值（不是空间坐标）', ha='center')
fig.subplots_adjust(hspace=.12)
save(fig, '04-points-parameters-knots')

# 5. Same curve; equally spaced u produce different traveled lengths.
u_marks = np.linspace(0, 1, 5)
locations = small_curve(u_marks)
derivative = small_curve.derivative()
def length_to(u):
    bps = [0] + ([.5] if u > .5 else []) + ([float(u)] if u > 0 else [])
    return sum(quad(lambda z: np.linalg.norm(derivative(z)), a, b,
                    epsabs=1e-12, epsrel=1e-12)[0] for a, b in zip(bps[:-1], bps[1:]))
lengths = [length_to(u) for u in u_marks]
fig, axes = plt.subplots(2, 1, figsize=(8.7, 6.6), gridspec_kw={'height_ratios':[1.2, 1]})
ax = axes[0]
ax.plot(vv[:, 0], vv[:, 1], color=BLUE, lw=2)
for u, p in zip(u_marks, locations):
    mark(ax, p, f'u = {u:g}', (-8, 12) if u < .5 else (8, 12))
for i in range(4):
    mid = small_curve((u_marks[i]+u_marks[i+1])/2)
    ax.text(mid[0], -.22, f'走 {lengths[i+1]-lengths[i]:.3f}', ha='center')
ax.text(2.1, 1.15, '相同 Δu = 0.25，对应不同的曲线距离', ha='center', fontsize=14)
ax.set(xlim=(-.4, 4.7), ylim=(-.6, 1.45))
clean(ax)
ax = axes[1]
dense_u = np.linspace(0, 1, 101)
ax.plot(dense_u, [length_to(u) for u in dense_u], color=BLUE, lw=2)
ax.scatter(u_marks, lengths, color=BLUE, s=25)
ax.set(xlabel='参数 u（无时间单位）', ylabel='累计弧长 S(u)（长度单位）', xlim=(0, 1))
ax.grid(alpha=.15)
fig.tight_layout(h_pad=1.5)
save(fig, '05-parameter-and-arc-length')

# 6. Illustrative error point, not a claim about the actual fitted error.
fig, ax = plt.subplots(figsize=(7.4, 6.1))
ax.add_patch(Circle(O, 150, fill=False, color=BLUE, lw=2))
ax.plot([-220, 230], [-150, -150], color=GRAY, lw=1.5)
q = np.array([60., -100.])
q_line = np.array([60., -150.])
q_circle = 150*q/np.linalg.norm(q)
ax.plot([0, q_circle[0]], [0, q_circle[1]], color=GRAY, ls=':', lw=1)
ax.plot([q[0], q_line[0]], [q[1], q_line[1]], color=GREEN, lw=2.5)
ax.plot([q[0], q_circle[0]], [q[1], q_circle[1]], color=ORANGE, lw=2.5)
mark(ax, O, '圆心 O', (-9, 5))
mark(ax, q, '示意点 Q', (12, 12))
mark(ax, q_line, '直线投影', (-20, -14), GREEN)
mark(ax, q_circle, '圆上投影', (18, 10), ORANGE)
ax.annotate('到直线：50', (60, -124), xytext=(-120, -80),
            arrowprops=dict(arrowstyle='->', color=GREEN), color=GREEN)
circle_error = float(150-np.linalg.norm(q))
ax.annotate(f'到圆：{circle_error:.2f}', (q+q_circle)/2, xytext=(130, -45),
            arrowprops=dict(arrowstyle='->', color=ORANGE), color=ORANGE)
ax.text(0, 190, '两种距离取较小值', ha='center', fontsize=15)
ax.text(0, -255, f'这个示意点的误差 = min(50, {circle_error:.2f}) = {circle_error:.2f}\n'
        '这里只解释定义；不是正式拟合的 0.04～0.05 误差结果。', ha='center')
ax.set(xlim=(-235, 255), ylim=(-290, 220))
clean(ax)
save(fig, '06-error-definition')

weights = BSpline(small_U, np.eye(6), 4)(.25)
data = dict(distance_D=D, theta_degrees=float(np.degrees(theta)), ct=ct, st=st,
            e=((A-O)/D).tolist(), normal=[float(-(A-O)[1]/D), float((A-O)[0]/D)],
            first_B_length=geo['s_ab'], second_B_length=geo['s_ab']+geo['s_circle'],
            theoretical_total_length=geo['s_total'], teaching_controls=small.tolist(),
            teaching_parameters=small_t.tolist(), teaching_knots=small_U.tolist(),
            weights_at_quarter=weights.tolist(), point_at_quarter=small_curve(.25).tolist(),
            u_marks=u_marks.tolist(), cumulative_lengths=lengths)
(OUT/'figure_data.json').write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
print(json.dumps(data, ensure_ascii=False, indent=2))
print('Generated six figures as PNG and SVG in docs/images.')
