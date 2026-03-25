"""
AdaptSim -- Publication-quality figures for Lancet Digital Health manuscript.

Figure 1: Cumulative alpha-spending curves (4 spending functions)
Figure 2: Type I error convergence by simulation count
Figure 3: Power vs effect size for different numbers of interim analyses

All computations use the same formulas as adaptsim.html.
"""

import sys
import io
import os
import math
import numpy as np
from scipy import stats

# UTF-8 stdout for Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

# -- House style ---------------------------------------------------------------
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'lines.linewidth': 1.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'mathtext.fontset': 'dejavuserif',
})

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# -- Colour palette (Lancet-style, colour-blind safe) --------------------------
COL_OBF     = '#0072B2'   # blue
COL_POCOCK  = '#D55E00'   # vermilion
COL_LD      = '#009E73'   # teal
COL_HSD     = '#CC79A7'   # rose
COL_NOMINAL = '#999999'   # grey for reference lines
COL_K = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00']  # for K=1..5


# ==============================================================================
# Core math (matches adaptsim.html exactly)
# ==============================================================================

def normal_cdf(x):
    return stats.norm.cdf(x)

def normal_quantile(p):
    return stats.norm.ppf(p)

def normal_pdf(x):
    return stats.norm.pdf(x)


def alpha_spend_obf(t, alpha):
    """O'Brien-Fleming spending function."""
    if t <= 0:
        return 0.0
    if t >= 1:
        return alpha
    za2 = normal_quantile(1 - alpha / 2)
    return 2 - 2 * normal_cdf(za2 / math.sqrt(t))


def alpha_spend_pocock(t, alpha):
    """Pocock spending function."""
    if t <= 0:
        return 0.0
    if t >= 1:
        return alpha
    return alpha * math.log(1 + (math.e - 1) * t)


def alpha_spend_landemets(t, alpha, rho=3):
    """Lan-DeMets spending (power family): alpha * t^rho."""
    if t <= 0:
        return 0.0
    if t >= 1:
        return alpha
    return min(alpha, alpha * t ** rho)


def alpha_spend_hsd(t, alpha, gamma=-4):
    """Hwang-Shih-DeCani spending function."""
    if t <= 0:
        return 0.0
    if t >= 1:
        return alpha
    if abs(gamma) < 1e-8:
        return alpha * t
    return alpha * (1 - math.exp(-gamma * t)) / (1 - math.exp(-gamma))


# ==============================================================================
# Grid-based boundary computation (efficient for K up to 10+)
# ==============================================================================
#
# Instead of recursively nested quadrature (O(M^K) for K looks with M quad points),
# we discretize the density on a fine grid and carry it forward at each look.
# This is O(K * M * N_grid) where N_grid ~ 1000, much faster.

N_GRID = 1201  # grid points for density propagation
Z_RANGE = (-6, 6)  # range for z-grid


def _make_grid():
    """Create z-grid and spacing."""
    z = np.linspace(Z_RANGE[0], Z_RANGE[1], N_GRID)
    dz = z[1] - z[0]
    return z, dz


def compute_boundaries_fast(alpha, info_fracs, spend_fn, two_sided=False):
    """
    Compute group-sequential efficacy boundaries using grid-based density
    propagation. Much faster than recursive quadrature for K >= 3.
    """
    K = len(info_fracs)
    z_grid, dz = _make_grid()
    boundaries = []
    alpha_spent_prev = 0.0

    # Initial density: standard normal PDF on the grid
    density = stats.norm.pdf(z_grid)

    for k in range(K):
        t = info_fracs[k]
        target_alpha = spend_fn(t, alpha)
        delta_alpha = target_alpha - alpha_spent_prev
        if delta_alpha < 1e-15:
            delta_alpha = 1e-15

        if k == 0:
            # First look: find z_eff from tail probability
            if two_sided:
                z_eff = normal_quantile(1 - delta_alpha / 2)
            else:
                z_eff = normal_quantile(1 - delta_alpha)

            # Record boundary
            boundaries.append({
                'look': k + 1,
                'info_frac': t,
                'z_eff': z_eff,
                'z_fut': -np.inf,
                'cumAlpha': target_alpha,
            })

            # Truncate density: zero out rejected region
            if two_sided:
                mask = (z_grid >= -z_eff) & (z_grid <= z_eff)
            else:
                mask = z_grid <= z_eff
            density = density * mask

        else:
            t_prev = info_fracs[k - 1]
            rho = math.sqrt(t_prev / t)
            sigma = math.sqrt(1 - t_prev / t)

            # Propagate density forward: Z_k | Z_{k-1}=z ~ N(rho*z, sigma^2)
            # new_density(z') = integral density(z) * phi((z' - rho*z)/sigma)/sigma dz
            # This is a convolution-like operation; we compute it via matrix multiply
            # on the grid.
            new_density = np.zeros(N_GRID)
            for i in range(N_GRID):
                zp = z_grid[i]
                cond_means = rho * z_grid  # conditional mean for each z in grid
                cond_pdf = stats.norm.pdf((zp - cond_means) / sigma) / sigma
                new_density[i] = np.sum(density * cond_pdf) * dz

            density = new_density

            # Find z_eff via bisection on the tail probability of the propagated density
            # P(Z_k > c_k, survived) = integral of density(z) for z > c_k
            # (for two-sided: also z < -c_k)

            # Use the grid to compute tail probability as a function of c
            def tail_prob(c):
                if two_sided:
                    mask = (z_grid > c) | (z_grid < -c)
                else:
                    mask = z_grid > c
                return np.sum(density[mask]) * dz

            # Bisection
            z_lo, z_hi = 0.0, 5.5
            for _ in range(60):
                z_mid = (z_lo + z_hi) / 2
                prob = tail_prob(z_mid)
                if prob > delta_alpha:
                    z_lo = z_mid
                else:
                    z_hi = z_mid
            z_eff = (z_lo + z_hi) / 2

            boundaries.append({
                'look': k + 1,
                'info_frac': t,
                'z_eff': z_eff,
                'z_fut': -np.inf,
                'cumAlpha': target_alpha,
            })

            # Truncate density: zero out rejected region
            if two_sided:
                mask = (z_grid >= -z_eff) & (z_grid <= z_eff)
            else:
                mask = z_grid <= z_eff
            density = density * mask

        alpha_spent_prev = target_alpha

    return boundaries


def compute_power_analytical(boundaries, info_fracs, theta, two_sided=False):
    """
    Compute power analytically using grid-based density propagation under
    drift theta.  Z_k | Z_{k-1}=z ~ N(rho*z + drift_inc, sigma^2)
    where drift_inc = theta*(sqrt(t_k) - rho*sqrt(t_{k-1})).
    """
    K = len(info_fracs)
    z_grid, dz = _make_grid()

    # Initial density under H1: Z_1 ~ N(theta*sqrt(t_1), 1)
    density = stats.norm.pdf(z_grid, loc=theta * math.sqrt(info_fracs[0]), scale=1.0)

    total_rejection = 0.0

    for k in range(K):
        t = info_fracs[k]
        b = boundaries[k]
        z_eff = b['z_eff']

        if k > 0:
            t_prev = info_fracs[k - 1]
            rho = math.sqrt(t_prev / t)
            sigma = math.sqrt(1 - t_prev / t)
            drift_inc = theta * (math.sqrt(t) - rho * math.sqrt(t_prev))

            # Propagate density
            new_density = np.zeros(N_GRID)
            for i in range(N_GRID):
                zp = z_grid[i]
                cond_means = rho * z_grid + drift_inc
                cond_pdf = stats.norm.pdf((zp - cond_means) / sigma) / sigma
                new_density[i] = np.sum(density * cond_pdf) * dz
            density = new_density

        # Rejection probability at this look
        if two_sided:
            rej_mask = (z_grid > z_eff) | (z_grid < -z_eff)
        else:
            rej_mask = z_grid > z_eff
        rej_prob = np.sum(density[rej_mask]) * dz
        total_rejection += rej_prob

        # Truncate for continuation
        if two_sided:
            cont_mask = (z_grid >= -z_eff) & (z_grid <= z_eff)
        else:
            cont_mask = z_grid <= z_eff
        density = density * cont_mask

    return total_rejection


# ==============================================================================
# Vectorized Monte Carlo simulation (for Figure 2 convergence)
# ==============================================================================

def monte_carlo_power(boundaries, info_fracs, theta, n_sims, seed=42,
                      two_sided=False):
    """
    Vectorized Monte Carlo: simulate n_sims trials simultaneously using numpy.
    Returns rejection rate under drift theta.
    """
    rng = np.random.RandomState(seed)
    K = len(info_fracs)

    active = np.ones(n_sims, dtype=bool)
    rejected = np.zeros(n_sims, dtype=bool)
    z_prev = np.zeros(n_sims)
    t_prev = 0.0

    for k in range(K):
        t = info_fracs[k]
        b = boundaries[k]
        n_active = int(active.sum())
        if n_active == 0:
            break

        if k == 0:
            z = rng.normal(loc=theta * math.sqrt(t), scale=1.0, size=n_sims)
        else:
            rho_k = math.sqrt(t_prev / t)
            sigma_k = math.sqrt(1 - t_prev / t)
            drift_inc = theta * (math.sqrt(t) - rho_k * math.sqrt(t_prev))
            noise = rng.normal(0, sigma_k, size=n_sims)
            z = rho_k * z_prev + drift_inc + noise

        # Efficacy
        if two_sided:
            eff_mask = active & (np.abs(z) >= b['z_eff'])
        else:
            eff_mask = active & (z >= b['z_eff'])
        rejected[eff_mask] = True
        active[eff_mask] = False

        # Futility
        if b['z_fut'] > -10 and k < K - 1:
            if two_sided:
                fut_mask = active & (np.abs(z) <= b['z_fut'])
            else:
                fut_mask = active & (z <= b['z_fut'])
            active[fut_mask] = False

        z_prev = z
        t_prev = t

    return rejected.sum() / n_sims


# ==============================================================================
# FIGURE 1: Alpha-Spending Curves
# ==============================================================================

def figure1():
    print("Generating Figure 1: Alpha-spending curves ...")
    alpha = 0.025  # one-sided

    t_vals = np.linspace(0, 1, 500)

    obf_vals     = [alpha_spend_obf(t, alpha)                for t in t_vals]
    pocock_vals  = [alpha_spend_pocock(t, alpha)              for t in t_vals]
    ld_vals      = [alpha_spend_landemets(t, alpha, rho=3)    for t in t_vals]
    hsd_vals     = [alpha_spend_hsd(t, alpha, gamma=-4)       for t in t_vals]

    fig, ax = plt.subplots(figsize=(5.5, 4.0))

    ax.plot(t_vals, obf_vals,    color=COL_OBF,    label="O'Brien-Fleming",        linestyle='-')
    ax.plot(t_vals, pocock_vals, color=COL_POCOCK,  label='Pocock',                 linestyle='--')
    ax.plot(t_vals, ld_vals,     color=COL_LD,      label=r'Lan-DeMets ($\rho$=3)', linestyle='-.')
    ax.plot(t_vals, hsd_vals,    color=COL_HSD,     label=r'HSD ($\gamma$=$-$4)',   linestyle=':',
            linewidth=2.2)

    # Reference: identity line (uniform spending)
    ax.plot([0, 1], [0, alpha], color=COL_NOMINAL, linewidth=0.8, linestyle='--',
            alpha=0.5, label='Uniform (reference)')

    # Horizontal line at nominal alpha
    ax.axhline(y=alpha, color=COL_NOMINAL, linewidth=0.6, linestyle=':', alpha=0.4)
    ax.text(0.02, alpha + 0.0008, r'$\alpha$ = 0.025', fontsize=8, color=COL_NOMINAL)

    ax.set_xlabel('Information fraction ($t$)')
    ax.set_ylabel(r'Cumulative $\alpha$ spent, $\alpha^*(t)$')
    ax.set_title('Figure 1. Alpha-spending functions for group-sequential designs',
                 fontsize=10, fontweight='bold', loc='left')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, alpha * 1.15)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_locator(MultipleLocator(0.005))

    ax.legend(loc='upper left', frameon=True, framealpha=0.9, edgecolor='#cccccc',
              borderpad=0.6)

    # Annotations: mark spending at t=0.5
    for fn, col, label_txt, yoff in [
        (alpha_spend_obf,     COL_OBF,    'OBF',    0.0004),
        (alpha_spend_pocock,  COL_POCOCK, 'Pocock', -0.0010),
    ]:
        val = fn(0.5, alpha)
        ax.plot(0.5, val, 'o', color=col, markersize=4, zorder=5)
        ax.annotate(f'{val:.4f}', xy=(0.5, val), xytext=(0.55, val + yoff),
                    fontsize=7, color=col,
                    arrowprops=dict(arrowstyle='-', color=col, lw=0.5))

    fig.tight_layout()

    for ext in ['png', 'pdf']:
        path = os.path.join(OUTDIR, f'figure1_alpha_spending.{ext}')
        fig.savefig(path, dpi=300)
        print(f"  Saved: {path}")
    plt.close(fig)


# ==============================================================================
# FIGURE 2: Type I Error Convergence
# ==============================================================================

def figure2():
    print("Generating Figure 2: Type I error convergence ...")
    alpha = 0.025
    info_fracs = [1/3, 2/3, 1.0]
    spend_fn = lambda t, a: alpha_spend_obf(t, a)

    boundaries = compute_boundaries_fast(alpha, info_fracs, spend_fn, two_sided=False)
    bnd_str = ', '.join(f"{b['z_eff']:.4f}" for b in boundaries)
    print(f"  Boundaries: [{bnd_str}]")

    sim_counts = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]

    n_reps = 20
    results = {n: [] for n in sim_counts}

    for n_sims in sim_counts:
        print(f"  Simulating n={n_sims:,} x {n_reps} replicates ...")
        for rep in range(n_reps):
            seed = 1000 * rep + n_sims
            rate = monte_carlo_power(boundaries, info_fracs, theta=0.0,
                                     n_sims=n_sims, seed=seed)
            results[n_sims].append(rate)

    # Plot
    fig, ax = plt.subplots(figsize=(5.5, 4.0))

    x_pos = np.arange(len(sim_counts))
    x_labels = []
    for n in sim_counts:
        if n >= 1000:
            x_labels.append(f'{n//1000}K')
        else:
            x_labels.append(str(n))

    bp_data = [results[n] for n in sim_counts]

    bp = ax.boxplot(bp_data, positions=x_pos, widths=0.5,
                    patch_artist=True, showfliers=True,
                    flierprops=dict(marker='o', markersize=3, markerfacecolor=COL_NOMINAL,
                                    markeredgecolor=COL_NOMINAL, alpha=0.5),
                    medianprops=dict(color='white', linewidth=1.5))

    for patch in bp['boxes']:
        patch.set_facecolor(COL_OBF)
        patch.set_alpha(0.7)
        patch.set_edgecolor(COL_OBF)

    # Nominal alpha line
    ax.axhline(y=alpha, color=COL_POCOCK, linewidth=1.2, linestyle='--',
               label=r'Nominal $\alpha$ = 0.025', zorder=0)

    # Tolerance band
    ax.axhspan(alpha - 0.002, alpha + 0.002, color=COL_LD, alpha=0.12,
               label=r'$\pm$ 0.002 tolerance', zorder=0)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Number of simulated trials')
    ax.set_ylabel('Empirical Type I error rate')
    ax.set_title('Figure 2. Type I error convergence by simulation count',
                 fontsize=10, fontweight='bold', loc='left')

    ax.set_ylim(0.010, 0.040)
    ax.yaxis.set_major_locator(MultipleLocator(0.005))

    ax.legend(loc='upper right', frameon=True, framealpha=0.9, edgecolor='#cccccc')

    # Annotation
    ax.annotate('Converges within\n0.002 of nominal',
                xy=(x_pos[sim_counts.index(10000)], alpha),
                xytext=(x_pos[sim_counts.index(10000)] - 2.5, 0.035),
                fontsize=8, color='#333333',
                arrowprops=dict(arrowstyle='->', color='#333333', lw=0.8),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='#cccccc', alpha=0.9))

    fig.tight_layout()

    for ext in ['png', 'pdf']:
        path = os.path.join(OUTDIR, f'figure2_type1_convergence.{ext}')
        fig.savefig(path, dpi=300)
        print(f"  Saved: {path}")
    plt.close(fig)


# ==============================================================================
# FIGURE 3: Power vs Effect Size for K Interim Analyses
# ==============================================================================

def figure3():
    print("Generating Figure 3: Power vs effect size by K interim analyses ...")
    alpha = 0.025

    # Theta range: drift parameter (= delta * sqrt(I_max))
    theta_vals = np.linspace(0, 3.5, 36)

    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    for K_idx, K in enumerate([1, 2, 3, 4, 5]):
        info_fracs = [(i + 1) / K for i in range(K)]
        spend_fn = lambda t, a: alpha_spend_obf(t, a)

        print(f"  Computing boundaries for K={K} ...")
        boundaries = compute_boundaries_fast(alpha, info_fracs, spend_fn,
                                              two_sided=False)
        bnd_str = ', '.join(f"{b['z_eff']:.4f}" for b in boundaries)
        print(f"    Boundaries: [{bnd_str}]")

        powers = []
        for theta in theta_vals:
            pw = compute_power_analytical(boundaries, info_fracs, theta,
                                           two_sided=False)
            powers.append(min(pw, 1.0))  # clamp to [0, 1]

        idx_2 = int(np.argmin(np.abs(theta_vals - 2.0)))
        print(f"    Power at theta~2.0: {powers[idx_2]:.4f}")

        label_suffix = ' (fixed-sample)' if K == 1 else ''
        lk = 'look' if K == 1 else 'looks'
        ax.plot(theta_vals, powers,
                color=COL_K[K_idx],
                linestyle=['-', '--', '-.', ':', (0, (3, 1, 1, 1))][K_idx],
                linewidth=[2.0, 1.8, 1.8, 2.0, 1.8][K_idx],
                label=f'K = {K} {lk}{label_suffix}')

    # Reference lines
    ax.axhline(y=0.8, color=COL_NOMINAL, linewidth=0.6, linestyle=':', alpha=0.5)
    ax.text(0.05, 0.815, '80% power', fontsize=7, color=COL_NOMINAL)

    ax.axhline(y=0.9, color=COL_NOMINAL, linewidth=0.6, linestyle=':', alpha=0.5)
    ax.text(0.05, 0.915, '90% power', fontsize=7, color=COL_NOMINAL)

    ax.set_xlabel(r'Standardised drift ($\theta$)')
    ax.set_ylabel('Power (probability of rejection)')
    ax.set_title("Figure 3. Power by effect size and number of interim analyses\n"
                 "(O'Brien-Fleming spending, equally spaced looks)",
                 fontsize=10, fontweight='bold', loc='left')

    ax.set_xlim(0, 3.5)
    ax.set_ylim(0, 1.02)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0%}' if y <= 1 else ''))

    ax.legend(loc='lower right', frameon=True, framealpha=0.9, edgecolor='#cccccc',
              borderpad=0.6, title="Interim analyses", title_fontsize=9)

    # Secondary x-axis: approx effect size for N_total=2000
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    delta_ticks = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
    ax2.set_xticks(delta_ticks)
    ax2.set_xticklabels([f'{d / np.sqrt(500):.3f}' for d in delta_ticks])
    ax2.set_xlabel(r'Approximate effect size ($\delta$, $N_{total}$ = 2000)', fontsize=9)
    ax2.tick_params(labelsize=8)
    ax2.spines['top'].set_visible(True)
    ax2.spines['top'].set_linewidth(0.5)
    ax2.spines['top'].set_color('#999999')

    fig.tight_layout()

    for ext in ['png', 'pdf']:
        path = os.path.join(OUTDIR, f'figure3_power_curves.{ext}')
        fig.savefig(path, dpi=300)
        print(f"  Saved: {path}")
    plt.close(fig)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    print(f"Output directory: {OUTDIR}")
    print("=" * 60)

    figure1()
    print()
    figure2()
    print()
    figure3()

    print()
    print("=" * 60)
    print("All 3 figures generated successfully.")
    print(f"Files in {OUTDIR}:")
    for f in sorted(os.listdir(OUTDIR)):
        if f.endswith(('.png', '.pdf')):
            fpath = os.path.join(OUTDIR, f)
            size_kb = os.path.getsize(fpath) / 1024
            print(f"  {f}  ({size_kb:.0f} KB)")
