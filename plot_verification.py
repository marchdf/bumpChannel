#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ========================================================================
#
# Some defaults variables
#
# ========================================================================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Times')
cmap_med = ['#F15A60', '#7AC36A', '#5A9BD4', '#FAA75B',
            '#9E67AB', '#CE7058', '#D77FB4', '#737373']
cmap = ['#EE2E2F', '#008C48', '#185AA9', '#F47D23',
        '#662C91', '#A21D21', '#B43894', '#010202']
dashseq = [(None, None), [10, 5], [10, 4, 3, 4], [
    3, 3], [10, 4, 3, 4, 3, 4], [3, 3], [3, 3]]
markertype = ['s', 'd', 'o', 'p', 'h']
# scale is used to adjust issues in pp.py area lengthscale without having
# to rerun
scale = 4.0/3.0

# ========================================================================
#
# Function definitions
#
# ========================================================================


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == '__main__':

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='A simple plot tool')
    parser.add_argument(
        '-s', '--show', help='Show the plots', action='store_true')
    args = parser.parse_args()

    # ========================================================================
    # NASA CFL3D output

    # Skin friction at different slices
    fname = os.path.join(os.path.abspath('nasa_data'), 'cf_cfl3d.dat')
    df = pd.read_csv(fname)

    subdf = df[np.fabs(df['xslice'] - 0.75) < 1e-14]
    plt.figure(0)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[0],
                 marker=markertype[0], mec=cmap[0], mfc=cmap[0], ms=10,
                 label='NASA CFL3D')

    subdf = df[np.fabs(df['xslice'] - 0.6321975) < 1e-14]
    plt.figure(1)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[0],
                 marker=markertype[0], mec=cmap[0], mfc=cmap[0], ms=10,
                 label='NASA CFL3D')

    subdf = df[np.fabs(df['xslice'] - 0.8678025) < 1e-14]
    plt.figure(2)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[0],
                 marker=markertype[0], mec=cmap[0], mfc=cmap[0], ms=10,
                 label='NASA CFL3D')

    # Other coefficients
    fname = os.path.join(os.path.abspath('nasa_data'), 'coeffs_cfl3d.dat')
    df = pd.read_csv(fname)

    plt.figure(3)
    plt.semilogx(df['h'], df['cd'], lw=2, color=cmap[0], marker=markertype[0],
                 mec=cmap[0], mfc=cmap[0], ms=10, label='NASA CFL3D')

    plt.figure(4)
    plt.semilogx(df['h'], df['cdp'], lw=2, color=cmap[0], marker=markertype[0],
                 mec=cmap[0], mfc=cmap[0], ms=10, label='NASA CFL3D')

    plt.figure(5)
    plt.semilogx(df['h'], df['cdv'], lw=2, color=cmap[0], marker=markertype[0],
                 mec=cmap[0], mfc=cmap[0], ms=10, label='NASA CFL3D')

    plt.figure(6)
    plt.semilogx(df['h'], df['cl'], lw=2, color=cmap[0], marker=markertype[0],
                 mec=cmap[0], mfc=cmap[0], ms=10, label='NASA CFL3D')

    # wall cf and cp
    fname = os.path.join(os.path.abspath('nasa_data'), 'wall_coeffs_cfl3d.dat')
    df = pd.read_csv(fname)

    plt.figure(7)
    p = plt.plot(df['x'], df['cf'], lw=2, color=cmap[0], label='NASA CFL3D')
    p[0].set_dashes(dashseq[0])

    plt.figure(8)
    p = plt.plot(df['x'], df['cp'], lw=2, color=cmap[0], label='NASA CFL3D')
    p[0].set_dashes(dashseq[0])

    # ======================================================================
    # NASA FUN3D output

    # Skin friction at different slices
    fname = os.path.join(os.path.abspath('nasa_data'), 'cf_fun3d.dat')
    df = pd.read_csv(fname)

    subdf = df[np.fabs(df['xslice'] - 0.75) < 1e-14]
    plt.figure(0)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[1],
                 marker=markertype[1], mec=cmap[1], mfc=cmap[1], ms=10,
                 label='NASA FUN3D')

    subdf = df[np.fabs(df['xslice'] - 0.6321975) < 1e-14]
    plt.figure(1)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[1],
                 marker=markertype[1], mec=cmap[1], mfc=cmap[1], ms=10,
                 label='NASA FUN3D')

    subdf = df[np.fabs(df['xslice'] - 0.8678025) < 1e-14]
    plt.figure(2)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[1],
                 marker=markertype[1], mec=cmap[1], mfc=cmap[1], ms=10,
                 label='NASA FUN3D')

    # Other coefficients
    fname = os.path.join(os.path.abspath('nasa_data'), 'coeffs_fun3d.dat')
    df = pd.read_csv(fname)

    plt.figure(3)
    plt.semilogx(df['h'], df['cd'], lw=2, color=cmap[1], marker=markertype[1],
                 mec=cmap[1], mfc=cmap[1], ms=10, label='NASA FUN3D')

    plt.figure(4)
    plt.semilogx(df['h'], df['cdp'], lw=2, color=cmap[1], marker=markertype[1],
                 mec=cmap[1], mfc=cmap[1], ms=10, label='NASA FUN3D')

    plt.figure(5)
    plt.semilogx(df['h'], df['cdv'], lw=2, color=cmap[1], marker=markertype[1],
                 mec=cmap[1], mfc=cmap[1], ms=10, label='NASA FUN3D')

    plt.figure(6)
    plt.semilogx(df['h'], df['cl'], lw=2, color=cmap[1], marker=markertype[1],
                 mec=cmap[1], mfc=cmap[1], ms=10, label='NASA FUN3D')

    # wall cf and cp
    fname = os.path.join(os.path.abspath('nasa_data'), 'wall_coeffs_fun3d.dat')
    df = pd.read_csv(fname)

    plt.figure(7)
    p = plt.plot(df['x'], df['cf'], lw=2, color=cmap[1], label='NASA FUN3D')
    p[0].set_dashes(dashseq[1])

    plt.figure(8)
    p = plt.plot(df['x'], df['cp'], lw=2, color=cmap[1], label='NASA FUN3D')
    p[0].set_dashes(dashseq[1])

    # ======================================================================
    # Nalu output

    # Skin friction at different slices
    fname = 'cf.dat'
    df = pd.read_csv(fname)

    subdf = df[np.fabs(df['xslice'] - 0.75) < 1e-14]
    plt.figure(0)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[2],
                 marker=markertype[1], mec=cmap[2], mfc=cmap[2], ms=10,
                 label='Nalu')

    subdf = df[np.fabs(df['xslice'] - 0.6321975) < 1e-14]
    plt.figure(1)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[2],
                 marker=markertype[1], mec=cmap[2], mfc=cmap[2], ms=10,
                 label='Nalu')

    subdf = df[np.fabs(df['xslice'] - 0.8678025) < 1e-14]
    plt.figure(2)
    plt.semilogx(subdf['h'], subdf['cf'], lw=2, color=cmap[2],
                 marker=markertype[1], mec=cmap[2], mfc=cmap[2], ms=10,
                 label='Nalu')

    # Other coefficients
    fname = 'coeffs.dat'
    df = pd.read_csv(fname)

    plt.figure(3)
    plt.semilogx(df['h'], df['cd']*scale, lw=2, color=cmap[2], marker=markertype[1],
                 mec=cmap[2], mfc=cmap[2], ms=10, label='Nalu')

    plt.figure(4)
    plt.semilogx(df['h'], df['cdp']*scale, lw=2, color=cmap[2], marker=markertype[1],
                 mec=cmap[2], mfc=cmap[2], ms=10, label='Nalu')

    plt.figure(5)
    plt.semilogx(df['h'], df['cdv']*scale, lw=2, color=cmap[2], marker=markertype[1],
                 mec=cmap[2], mfc=cmap[2], ms=10, label='Nalu')

    plt.figure(6)
    plt.semilogx(df['h'], df['cl']*scale, lw=2, color=cmap[2], marker=markertype[1],
                 mec=cmap[2], mfc=cmap[2], ms=10, label='Nalu')

    # wall cf and cp
    df = pd.read_csv('1409x641/results/wall_coeffs.dat')

    plt.figure(7)
    p = plt.plot(df['x'], df['cf'], lw=2, color=cmap[2], label='Nalu')
    p[0].set_dashes(dashseq[1])

    plt.figure(8)
    p = plt.plot(df['x'], df['cp'], lw=2, color=cmap[2], label='Nalu')
    p[0].set_dashes(dashseq[1])

    # ========================================================================
    # Format the plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_f$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('cf_0.75.pdf', format='pdf')
    plt.savefig('cf_0.75.png', format='png')

    plt.figure(1)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_f$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('cf_0.63.pdf', format='pdf')
    plt.savefig('cf_0.63.png', format='png')

    plt.figure(2)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_f$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('cf_0.86.pdf', format='pdf')
    plt.savefig('cf_0.86.png', format='png')

    plt.figure(3)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_d$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('cd.pdf', format='pdf')
    plt.savefig('cd.png', format='png')

    plt.figure(4)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_{d,p}$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('cdp.pdf', format='pdf')
    plt.savefig('cdp.png', format='png')

    plt.figure(5)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_{d,v}$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('cdv.pdf', format='pdf')
    plt.savefig('cdv.png', format='png')

    plt.figure(6)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_l$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('cl.pdf', format='pdf')
    plt.savefig('cl.png', format='png')

    plt.figure(7)
    ax = plt.gca()
    plt.xlabel(r"$x~[m]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_f$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.xlim([0.0, 1.5])
    plt.ylim([0.0, 0.008])
    plt.tight_layout()
    plt.savefig('wall_cf.pdf', format='pdf')
    plt.savefig('wall_cf.png', format='png')

    plt.figure(8)
    ax = plt.gca()
    plt.xlabel(r"$x~[m]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_p$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.xlim([0.0, 1.5])
    plt.ylim([-0.8, 0.4])
    plt.tight_layout()
    plt.savefig('wall_cp.pdf', format='pdf')
    plt.savefig('wall_cp.png', format='png')

    if args.show:
        plt.show()
