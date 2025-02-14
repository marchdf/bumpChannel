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
plt.rc("text", usetex=True)
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]

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
if __name__ == "__main__":

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(description="A simple plot tool")
    parser.add_argument("-s", "--show", help="Show the plots", action="store_true")
    args = parser.parse_args()
    tol = 1e-14

    # ========================================================================
    # NASA CFL3D output

    # Skin friction at different slices
    fname = os.path.join(os.path.abspath("nasa_data"), "cf_cfl3d.dat")
    df = pd.read_csv(fname)

    subdf = df[np.fabs(df["xslice"] - 0.75) < tol]
    plt.figure(0)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    subdf = df[np.fabs(df["xslice"] - 0.6321975) < tol]
    plt.figure(1)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    subdf = df[np.fabs(df["xslice"] - 0.8678025) < tol]
    plt.figure(2)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    # Other coefficients
    fname = os.path.join(os.path.abspath("nasa_data"), "coeffs_cfl3d.dat")
    df = pd.read_csv(fname)

    plt.figure(3)
    plt.semilogx(
        df["h"],
        df["cd"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    plt.figure(4)
    plt.semilogx(
        df["h"],
        df["cdp"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    plt.figure(5)
    plt.semilogx(
        df["h"],
        df["cdv"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    plt.figure(6)
    plt.semilogx(
        df["h"],
        df["cl"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    # wall cf and cp
    fname = os.path.join(os.path.abspath("nasa_data"), "wall_coeffs_cfl3d.dat")
    df = pd.read_csv(fname)

    plt.figure(7)
    p = plt.plot(df["x"], df["cf"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    plt.figure(8)
    p = plt.plot(df["x"], df["cp"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    # profiles
    fname = os.path.join(os.path.abspath("nasa_data"), "mut_profile_cfl3d.dat")
    df = pd.read_csv(fname, delim_whitespace=True)

    plt.figure(9)
    p = plt.plot(df["mut"], df["y"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    fname = os.path.join(os.path.abspath("nasa_data"), "k_sdr_profile_cfl3d.dat")
    df = pd.read_csv(fname, delim_whitespace=True)

    plt.figure(10)
    p = plt.loglog(df["tke"], df["y"] - 0.05, lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    plt.figure(11)
    p = plt.loglog(df["sdr"], df["y"] - 0.05, lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    # ======================================================================
    # NASA FUN3D output

    # Skin friction at different slices
    fname = os.path.join(os.path.abspath("nasa_data"), "cf_fun3d.dat")
    df = pd.read_csv(fname)

    subdf = df[np.fabs(df["xslice"] - 0.75) < tol]
    plt.figure(0)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    subdf = df[np.fabs(df["xslice"] - 0.6321975) < tol]
    plt.figure(1)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    subdf = df[np.fabs(df["xslice"] - 0.8678025) < tol]
    plt.figure(2)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    # Other coefficients
    fname = os.path.join(os.path.abspath("nasa_data"), "coeffs_fun3d.dat")
    df = pd.read_csv(fname)

    plt.figure(3)
    plt.semilogx(
        df["h"],
        df["cd"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    plt.figure(4)
    plt.semilogx(
        df["h"],
        df["cdp"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    plt.figure(5)
    plt.semilogx(
        df["h"],
        df["cdv"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    plt.figure(6)
    plt.semilogx(
        df["h"],
        df["cl"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    # wall cf and cp
    fname = os.path.join(os.path.abspath("nasa_data"), "wall_coeffs_fun3d.dat")
    df = pd.read_csv(fname)

    plt.figure(7)
    p = plt.plot(df["x"], df["cf"], lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])

    plt.figure(8)
    p = plt.plot(df["x"], df["cp"], lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])

    # profiles
    fname = os.path.join(os.path.abspath("nasa_data"), "mut_profile_fun3d.dat")
    df = pd.read_csv(fname, delim_whitespace=True)

    plt.figure(9)
    p = plt.plot(df["mut"], df["y"], lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])

    fname = os.path.join(os.path.abspath("nasa_data"), "k_sdr_profile_fun3d.dat")
    df = pd.read_csv(fname, delim_whitespace=True)

    plt.figure(10)
    p = plt.loglog(df["tke"], df["y"] - 0.05, lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])

    plt.figure(11)
    p = plt.loglog(df["sdr"], df["y"] - 0.05, lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])

    # ======================================================================
    # Nalu output

    # Skin friction at different slices
    fname = "cf.dat"
    df = pd.read_csv(fname, index_col=0)
    df.drop("89x41", inplace=True)

    subdf = df[np.fabs(df["xslice"] - 0.75) < tol]
    plt.figure(0)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[2],
        marker=markertype[1],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    subdf = df[np.fabs(df["xslice"] - 0.6321975) < tol]
    plt.figure(1)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[2],
        marker=markertype[1],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    subdf = df[np.fabs(df["xslice"] - 0.8678025) < tol]
    plt.figure(2)
    plt.semilogx(
        subdf["h"],
        subdf["cf"],
        lw=2,
        color=cmap[2],
        marker=markertype[1],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    # Other coefficients
    fname = "coeffs.dat"
    df = pd.read_csv(fname, index_col=0)
    df.drop("89x41", inplace=True)

    plt.figure(3)
    plt.semilogx(
        df["h"],
        df["cd"],
        lw=2,
        color=cmap[2],
        marker=markertype[1],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    plt.figure(4)
    plt.semilogx(
        df["h"],
        df["cdp"],
        lw=2,
        color=cmap[2],
        marker=markertype[1],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    plt.figure(5)
    plt.semilogx(
        df["h"],
        df["cdv"],
        lw=2,
        color=cmap[2],
        marker=markertype[1],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    plt.figure(6)
    plt.semilogx(
        df["h"],
        df["cl"],
        lw=2,
        color=cmap[2],
        marker=markertype[1],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    # wall cf and cp
    fdir = "1409x641"
    df = pd.read_csv(os.path.join(fdir, "results", "wall_coeffs.dat"))

    plt.figure(7)
    p = plt.plot(df["x"], df["cf"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    plt.figure(8)
    p = plt.plot(df["x"], df["cp"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    # profiles
    df = pd.read_csv(os.path.join(fdir, "results", "profiles.dat"))

    plt.figure(9)
    p = plt.plot(df["nut"], df["z"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    plt.figure(10)
    p = plt.loglog(df["tke"], df["z"] - 0.05, lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    plt.figure(11)
    p = plt.loglog(df["sdr"], df["z"] - 0.05, lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    # ========================================================================
    # Format the plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_f$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cf_0.75.pdf", format="pdf")
    plt.savefig("cf_0.75.png", format="png")

    plt.figure(1)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_f$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cf_0.63.pdf", format="pdf")
    plt.savefig("cf_0.63.png", format="png")

    plt.figure(2)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_f$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cf_0.86.pdf", format="pdf")
    plt.savefig("cf_0.86.png", format="png")

    plt.figure(3)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_d$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cd.pdf", format="pdf")
    plt.savefig("cd.png", format="png")

    plt.figure(4)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_{d,p}$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cdp.pdf", format="pdf")
    plt.savefig("cdp.png", format="png")

    plt.figure(5)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_{d,v}$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cdv.pdf", format="pdf")
    plt.savefig("cdv.png", format="png")

    plt.figure(6)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_l$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cl.pdf", format="pdf")
    plt.savefig("cl.png", format="png")

    plt.figure(7)
    ax = plt.gca()
    plt.xlabel(r"$x~[m]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_f$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.xlim([0.0, 1.5])
    plt.ylim([0.0, 0.008])
    plt.tight_layout()
    plt.savefig("wall_cf.pdf", format="pdf")
    plt.savefig("wall_cf.png", format="png")

    plt.figure(8)
    ax = plt.gca()
    plt.xlabel(r"$x~[m]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_p$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.xlim([0.0, 1.5])
    plt.ylim([-0.8, 0.4])
    plt.tight_layout()
    plt.savefig("wall_cp.pdf", format="pdf")
    plt.savefig("wall_cp.png", format="png")

    plt.figure(9)
    ax = plt.gca()
    plt.xlabel(r"$\mu_t / \mu_\infty$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y~[m]$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.ylim([0.05, 0.065])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("nut.pdf", format="pdf")
    plt.savefig("nut.png", format="png")

    plt.figure(10)
    ax = plt.gca()
    plt.xlabel(r"$k / a_\infty^2$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y - 0.05~[m]$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.xlim([1e-11, 1e-3])
    plt.ylim([1e-7, 1e-1])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("tke.pdf", format="pdf")
    plt.savefig("tke.png", format="png")

    plt.figure(11)
    ax = plt.gca()
    plt.xlabel(
        r"$\omega \mu_\infty / (\rho_\infty a_\infty^2)$",
        fontsize=22,
        fontweight="bold",
    )
    plt.ylabel(r"$y - 0.05~[m]$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.xlim([1e-9, 1e0])
    plt.ylim([1e-7, 1e-1])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("sdr.pdf", format="pdf")
    plt.savefig("sdr.png", format="png")

    if args.show:
        plt.show()
