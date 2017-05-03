#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import os
import numpy as np
import pandas as pd
import netCDF4 as nc
import glob
import yaml
import scipy.integrate as spi


# ========================================================================
#
# Function definitions
#
# ========================================================================
def parse_ic(fname):
    """Parse the Nalu yaml input file for the initial conditions"""
    with open(fname, 'r') as stream:
        try:
            dat = yaml.load(stream)
            u0 = float(dat['realms'][0]['initial_conditions']
                       [0]['value']['velocity'][0])
            rho0 = float(dat['realms'][0]['material_properties']
                         ['specifications'][0]['value'])
            mu = float(dat['realms'][0]['material_properties']
                       ['specifications'][1]['value'])

            return u0, rho0, mu

        except yaml.YAMLError as exc:
            print(exc)


# ========================================================================
def get_wall_values(enames):
    """Get wall values from Exodus file."""
    data = np.array([]).reshape(0, 3)
    lst = []

    # Loop on variables
    for ename in enames:
        dat = nc.Dataset(ename)

        ssn = ["%s" % nc.chartostring(ss)
               for ss in dat.variables['ss_names'][:]]
        vn = ["%s" % nc.chartostring(nn)
              for nn in dat.variables['name_nod_var'][:]]
        ### JAM: removed the "b'*'" from the begging of these arguments, 
        ###      not sure why they were there?
        idx_ss = ssn.index("bottomwall")
        idx_tw = vn.index("tau_wall")
        idx_p = vn.index("pressure")

        #######***************************************************************************##############
        ### JAM: This is the replacement code

        listSize = dat.variables['coordx'].size
        dArray = []
        # Loop through the dataset to extract the desired variables
        for i in range(0,listSize):
            # tau_wall should only be defined on the wall, so use that as a "flag" to
            # identify which are wall values
            if (dat.variables['vals_nod_var{0:d}'.format(idx_tw + 1)][-1][i] > 0.0):
                dArray.append([dat.variables['coordx'][i],dat.variables['coordz'][i],
                   dat.variables['vals_nod_var{0:d}'.format(idx_tw + 1)][-1][i],
                   dat.variables['vals_nod_var{0:d}'.format(idx_p + 1)][-1][i]])

        # To dataframe
        # Transpose the dataArray and append to list (assuming its non-empty)
        if (len(dArray) > 0):
            tdArray = map(list, zip(*dArray))
            df = pd.DataFrame(data=np.vstack((tdArray[0],tdArray[2],tdArray[3])).T,
                          columns=['x', 'tau_wall', 'pressure'])
            lst.append(df)
        #######***************************************************************************##############


        #######***************************************************************************##############
        ### JAM: Commented this out as it seems to not capture all the wall points,
        ###      I have replaced it with the above method for capturing wall values
        """        
        # Get wall x coordinates
        try:
            wall_elem_idx = dat.variables[
                'elem_ss{0:d}'.format(idx_ss + 1)][:] - 1
        except:
            continue
        wall_connect1 = dat.variables['connect1'][wall_elem_idx].flatten() - 1
        wall_coordx = dat.variables['coordx'][wall_connect1]
        wall_coordz = dat.variables['coordz'][wall_connect1]

        actual_idx = np.where(wall_coordz < 1e-14)
        wall_x = wall_coordx[actual_idx]
        wall_connect1 = wall_connect1[actual_idx]

        wall_x, unique_idx = np.unique(wall_x, return_index=True)
        wall_connect1 = wall_connect1[unique_idx]

        # Get wall variables
        if (wall_connect1.size > 0):
            tau_wall = dat.variables[
                'vals_nod_var{0:d}'.format(idx_tw + 1)][-1, wall_connect1]
            pressure = dat.variables[
                'vals_nod_var{0:d}'.format(idx_p + 1)][-1, wall_connect1]

        # To dataframe
        df = pd.DataFrame(data=np.vstack((wall_x, tau_wall, pressure)).T,
                          columns=['x', 'tau_wall', 'pressure'])
        lst.append(df)
        """
        #######***************************************************************************##############

    # Save
    dfw = pd.concat(lst, ignore_index=True)
    dfw = dfw.sort_values(by=['x'])
    dfw = dfw.drop_duplicates(subset='x')
    return dfw.reset_index()


# ========================================================================
def get_ux_front(enames):
    """Get ux on the front from Exodus file."""

    lst = []

    for ename in enames:
        dat = nc.Dataset(ename)

        ssn = ["%s" % nc.chartostring(ss)
               for ss in dat.variables['ss_names'][:]]
        vn = ["%s" % nc.chartostring(nn)
              for nn in dat.variables['name_nod_var'][:]]

        ### JAM : Removed b... see above note
        idx_front = ssn.index("front")
        idx_ux = vn.index("velocity_x")

        try:
            front_elem_idx = dat.variables[
                'elem_ss{0:d}'.format(idx_front + 1)][:] - 1
        except:
            continue

        front_connect1 = dat.variables['connect1'][
            front_elem_idx].flatten() - 1
        front_coordx = dat.variables['coordx'][front_connect1]
        front_coordy = dat.variables['coordy'][front_connect1]
        front_coordz = dat.variables['coordz'][front_connect1]

        actual_idx = np.where(front_coordy >= 0.0)
        front_x = front_coordx[actual_idx]
        front_y = front_coordy[actual_idx]
        front_z = front_coordz[actual_idx]
        front_connect1 = front_connect1[actual_idx]

        front_ux = dat.variables[
            'vals_nod_var{0:d}'.format(idx_ux + 1)][-1, front_connect1]

        colnames = ['x', 'z', 'ux']
        df = pd.DataFrame(data=np.vstack((front_x, front_z, front_ux)).T,
                          columns=colnames)
        lst.append(df)

    # Save
    dfu = pd.concat(lst, ignore_index=True)
    dfu.x[np.fabs(dfu['x']) < 1e-14] = 0.0  # make true zeros
    dfu = dfu.sort_values(by=['x', 'z'])
    dfu = dfu.drop_duplicates(subset=['x', 'z'])
    dfu = dfu[dfu['x'] >= 0]  # remove everything before the plate
    return dfu.reset_index()


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == '__main__':

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Postprocess Nalu data')
    args = parser.parse_args()

    # ========================================================================
    # Setup
    ppdirs = ['89x41', "177x81", "353x161"]# "705x321", "1409x641"]
    fdirs = [os.path.abspath(fdir) for fdir in ppdirs]
    rdirs = [os.path.join(fdir, 'results') for fdir in fdirs]
    ocname = os.path.join(os.path.abspath('.'), 'coeffs.dat')
    ocfname = os.path.join(os.path.abspath('.'), 'cf.dat')
    colnames = ['N2', 'h', 'cd', 'cdp', 'cdv', 'cl']
    dfc = pd.DataFrame(index=ppdirs, columns=colnames)
    dfcf = pd.DataFrame(index=ppdirs, columns=['ppdir', 'N2', 'h', 'xslice', 'cf', 'cp', 'cd', 'cl'])

    # ========================================================================
    # Post-process
    cnt = 0
    for ppdir, fdir, rdir in zip(ppdirs, fdirs, rdirs):
        print('Post-processing directory:', ppdir)
        fname = os.path.join(rdir, 'bumpChannel.dat')
        yname = os.path.join(fdir, 'bumpChannel.i')
        enames = glob.glob(os.path.join(rdir, 'bumpChannel.e*'))
        owname = os.path.join(rdir, 'wall_coeffs.dat')
        ouname = os.path.join(rdir, 'yp_up.dat')

        # Derived quantities
        L = 2.0
        W = 1.0
        area = L * W
        u0, rho0, mu = parse_ic(yname)
        dynPres = rho0 * 0.5 * u0 * u0

        # ---------------------------------------------
        # Get wall values, coefficients, etc
        dfw = get_wall_values(enames)

        # Get integrated wall values
        df = pd.read_csv(fname, delim_whitespace=True)

        # Calculate coefficients
        df['cl'] = (df['Fpy'] + df['Fvy']) / (dynPres * area)
        df['cd'] = (df['Fpx'] + df['Fvx']) / (dynPres * area)
        df['cdp'] = df['Fpx'] / (dynPres * area)
        df['cdv'] = df['Fvx'] / (dynPres * area)
        dfw['cf'] = dfw['tau_wall'] / dynPres
        dfw['cp'] = dfw['pressure'] / dynPres
        xslice1 = 0.75
        cf_slice1 = np.interp(xslice1, dfw['x'], dfw['cf'])
        cp_slice1 = np.interp(xslice1, dfw['x'], dfw['cp'])
        xslice2 = 0.6321075
        cf_slice2 = np.interp(xslice2, dfw['x'], dfw['cf'])
        cp_slice2 = np.interp(xslice2, dfw['x'], dfw['cp'])
        xslice3 = 0.8678025
        cf_slice3 = np.interp(xslice3, dfw['x'], dfw['cf'])
        cp_slice3 = np.interp(xslice3, dfw['x'], dfw['cp'])

        # ---------------------------------------------
        # Also calculate Re_theta and u+
        dfu = get_ux_front(enames)

        # Reshapes
        res = [int(s) for s in ppdir.split('x')]
        nx, nz = int(dfu.shape[0] / res[1]), res[1]
        x = dfu['x'].values.reshape((nx, nz))
        z = dfu['z'].values.reshape((nx, nz))
        ux = dfu['ux'].values.reshape((nx, nz))

        # Integrate to get theta and other quantities (and save)
        # The integral bound is where u is less than 99.5% of the freestream
        # velocity

        ##########**********JAM: FIXME: This part is commented out casue it is
        #                               causing issues with the bumpChannel, it
        #                               doens't necessarily seemed needed either?
        #theta = []
        #for i in range(ux.shape[0]):
        #    idx = ux[i, :] < 0.995 * u0
        #    theta.append(spi.simps(ux[i, idx] / u0 *
        #                           (1 - ux[i, idx] / u0), z[i, idx]))
        #
        #dfw['theta'] = theta
        #dfw['retheta'] = rho0 * u0 * dfw['theta'] / mu

        ## Get y+ and u+ at Re-theta = 10000 (or next best thing)
        #idx = np.argmin(np.fabs(dfw['retheta'] - 10000))
        #up = ux[idx, :] / np.sqrt(dfw['tau_wall'][idx] / rho0)
        #yp = z[idx, :] * np.sqrt(dfw['tau_wall'][idx] / rho0) / (mu / rho0)

        #odf = pd.DataFrame(data=np.vstack((z[idx, :], yp, up)).T,
        #                   columns=['z', 'yp', 'up'])
        ##########********** END FIXME

        # ---------------------------------------------
        # Save
        res = [float(s) for s in ppdir.split('x')]
        N2 = (res[0] - 1) * (res[1] - 1)
        h = np.sqrt(1. / N2)
        dfc.loc[ppdir] = [N2, h, df['cd'].iloc[-1],
                          df['cdp'].iloc[-1], df['cdv'].iloc[-1],
                          df['cl'].iloc[-1]]
        cnt = cnt + 1
        dfcf.loc[cnt] = [ppdir, N2, h, xslice1, cf_slice1, cp_slice1,
                           df['cd'].iloc[-1], df['cl'].iloc[-1]]
        cnt = cnt + 1
        dfcf.loc[cnt] = [ppdir, N2, h, xslice2, cf_slice2, cp_slice2,
                           df['cd'].iloc[-1], df['cl'].iloc[-1]]
        cnt = cnt + 1
        dfcf.loc[cnt] = [ppdir, N2, h, xslice3, cf_slice3, cp_slice3,
                           df['cd'].iloc[-1], df['cl'].iloc[-1]]
        dfw.to_csv(owname, index=False)

        # JAM: FIXME: commented out due to association with above commented out part
        #odf.to_csv(ouname, index=False)

    # Save the coefficients in a convenient table
    dfc.to_csv(ocname)
    dfcf.sort_values(by=['xslice', 'h'])
    dfcf.to_csv(ocfname)
