import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import pickle

import h5py

from fullmc import *

if __name__ == '__main__':
    work_dir = os.path.join('out', os.path.splitext(os.path.basename(__file__))[0])

    run = True
    # run = False

    debug = 0
    # debug = 1

    # nx = 30
    # ny = 30
    # nz = 6
    nx = 25
    ny = 1
    nz = 25
    dx = 1000.0
    dy = 1000.0
    dz = 1000.0

    # transfermode = 0 # ICA
    transfermode = 1 # 3D

    # swlw = 0 # LW
    # swlw = 1 # SW
    # source = 0 # TOA Direct
    # source = 1 # TOA Diffuse (Lambertian)
    # source = 2 # Volumetric within the atmosphere

    # source = 2 # Volumetric within the atmosphere
    # swlw = 0 # LW
    # nphoton = 2000

    # source = 0 # TOA Direct
    # source = 1 # TOA Diffuse (Lambertian)
    source = 4 # TOA detector
    swlw = 1 # SW
    # nphoton = 1
    # nphoton = 2e1
    # nphoton = 1e2
    nphoton = 5e3
    # nphoton = 2e4
    # nphoton = 5e4

    # solmu = 1.0
    # solmu = np.sqrt(2.0) / 2.0
    solmu = np.cos(np.radians(25.0))
    # solmu = np.cos(np.radians(30.0))
    # solmu = 0.5
    solphi = 0.5 * np.pi

    viewmu = 1.0
    viewphi = 0.0
    # viewmu = np.cos(np.radians(15.0))
    # viewphi = np.radians(-90.0)

    # Ncpu = None
    Ncpu = 8
    # Ncpu = 32

    # taua_values = np.logspace(-3, 0, 4)
    # taua_values = np.logspace(-3.2, 0.2, 18)
    # taua_values = np.linspace(0.01, 2, 4)
    ###

    with h5py.File('large_dat/o2band_benchmark.h5', 'r') as f:
        wvl = f['o2a/wvl'][:]
        taua_values_wv = f['o2a/optical_thickness/o2_column'][:]
        taua_z_wv = f['o2a/optical_thickness/o2_layer'][0:nz, :]
        tau_rayleigh_z_wv = f['o2a/optical_thickness/rayleigh_layer'][0:nz, :]
    
    index_sort_ = np.argsort(taua_values_wv)
    ratio = 0.9
    index_sort = index_sort_[:int(len(index_sort_)*ratio)] # only use the first 3/4 of the sorted values to avoid too large absorption optical depth
    wvl_all = wvl[index_sort]
    taua_values_all = taua_values_wv[index_sort]
    taua_z_all = taua_z_wv[:, index_sort]
    taur_z_all = tau_rayleigh_z_wv[:, index_sort]

    index_sparse = np.arange(0, len(taua_values_all), 1000)

    wvls = wvl_all[index_sparse]
    taua_values = taua_values_all[index_sparse]
    taua_z = taua_z_all[:, index_sparse]
    taur_z = taur_z_all[:, index_sparse]


    radimg1_arr = np.full((len(taua_values), nx, ny, 2), np.nan, dtype=np.float64)
    radimg2_arr = np.full((len(taua_values), nx, ny, 2), np.nan, dtype=np.float64)

    if run:
        for iabsv, taua_value in enumerate(taua_values):
            print('iabsv = %d/%d, taua_value = %.3e' % (iabsv+1, len(taua_values), taua_value))

            fmc = FullMC(
                    nx=nx, ny=ny, nz=nz,
                    dx=dx, dy=dy, dz=dz,
                    transfermode=transfermode,
                    source=source,
                    swlw=swlw,
                    solmu=solmu, solphi=solphi,
                    viewmu=viewmu, viewphi=viewphi,
                    nphoton=nphoton,
                    debug=debug,
                    Ncpu=Ncpu,
                    wrkdir=work_dir,
                )

            ksca = np.zeros((nx, ny, nz), dtype=np.float64)
            kabs = np.zeros((nx, ny, nz), dtype=np.float64)
            gparam = np.zeros((nx, ny, nz), dtype=np.float64)
            
            # ksca[:, :, :] = 1e-6
            # ksca[nx//5*2:nx//5*3, 0:ny, 4:6] = 5e-1
            # kabs[:, :, :] = taua_value / (dz * nz)
            # gparam[:, :, :] = 0.0001
            # gparam[nx//5*2:nx//5*3, 0:ny, 4:6] = 0.85
            ksca[:, :, :] = taur_z[:, iabsv][np.newaxis, np.newaxis, :] / dz
            ksca[nx//5*2:nx//5*3, 0:ny, 4:6] += 1e-1
            kabs[:, :, :] = taua_z[:, iabsv][np.newaxis, np.newaxis, :] / dz
            gparam[:, :, :] = 0.0001
            gparam[nx//5*2:nx//5*3, 0:ny, 4:6] = 0.85

            fmc.kext[:, :, :] = kabs + ksca
            fmc.kabs[:, :, :] = kabs
            fmc.gparam[:, :, :] = gparam
            fmc.bplnk[:, :, :] = 0.0
            fmc.bgrnd[:, :] = 0.0
            fmc.galb[:, :] = 0.3

            # fmc.read_atmtxt('dat/mod/les_mod_02.txt')

            fmc.run_mc()
            radimg1 = fmc.read_result(kind='img')

            fmc2 = fmc
            fmc2.transfermode = 0 # ICA

            fmc2.run_mc()
            radimg2 = fmc2.read_result(kind='img')

            radimg1_arr[iabsv, :, :, :] = radimg1[:, :, :]
            radimg2_arr[iabsv, :, :, :] = radimg2[:, :, :]
        
        pkl_dir = os.path.join(work_dir, 'pkl')
        if not os.path.exists(pkl_dir):
            os.makedirs(pkl_dir)
        with open(os.path.join(pkl_dir, 'pkl_01_radtoa_abs.pkl'), 'wb') as f:
            pickle.dump({
                'nx': nx, 'ny': ny, 'nz': nz,
                'dx': dx, 'dy': dy, 'dz': dz,
                'transfermode': transfermode,
                'source': source,
                'swlw': swlw,
                'solmu': solmu, 'solphi': solphi,
                'viewmu': viewmu, 'viewphi': viewphi,
                'nphoton': nphoton,
                'Ncpu': Ncpu,
                'taua_values': taua_values,
                'taua_z': taua_z,
                'taur_z': taur_z,
                'wvl': wvl,
                'radimg1_arr': radimg1_arr,
                'radimg2_arr': radimg2_arr,
                'fmc': fmc,
            }, f)
    else:
        pkl_dir = os.path.join(work_dir, 'pkl')
        with open(os.path.join(pkl_dir, 'pkl_01_radtoa_abs.pkl'), 'rb') as f:
            data = pickle.load(f)
            nx = data['nx']
            ny = data['ny']
            nz = data['nz']
            dx = data['dx']
            dy = data['dy']
            dz = data['dz']
            transfermode = data['transfermode']
            source = data['source']
            swlw = data['swlw']
            solmu = data['solmu']
            solphi = data['solphi']
            viewmu = data['viewmu']
            viewphi = data['viewphi']
            nphoton = data['nphoton']
            Ncpu = data['Ncpu']
            taua_values = data['taua_values']
            taua_z = data['taua_z']
            taur_z = data['taur_z']
            wvl = data['wvl']
            radimg1_arr = data['radimg1_arr']
            radimg2_arr = data['radimg2_arr']
            fmc = data['fmc']
    taua_slant_values = taua_values * (1.0 / solmu + 1.0 / viewmu)
    lntrans1_arr = np.log(np.pi*radimg1_arr)
    lntrans2_arr = np.log(np.pi*radimg2_arr)
    for ix in range(nx):
        x1, y1 = taua_slant_values, lntrans1_arr[:, ix, 0, 0]
        x2, y2 = taua_slant_values, lntrans2_arr[:, ix, 0, 0]
        try:
            coeff1 = np.polyfit(x1[x1 < 5.0], y1[x1 < 5.0], 2)
            coeff2 = np.polyfit(x2[x2 < 5.0], y2[x2 < 5.0], 2)
        except:
            continue
        x_arr = np.linspace(0.0, 5.0, 100)
        y_arr1 = np.polyval(coeff1, x_arr)
        y_arr2 = np.polyval(coeff2, x_arr)
        fig = plt.figure(figsize=(5, 7))
        ax = fig.add_subplot(2, 1, 1)
        ax.set_ylabel('Slant absorption optical depth')
        ax.set_ylabel('ln(transmittance)')
        ax.scatter(taua_slant_values, lntrans1_arr[:, ix, 0, 0], s=4, marker='o', label='3D', color='tab:blue')
        ax.scatter(taua_slant_values, lntrans2_arr[:, ix, 0, 0], s=4, marker='s', label='ICA', color='tab:orange')
        ax.plot(x_arr, y_arr1, '--', color='tab:blue', label='3D fit (%.1e x^2 + %.1e x + %.1e)' % (coeff1[0], coeff1[1], coeff1[2]))
        ax.plot(x_arr, y_arr2, '--', color='tab:orange', label='ICA fit (%.1e x^2 + %.1e x + %.1e)' % (coeff2[0], coeff2[1], coeff2[2]))
        ax.set_title('Radiance vs Absorption Optical Depth at X={:.1f} km'.format(ix*dx*1e-3))
        ax.legend()
        ax.set_xlim(0.0, 5.0)
        ax.set_ylim(min(np.min(y_arr1[x_arr < 5.0]), np.min(y_arr2[x_arr < 5.0])), 0.0)
        # ax.set_xlim(0.0, taua_slant_values.max() * 1.05)
        ax2 = fig.add_subplot(2, 1, 2)
        # ax2.set_xlabel('Column absorption optical depth')
        # ax2.set_ylabel('Ratio (3D/ICA)')
        # ratio = lntrans1_arr[:, ix, 0, 0] / lntrans2_arr[:, ix, 0, 0]
        # ax2.plot(taua_slant_values, ratio, 'o-', label='Ratio (3D/ICA)', color='tab:green', markersize=4)
        # ax2.set_xlim(taua_slant_values[0] * 0.9, taua_slant_values[-1] * 1.1)
        ax2.set_xlabel('Wavelength (nm)')
        ax2.set_ylabel('Radiance (W/m^2/nm/sr)')
        ax2.scatter(wvls, radimg1_arr[:, ix, 0, 0], label='3D', color='tab:blue', s=4)
        ax2.scatter(wvls, radimg2_arr[:, ix, 0, 0], label='ICA', color='tab:orange', s=4)
        ax2.set_xlim(wvl[0], wvl[-1])
        ax2.legend()
        fig.tight_layout()
        fig.savefig(f'{work_dir}/01_radtoa_radiance_vs_absorption_x{ix:02d}.png'.format(ix), dpi=300, bbox_inches='tight')
        plt.close(fig)
    
    for iabs in range(len(taua_values)):
        plot_data = np.pi*radimg1_arr[iabs, :, 0, 0]
        plot_data2 = np.pi*radimg2_arr[iabs, :, 0, 0]

        fig = plt.figure(figsize=(5, 3.2))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(np.arange(nx) * dx * 1e-3, plot_data, 'o-', label='3D', color='tab:blue', markersize=4)
        ax.plot(np.arange(nx) * dx * 1e-3, plot_data2, 's--', label='ICA', color='tab:orange', markersize=4)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Reflectance')
        ax.set_title('Reflectance vs X at wvl = {:.5f} nm, Absorption Optical Depth = {:.3e}'.format(wvls[iabs], taua_values[iabs]))
        ax.legend()
        fig.savefig(f'{work_dir}/01_radtoa_radiance_vs_x_{iabs:05d}.png', dpi=300, bbox_inches='tight')
        plt.close(fig)

    plot_data1 = (fmc.kext[:, 0, :] - fmc.kabs[:, 0, :])
    xx, zz = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, nz * 80.0, nz))

    fig = plt.figure(figsize=(5, 3.2))
    ax = fig.add_subplot(1, 2, 1)
    ax.set_aspect('equal')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Z (km)')
    pcm = ax.pcolormesh(xx*1e-3, zz*1e-3, plot_data1.T, cmap='cividis', vmin=0., vmax=np.max(plot_data1))
    ax.set_title('Scattering coefficient')
    fig.colorbar(pcm, ax=ax, orientation='horizontal', label='(1/m)', pad=0.22, shrink=0.45)
    fig.savefig(f'{work_dir}/01_radtoa_optical_depth.png', dpi=300, bbox_inches='tight')
    