import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import pickle

from fullmc import *

if __name__ == '__main__':
    work_dir = os.path.join('out', os.path.splitext(os.path.basename(__file__))[0])

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
    nphoton = 1e2
    # nphoton = 5e3
    # nphoton = 2e4
    # nphoton = 5e4

    # solmu = 1.0
    # solmu = np.sqrt(2.0) / 2.0
    # solmu = np.cos(np.radians(25.0))
    solmu = np.cos(np.radians(45.0))
    # solmu = 0.5
    solphi = 0.5 * np.pi

    # viewmu = 1.0
    # viewphi = 0.0
    viewmu = np.cos(np.radians(45.0))
    viewphi = np.radians(-90.0)

    # Ncpu = None
    Ncpu = 8
    # Ncpu = 32

    # taua_values = np.logspace(-3, 0, 4)
    # taua_values = np.logspace(-3.2, 0.2, 18)
    # taua_values = np.linspace(0.01, 2, 4)
    taua_values = np.linspace(0.01, 2, 10)
    # taua_values = np.linspace(0.01, 2, 200)
    radimg1_arr = np.full((len(taua_values), nx, ny, 2), np.nan, dtype=np.float64)
    radimg2_arr = np.full((len(taua_values), nx, ny, 2), np.nan, dtype=np.float64)

    # for iabsv, taua_value in enumerate(taua_values):
    #     print('iabsv = %d/%d, taua_value = %.3e' % (iabsv+1, len(taua_values), taua_value))

    #     fmc = FullMC(
    #             nx=nx, ny=ny, nz=nz,
    #             dx=dx, dy=dy, dz=dz,
    #             transfermode=transfermode,
    #             source=source,
    #             swlw=swlw,
    #             solmu=solmu, solphi=solphi,
    #             viewmu=viewmu, viewphi=viewphi,
    #             nphoton=nphoton,
    #             debug=debug,
    #             Ncpu=Ncpu,
    #             wrkdir=work_dir,
    #         )

    #     ksca = np.zeros((nx, ny, nz), dtype=np.float64)
    #     kabs = np.zeros((nx, ny, nz), dtype=np.float64)
    #     gparam = np.zeros((nx, ny, nz), dtype=np.float64)
        
    #     ksca[:, :, :] = 1e-6
    #     ksca[nx//5*2:nx//5*3, 0:ny, 5:7] = 5e-1
    #     kabs[:, :, :] = taua_value / (dz * nz)
    #     gparam[:, :, :] = 0.0001
    #     gparam[nx//5*2:nx//5*3, 0:ny, 5:7] = 0.85

    #     fmc.ksca[:, :, :] = ksca
    #     fmc.kabs[:, :, :, 0] = kabs
    #     fmc.gparam[:, :, :] = gparam
    #     fmc.bplnk[:, :, :] = 0.0
    #     fmc.bgrnd[:, :] = 0.0
    #     fmc.galb[:, :] = 0.3

    #     # fmc.read_atmtxt('dat/mod/les_mod_02.txt')

    #     fmc.run_mc()
    #     radimg1 = fmc.read_result(kind='img')

    #     fmc2 = fmc
    #     fmc2.transfermode = 0 # ICA

    #     fmc2.run_mc()
    #     radimg2 = fmc2.read_result(kind='img')

    #     radimg1_arr[iabsv, :, :, :] = radimg1[:, :, :]
    #     radimg2_arr[iabsv, :, :, :] = radimg2[:, :, :]
    
    fmc = FullMC(
            nx=nx, ny=ny, nz=nz,
            dx=dx, dy=dy, dz=dz,
            ng=len(taua_values),
            transfermode=transfermode,
            source=source,
            swlw=swlw,
            solmu=solmu, solphi=solphi,
            viewmu=viewmu, viewphi=viewphi,
            nphoton=nphoton,
            debug=debug,
            homogenize=False,
            Ncpu=Ncpu,
            wrkdir=work_dir,
        )

    ksca = np.zeros((nx, ny, nz), dtype=np.float64)
    kabs = np.zeros((nx, ny, nz, len(taua_values)), dtype=np.float64)
    gparam = np.zeros((nx, ny, nz), dtype=np.float64)
    
    ksca[:, :, :] = 1e-6
    ksca[nx//5*2:nx//5*3, 0:ny, 5:7] = 5e-1
    kabs[:, :, :, :] = taua_values[None, None, None, :] / (dz * nz)
    gparam[:, :, :] = 0.0001
    gparam[nx//5*2:nx//5*3, 0:ny, 5:7] = 0.85

    # fmc.kext[:, :, :] = kabs + ksca
    fmc.kabs[:, :, :, :] = kabs
    fmc.ksca[:, :, :] = ksca
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

    radimg1_arr[:, :, :, :] = radimg1[:, :, :, :].transpose(2, 0, 1, 3)
    radimg2_arr[:, :, :, :] = radimg2[:, :, :, :].transpose(2, 0, 1, 3)
    
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
            'radimg1_arr': radimg1_arr,
            'radimg2_arr': radimg2_arr,
        }, f)

    taua_slant_values = taua_values * (1.0 / solmu + 1.0 / viewmu)
    lntrans1_arr = np.log(np.pi*radimg1_arr)
    lntrans2_arr = np.log(np.pi*radimg2_arr)
    for ix in range(nx):
        fig = plt.figure(figsize=(5, 7))
        ax = fig.add_subplot(2, 1, 1)
        # ax.set_xscale('log')
        # ax.set_xlabel('Absorption Coefficient (1/m)')
        ax.set_ylabel('ln(transmittance)')
        ax.plot(taua_slant_values, lntrans1_arr[:, ix, 0, 0], 'o-', label='3D', color='tab:blue', markersize=4)
        ax.plot(taua_slant_values, lntrans2_arr[:, ix, 0, 0], 's-', label='ICA', color='tab:orange', markersize=4)
        ax.set_title('Radiance vs Absorption Optical Depth at X={:.1f} km'.format(ix*dx*1e-3))
        ax.legend()
        ax.set_xlim(taua_slant_values[0] * 0.9, taua_slant_values[-1] * 1.1)
        # ax.set_ylim(0.0, None)
        ax2 = fig.add_subplot(2, 1, 2)
        # ax2.set_xscale('log')
        ax2.set_xlabel('Column absorption optical depth')
        ax2.set_ylabel('Ratio (3D/ICA)')
        ratio = lntrans1_arr[:, ix, 0, 0] / lntrans2_arr[:, ix, 0, 0]
        ax2.plot(taua_slant_values, ratio, 'o-', label='Ratio (3D/ICA)', color='tab:green', markersize=4)
        # ax2.axhline(1.0, color='gray', linestyle='--')
        ax2.set_xlim(taua_slant_values[0] * 0.9, taua_slant_values[-1] * 1.1)
        # ax2.set_ylim(0.0, None)
        fig.tight_layout()
        fig.savefig(f'{work_dir}/01_radtoa_radiance_vs_absorption_x{ix:02d}.png'.format(ix), dpi=300, bbox_inches='tight')
        plt.close(fig)
    
    plot_data = np.pi*radimg1_arr[0, :, 0, 0]
    plot_data2 = np.pi*radimg2_arr[0, :, 0, 0]

    fig = plt.figure(figsize=(5, 3.2))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(np.arange(nx) * dx * 1e-3, plot_data, 'o-', label='3D', color='tab:blue', markersize=4)
    ax.plot(np.arange(nx) * dx * 1e-3, plot_data2, 's--', label='ICA', color='tab:orange', markersize=4)
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Reflectance')
    ax.set_title('Reflectance vs X at Absorption Optical Depth = {:.3e}'.format(taua_values[0]))
    ax.legend()
    fig.savefig(f'{work_dir}/01_radtoa_radiance_vs_x.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    plot_data1 = fmc.ksca[:, 0, :]
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
    