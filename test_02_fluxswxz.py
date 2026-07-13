import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from matplotlib.collections import LineCollection
import glob
import matplotlib
import copy

from fullmc import *

if __name__ == '__main__':
    work_dir = os.path.join('out', os.path.splitext(os.path.basename(__file__))[0])

    debug = 0
    # debug = 1

    # nx = 15
    # ny = 15
    # nz = 5
    nx = 28
    ny = 1
    nz = 28
    # nx = 25
    # ny = 1
    # nz = 16
    dx = 100.0
    dy = 100.0
    dz = 80.0

    transfermode = 0 # ICA
    # transfermode = 1 # 3D

    # swlw = 0 # LW
    # swlw = 1 # SW
    # source = 0 # TOA Direct
    # source = 1 # TOA Diffuse (Lambertian)
    # source = 2 # Volumetric within the atmosphere

    # source = 2 # Volumetric within the atmosphere
    # swlw = 0 # LW
    # nphoton = 2000

    source = 0 # TOA Direct
    # source = 1 # TOA Diffuse (Lambertian)
    # source = 4 # TOA detector
    swlw = 1 # SW
    # nphoton = 2e1
    nphoton = 3e2
    # nphoton = 5e3
    # nphoton = 2e4

    # solmu = 1.0
    # solmu = np.sqrt(2.0) / 2.0
    solmu = np.cos(np.radians(30.0))
    # solmu = 0.5
    solphi = 0.5 * np.pi

    viewmu = 1.0
    viewphi = 0.0

    # Ncpu = None
    Ncpu = 8

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
            homogenize=False,
            Ncpu=Ncpu,
            wrkdir=work_dir,
        )


    # kext = np.zeros((nx, ny, nz), dtype=np.float64)
    # kext[:, :, :] = 10.**(-4.5)
    # kext[-8:-3, 0:ny, -10:-3] = 10.**(-1.5)
    # kabs = np.zeros((nx, ny, nz), dtype=np.float64)
    # kabs[:, :, :] = 10.**(-4.6)
    # kabs[-8:-3, 0:ny, -10:-3] = 10.**(-2.)
    # gparam = np.zeros((nx, ny, nz), dtype=np.float64)
    # gparam[:, :, :] = 0.0001
    # gparam[-8:-3, 0:ny, -6:-3] = 0.01
    # bplnk = np.zeros((nx, ny, nz), dtype=np.float64)
    # bplnk[:, :, :] = 0.9
    # bgrnd = np.zeros((nx, ny), dtype=np.float64)
    # bgrnd[:, :] = 1.2
    # galb = np.zeros((nx, ny), dtype=np.float64)
    # galb[:, :] = 0.01

    # kext = np.zeros((nx, ny, nz), dtype=np.float64)
    # kext[:, :, :] = 1e-5
    # kext[4:8, 4:8, 2:4] = 5e-1
    # kabs = np.zeros((nx, ny, nz), dtype=np.float64)
    # kabs[:, :, :] = 1e-6
    # kabs[4:8, 4:8, 2:4] = 1e-6
    # gparam = np.zeros((nx, ny, nz), dtype=np.float64)
    # gparam[:, :, :] = 0.0001
    # gparam[4:8, 4:8, 2:4] = 0.85
    fmc.bplnk[:, :, :] = 0.0
    fmc.bgrnd[:, :] = 0.0
    fmc.galb[:, :] = 0.0

    # file_les = 'dat/mod/les_mod_01.txt'
    # files_les = glob.glob('dat/mod/lescrop_01_24-01_*.txt')
    # files_les = glob.glob('dat/mod/lescrop_01_24*.txt')
    # files_les = sorted(glob.glob('dat/mod/lescrop_01_*.txt'))
    # files_les = sorted(glob.glob('dat/mod/lescrop_02_*.txt'))
    files_les = sorted(glob.glob('dat/mod/lescrop_02_24*.txt'))

    data_plot = np.zeros((nx, nz), dtype=np.float64)
    data_plot2 = np.zeros((nx, nz), dtype=np.float64)
    

    for file_les in files_les:
        bnd = int(os.path.basename(file_les).split('_')[2].split('-')[0])
        if bnd < 15:
            continue
        
        print(file_les)
        fmc.read_atmtxt(file_les)

        absscale = np.exp(-np.mean(fmc.kabs[:, :, -1])*8.e3)
        print('Absorption scale factor: {:.6e}'.format(absscale))
        wgt = float(os.path.basename(file_les).split('_')[-2])
        fsol = float(os.path.basename(file_les).split('_')[-1].split('.txt')[0])

        fmc.transfermode = 0 # ICA
        fmc.run_mc()

        radirr1= fmc.read_result(kind='irr')
        radconv1 = fmc.read_result(kind='conv')


        fmc2 = copy.deepcopy(fmc)
        fmc2.transfermode = 1 # 3D
        fmc2.run_mc()
        radirr2 = fmc2.read_result(kind='irr')
        radconv2 = fmc2.read_result(kind='conv')

        indy = 0
        # data_plot = radirr1[:, indy, :, 0, 0, 0] + radirr1[:, indy, :, 1, 0, 0]
        # data_plot2 = radirr2[:, indy, :, 0, 0, 0] + radirr2[:, indy, :, 1, 0, 0]
        data_plot += (radirr1[:, indy, :, 0, 0, 0] + radirr1[:, indy, :, 1, 0, 0]) * fsol * wgt * absscale
        data_plot2 += (radirr2[:, indy, :, 0, 0, 0] + radirr2[:, indy, :, 1, 0, 0]) * fsol * wgt * absscale
        # heating_factor = 86400.0 / (1004.0 * 1.1)  # to convert W/m^3 to degC/day
        # data_plot += radconv1[:, indy, :, 0] * fsol * wgt * absscale * heating_factor
        # data_plot2 += radconv2[:, indy, :, 0] * fsol * wgt * absscale * heating_factor
        # heating_factor2 = 86400.0 / 4e6  # to convert W/m^2 to degC/day

    print('Mean surface irradiance (ICA): {:.6e} W/m^2'.format(np.mean(data_plot)))
    print('Mean surface irradiance (3D RT): {:.6e} W/m^2'.format(np.mean(data_plot2)))

    with open('dat/mod/lescrop_02_24-01_1.000e+00_4.409e+02.txt', 'r') as fh:
        line = fh.readline()
        nx_file, ny_file, nz_file = [int(x) for x in line.strip().split()]
        assert nx_file == nx
        assert ny_file == ny
        assert nz_file == nz
        kext = np.zeros((nx, ny, nz), dtype=np.float64)
        kabs = np.zeros((nx, ny, nz), dtype=np.float64)
        gparam = np.zeros((nx, ny, nz), dtype=np.float64)
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    line = fh.readline()
                    parts = line.strip().split()
                    kext[ix, iy, iz] = float(parts[0])
                    kabs[ix, iy, iz] = float(parts[1])
                    gparam[ix, iy, iz] = float(parts[2])

    data_max = max(np.max(data_plot), np.max(data_plot2))

    mask = (kext[:, indy, :] >= 3e-4)  # shape (nx, nz)
    x_edges = np.linspace(0.0, nx * dx * 1e-3, nx)
    z_edges = np.linspace(0.0, nz * dz * 1e-3, nz)
    segs = []
    for i in range(nx):
        for j in range(nz):
            if not mask[i, j]:
                continue
            x0, x1 = x_edges[i], x_edges[i + 1]
            z0, z1 = z_edges[j], z_edges[j + 1]
            if i == 0 or not mask[i - 1, j]:
                segs.append([(x0, z0), (x0, z1)])
            if i == nx - 1 or not mask[i + 1, j]:
                segs.append([(x1, z0), (x1, z1)])
            if j == 0 or not mask[i, j - 1]:
                segs.append([(x0, z0), (x1, z0)])
            if j == nz - 1 or not mask[i, j + 1]:
                segs.append([(x0, z1), (x1, z1)])

    xx, zz = np.meshgrid(np.linspace(0.0, nx * dx, nx) + dx/2., np.linspace(0.0, nz * dz, nz) + dz/2.)



    norm = mcolors.Normalize(vmin=0., vmax=data_max)
    # norm = mcolors.LogNorm(vmin=data_max*1e-4, vmax=data_max)

    # fig = plt.figure(figsize=(7, 3))
    fig = plt.figure(figsize=(3.5, 3.5))
    # ax = fig.add_subplot(1, 2, 1)
    # ax = fig.add_axes([0.08, 0.15, 0.35, 0.75])
    ax = fig.add_axes([0.08, 0.50, 0.72, 0.4])
    m = ax.pcolormesh(xx * 1e-3, zz * 1e-3, data_plot.T, cmap='Blues_r', norm=norm)
    lc = LineCollection(segs, colors='k', linewidths=1.0)
    ax.add_collection(lc)
    # ax.set_aspect('equal')
    # ax.set_xlabel('X (km)')
    ax.set_ylabel('Z (km)')
    ax.set_xticklabels([])
    # ax.set_title('(a) ICA')
    ax.text(0.03, 0.95, '(a) ICA', transform=ax.transAxes, fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='square', facecolor='white', alpha=0.8, edgecolor='none'))

    # fig.colorbar(m, ax=ax, label='Irradiance (W/m²)', orientation='horizontal', pad=0.2, shrink=0.65)
    # ax_cb = fig.add_axes([0.87, 0.2, 0.02, 0.5])  # [left, bottom, width, height] in figure coords
    ax_cb = fig.add_axes([0.85, 0.2, 0.02, 0.5])  # [left, bottom, width, height] in figure coords
    cbar = fig.colorbar(m, cax=ax_cb, label='Downward irradiance (W/m²)', orientation='vertical')
    # cbar.formatter.set_useMathText(True)
    # ax2 = fig.add_subplot(1, 2, 2)
    # ax2 = fig.add_axes([0.48, 0.15, 0.35, 0.75])
    ax2 = fig.add_axes([0.08, 0.05, 0.72, 0.4])
    m2 = ax2.pcolormesh(xx * 1e-3, zz * 1e-3, data_plot2.T, cmap='Blues_r', norm=norm)
    lc = LineCollection(segs, colors='k', linewidths=1.0, label='Cloud')
    ax2.add_collection(lc)
    # ax2.set_aspect('equal')
    ax2.set_xlabel('X (km)')
    ax2.set_ylabel('Z (m)')
    # ax2.set_yticklabels([])
    # ax2.set_title('(b) 3D')
    ax2.text(0.03, 0.95, '(b) 3D', transform=ax2.transAxes, fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='square', facecolor='white', alpha=0.8, edgecolor='none'))

    # create a separate axis below the panels for the legend
    handles, labels = ax2.get_legend_handles_labels()
    ax_leg = fig.add_axes([0.82, 0.75, 0.12, 0.12])  # [left, bottom, width, height] in figure coords
    ax_leg.axis('off')
    patch = matplotlib.patches.Patch(facecolor='none', edgecolor='k', label=(labels[0] if labels else 'Cloud'))
    ax_leg.legend([patch], [patch.get_label()], loc='center left', ncol=1, frameon=False)

    fig.text(0.92, 0.0, 'Solar zenith angle: 30°', ha='center', va='top', fontsize=10)

    fig.savefig(f'{work_dir}/02_flux_xz_comparison.png', dpi=300, bbox_inches='tight')
