import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
import os

from fullmc import *

if __name__ == '__main__':
    work_dir = os.path.join('out', os.path.splitext(os.path.basename(__file__))[0])

    debug = 0
    # debug = 1

    # nx = 30
    # ny = 30
    # nz = 6
    # nx = 28
    # ny = 28
    # nz = 28
    nx = 28
    ny = 1
    nz = 28
    dx = 100.0
    dy = 100.0
    dz = 80.0

    # derivative = 0 # no derivative
    derivative = 1 # calculate derivative w.r.t. ksca

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
    solmu = np.cos(np.radians(30.0))
    # solmu = 0.5
    solphi = 0.5 * np.pi

    # viewmu = 1.0
    # viewphi = 0.0
    viewmu = np.cos(np.radians(15.0))
    viewphi = np.radians(90.0)

    # Ncpu = None
    Ncpu = 8

    fmc = FullMC(
            nx=nx, ny=ny, nz=nz,
            dx=dx, dy=dy, dz=dz,
            transfermode=transfermode,
            derivative=derivative,
            source=source,
            swlw=swlw,
            solmu=solmu, solphi=solphi,
            viewmu=viewmu, viewphi=viewphi,
            nphoton=nphoton,
            debug=debug,
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
    # kext[8:16, 8:16, 2:4] = 5e-1
    # kabs = np.zeros((nx, ny, nz), dtype=np.float64)
    # kabs[:, :, :] = 1e-6
    # kabs[8:16, 8:16, 2:4] = 1e-6
    # gparam = np.zeros((nx, ny, nz), dtype=np.float64)
    # gparam[:, :, :] = 0.0001
    # gparam[8:16, 8:16, 2:4] = 0.85
    # kext = np.zeros((nx, ny, nz), dtype=np.float64)
    # kext[:, :, :] = 1e-5
    # kext[nx//5:nx//5*4, 0:ny, 5:7] = 5e-1
    # kabs = np.zeros((nx, ny, nz), dtype=np.float64)
    # kabs[:, :, :] = 1e-6
    # kabs[nx//5:nx//5*4, 0:ny, 5:7] = 1e-6
    # gparam = np.zeros((nx, ny, nz), dtype=np.float64)
    # gparam[:, :, :] = 0.0001
    # gparam[nx//5:nx//5*4, 0:ny, 5:7] = 0.85
    # fmc.kext[:, :, :] = kext
    # fmc.kabs[:, :, :] = kabs
    # fmc.gparam[:, :, :] = gparam

    fmc.bplnk[:, :, :] = 0.0
    fmc.bgrnd[:, :] = 0.0
    fmc.galb[:, :] = 0.3

    fmc.read_atmtxt('dat/mod/les_mod_01.txt')
    # fmc.read_atmtxt('dat/mod/les_mod_02.txt')

    fmc.run_mc()
    radimg1 = fmc.read_result(kind='img')
    radtracer1 = fmc.read_result(kind='tracer1')

    # radimg2 = radimg1

    fmc2 = fmc
    fmc2.transfermode = 0 # ICA

    fmc2.run_mc()
    radimg2 = fmc2.read_result(kind='img')
    radtracer2 = fmc2.read_result(kind='tracer1')

    data_plot = np.pi * radimg1[:, :, 0]
    data_plot2 = np.pi * radimg2[:, :, 0]
    data_max = max(np.max(data_plot), np.max(data_plot2))
    diff = data_plot2 - data_plot
    diff[data_plot2 < np.max(fmc.galb) * 1.1] = np.nan

    xx, yy = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, ny * 100.0, ny))

    fig = plt.figure(figsize=(5, 3.2))
    # fig = plt.figure(figsize=(7, 3.2))
    # ax2 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 2, 1)
    m2 = ax2.imshow(data_plot2, cmap='Greys_r', vmin=0., vmax=data_max)
    # m2 = ax2.pcolormesh(xx*1e-3, yy*1e-3, data_plot2.T, cmap='Greys_r', vmin=0., vmax=data_max)
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (km)')
    ax2.set_ylabel('Y (km)')
    ax2.set_title('(a) IPA')
    # ax = fig.add_subplot(1, 3, 2)
    ax = fig.add_subplot(1, 2, 2)
    m = ax.imshow(data_plot, cmap='Greys_r', vmin=0., vmax=data_max)
    # m = ax.pcolormesh(xx*1e-3, yy*1e-3, data_plot.T, cmap='Greys_r', vmin=0., vmax=data_max)
    ax.set_aspect('equal')
    ax.set_xlabel('X (km)')
    # ax.set_ylabel('Y (km)')
    ax.set_yticklabels([])
    ax.set_title('(b) 3D')
    fig.colorbar(m, ax=[ax2, ax], orientation='horizontal', label='Reflectance', pad=0.22, shrink=0.45)
    # ax3 = fig.add_subplot(1, 3, 3)
    # ax3.pcolormesh(xx*1e-3, yy*1e-3, 0.75*np.isnan(diff).T.astype(float), cmap='Greys_r', vmin=0., vmax=1.)
    # m3 = ax3.pcolormesh(xx*1e-3, yy*1e-3, diff.T, cmap='seismic', vmin=-np.nanmax(np.abs(diff)), vmax=np.nanmax(np.abs(diff)))
    # ax3.set_aspect('equal')
    # ax3.set_xlabel('X (km)')
    # # ax3.set_ylabel('Y (km)')
    # ax3.set_yticklabels([])
    # ax3.set_title('(c) IPA-3D')
    # fig.colorbar(m3, ax=ax3, orientation='horizontal', label='Difference in Reflectance', pad=0.22)

    fig.savefig(f'{work_dir}/01_radtoa_toa_detector.png', dpi=300, bbox_inches='tight')

    mask = fmc.ksca[:, 0, :] > 1e-4
    x_edges = np.linspace(0.0, nx * dx * 1e-3, nx + 1)
    z_edges = np.linspace(0.0, nz * dz * 1e-3, nz + 1)

    cloud_segs = []
    for i in range(nx):
        for j in range(nz):
            if not mask[i, j]:
                continue
            x0, x1 = x_edges[i], x_edges[i + 1]
            z0, z1 = z_edges[j], z_edges[j + 1]
            if i == 0 or not mask[i - 1, j]:
                cloud_segs.append([(x0, z0), (x0, z1)])
            if i == nx - 1 or not mask[i + 1, j]:
                cloud_segs.append([(x1, z0), (x1, z1)])
            if j == 0 or not mask[i, j - 1]:
                cloud_segs.append([(x0, z0), (x1, z0)])
            if j == nz - 1 or not mask[i, j + 1]:
                cloud_segs.append([(x0, z1), (x1, z1)])

    data_max = max(
        np.max(radtracer1[:, :, :, 0, :, 0] / radimg1[np.newaxis, np.newaxis, :, :, 0]), 
        np.max(radtracer2[:, :, :, 0, :, 0] / radimg2[np.newaxis, np.newaxis, :, :, 0])
    )
    for ix in range(nx):
        for iy in range(ny):
            data_plot = radtracer1[ix, iy, :, 0, :, 0] / radimg1[ix, iy, 0]
            data_plot2 = radtracer2[ix, iy, :, 0, :, 0] / radimg2[ix, iy, 0]
            xx, zz = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, nz * 80.0, nz))

            fig = plt.figure(figsize=(5, 3.2))
            ax = fig.add_subplot(1, 2, 1)
            lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
            ax.add_collection(lc)
            m = ax.pcolormesh(xx*1e-3, zz*1e-3, data_plot2.T, cmap='viridis', vmin=0., vmax=data_max)
            ax.set_aspect('equal')
            ax.set_xlabel('X (km)')
            ax.set_ylabel('Z (km)')
            ax.set_title('(a) IPA')

            ax = fig.add_subplot(1, 2, 2)
            lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
            ax.add_collection(lc)
            m = ax.pcolormesh(xx*1e-3, zz*1e-3, data_plot.T, cmap='viridis', vmin=0., vmax=data_max)
            ax.set_aspect('equal')
            ax.set_xlabel('X (km)')
            ax.set_yticklabels([])
            ax.set_title('(b) 3D')
            fig.colorbar(m, ax=[fig.axes[0], fig.axes[1]], orientation='horizontal', label='Sensitivity', pad=0.22, shrink=0.45)
            fig.suptitle(f'Sensitivity (X={ix}, Y={iy})')

            fig.savefig(f'{work_dir}/01_radtoa_tracer1_x{ix:05d}_y{iy:05d}.png', dpi=300, bbox_inches='tight')

            plt.close(fig)
    
    plot_data = (fmc.ksca[:, 0, :] + fmc.kabs[:, 0, :, 0]) * dz

    fig = plt.figure(figsize=(5, 3.2))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect('equal')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Z (km)')
    pcm = ax.pcolormesh(xx*1e-3, zz*1e-3, plot_data.T, cmap='Greys_r', vmin=0., vmax=np.max(plot_data))
    ax.set_title('Optical Depth')
    fig.colorbar(pcm, ax=ax, orientation='horizontal', label='Optical Depth', pad=0.22, shrink=0.45)
    fig.savefig(f'{work_dir}/01_radtoa_optical_depth.png', dpi=300, bbox_inches='tight')