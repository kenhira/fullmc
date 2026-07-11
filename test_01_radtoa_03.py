import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
import os
import copy
import pickle

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
    ny = 28
    nz = 28
    dx = 100.0
    dy = 100.0
    dz = 80.0

    derivative = 0 # no derivative
    # derivative = 1 # calculate derivative w.r.t. ksca

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
    # nphoton = 1e3
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


    fmc.bplnk[:, :, :] = 0.0
    fmc.bgrnd[:, :] = 0.0
    fmc.galb[:, :] = 0.0
    fmc.read_atmtxt('dat/mod/les_mod_02.txt')
    fmc.derivative = 0

    fmc1 = copy.deepcopy(fmc)
    fmc1.viewmu = np.cos(np.radians(30.0))
    fmc1.viewphi = np.radians(90.0)
    fmc1.run_mc()
    radimg1 = fmc1.read_result(kind='img')
    # radtracer1 = fmc1.read_result(kind='tracer1')

    fmc2 = copy.deepcopy(fmc)
    fmc2.viewmu = np.cos(np.radians(0.0))
    fmc2.viewphi = np.radians(0.0)
    fmc2.run_mc()
    radimg2 = fmc2.read_result(kind='img')
    # radtracer2 = fmc2.read_result(kind='tracer1')

    fmc3 = copy.deepcopy(fmc)
    fmc3.viewmu = np.cos(np.radians(30.0))
    fmc3.viewphi = np.radians(-90.0)
    fmc3.run_mc()
    radimg3 = fmc3.read_result(kind='img')
    # radtracer3 = fmc3.read_result(kind='tracer1')

    fmc_i = copy.deepcopy(fmc)
    fmc_i.derivative = 1
    mask = fmc.ksca[:, :, :] > 1e-4
    fmc_i.ksca[mask] = fmc.ksca[mask] * 0.5
    fmc_i.ksca[~mask] = fmc.ksca[~mask]

    fmc_i1 = copy.deepcopy(fmc_i)
    fmc_i1.viewmu = np.cos(np.radians(30.0))
    fmc_i1.viewphi = np.radians(90.0)
    fmc_i1.run_mc()
    radimg_i1 = fmc_i1.read_result(kind='img')
    radtracer_i1 = fmc_i1.read_result(kind='tracer1')

    fmc_i2 = copy.deepcopy(fmc_i)
    fmc_i2.viewmu = np.cos(np.radians(0.0))
    fmc_i2.viewphi = np.radians(0.0)
    fmc_i2.run_mc()
    radimg_i2 = fmc_i2.read_result(kind='img')
    radtracer_i2 = fmc_i2.read_result(kind='tracer1')

    fmc_i3 = copy.deepcopy(fmc_i)
    fmc_i3.viewmu = np.cos(np.radians(30.0))
    fmc_i3.viewphi = np.radians(-90.0)
    fmc_i3.run_mc()
    radimg_i3 = fmc_i3.read_result(kind='img')
    radtracer_i3 = fmc_i3.read_result(kind='tracer1')

    update = np.full((nx, ny, 3, nx, ny, nz), np.nan)
    imgdiff1 = radimg_i1[:, :, 0] - radimg1[:, :, 0]
    imgdiff2 = radimg_i2[:, :, 0] - radimg2[:, :, 0]
    imgdiff3 = radimg_i3[:, :, 0] - radimg3[:, :, 0]
    gain1 = np.where(radtracer_i1[:, :, :, :, :, 0] < 0.005, np.nan, -1.0 / radtracer_i1[:, :, :, :, :, 0])
    gain2 = np.where(radtracer_i2[:, :, :, :, :, 0] < 0.005, np.nan, -1.0 / radtracer_i2[:, :, :, :, :, 0])
    gain3 = np.where(radtracer_i3[:, :, :, :, :, 0] < 0.005, np.nan, -1.0 / radtracer_i3[:, :, :, :, :, 0])
    update[:, :, 0, :, :, :] = imgdiff1[:, :, np.newaxis, np.newaxis, np.newaxis] * gain1 * fmc_i.ksca[np.newaxis, np.newaxis, np.newaxis, :, :, :]
    update[:, :, 1, :, :, :] = imgdiff2[:, :, np.newaxis, np.newaxis, np.newaxis] * gain2 * fmc_i.ksca[np.newaxis, np.newaxis, np.newaxis, :, :, :]
    update[:, :, 2, :, :, :] = imgdiff3[:, :, np.newaxis, np.newaxis, np.newaxis] * gain3 * fmc_i.ksca[np.newaxis, np.newaxis, np.newaxis, :, :, :]

    data_max = max(
        np.max(radimg1[:, :, 0]), 
        np.max(radimg2[:, :, 0]),
        np.max(radimg3[:, :, 0])
    )

    fig = plt.figure(figsize=(5, 5))
    ax1 = fig.add_subplot(1, 3, 1)
    m1 = ax1.imshow(radimg1[:, :, 0], cmap='Greys_r', vmin=0., vmax=data_max)
    ax1.set_aspect('equal')
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_title('(a) Angle 1')
    ax2 = fig.add_subplot(1, 3, 2)
    m2 = ax2.imshow(radimg2[:, :, 0], cmap='Greys_r', vmin=0., vmax=data_max)
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (km)')
    # ax2.set_ylabel('Y (km)')
    ax2.set_yticklabels([])
    ax2.set_title('(b) Angle 2')
    ax3 = fig.add_subplot(1, 3, 3)
    m3 = ax3.imshow(radimg3[:, :, 0], cmap='Greys_r', vmin=0., vmax=data_max)
    ax3.set_aspect('equal')
    ax3.set_xlabel('X (km)')
    # ax3.set_ylabel('Y (km)')
    ax3.set_yticklabels([])
    ax3.set_title('(c) Angle 3')
    fig.colorbar(m1, ax=[ax1,ax2, ax3], orientation='horizontal', label='Radiance', pad=0.22, shrink=0.45)
    fig.suptitle('True radiance')
    fig.savefig(f'{work_dir}/01_radtoa_toa_detector_true.png', dpi=300, bbox_inches='tight')

    fig2 = plt.figure(figsize=(5, 5))
    ax1 = fig2.add_subplot(1, 3, 1)
    m1 = ax1.imshow(radimg_i1[:, :, 0], extent=[0, ny * dy * 1e-3, 0, nx * dx * 1e-3,], cmap='Greys_r', vmin=0., vmax=data_max)
    ax1.set_aspect('equal')
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_title('(a) Angle 1')
    ax2 = fig2.add_subplot(1, 3, 2)
    m2 = ax2.imshow(radimg_i2[:, :, 0], extent=[0, ny * dy * 1e-3, 0, nx * dx * 1e-3], cmap='Greys_r', vmin=0., vmax=data_max)
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (km)')
    # ax2.set_ylabel('Y (km)')
    ax2.set_yticklabels([])
    ax2.set_title('(b) Angle 2')
    ax3 = fig2.add_subplot(1, 3, 3)
    m3 = ax3.imshow(radimg_i3[:, :, 0], extent=[0, ny * dy * 1e-3, 0, nx * dx * 1e-3], cmap='Greys_r', vmin=0., vmax=data_max)
    ax3.set_aspect('equal')
    ax3.set_xlabel('X (km)')
    # ax3.set_ylabel('Y (km)')
    ax3.set_yticklabels([])
    ax3.set_title('(c) Angle 3')
    fig2.colorbar(m1, ax=[ax1,ax2, ax3], orientation='horizontal', label='Radiance', pad=0.22, shrink=0.45)
    fig2.suptitle('Radiance simulated with first-guess clouds')
    fig2.savefig(f'{work_dir}/01_radtoa_toa_detector_guess.png', dpi=300, bbox_inches='tight')

    fig3 = plt.figure(figsize=(5, 5))
    ax1 = fig3.add_subplot(1, 3, 1)
    m1 = ax1.imshow(imgdiff1[:, :], extent=[0, ny * dy * 1e-3, 0, nx * dx * 1e-3], cmap='RdBu_r', vmin=-0.5*data_max, vmax=0.5*data_max)
    ax1.set_aspect('equal')
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_title('(a) Angle 1')
    ax2 = fig3.add_subplot(1, 3, 2)
    m2 = ax2.imshow(imgdiff2[:, :], extent=[0, ny * dy * 1e-3, 0, nx * dx * 1e-3], cmap='RdBu_r', vmin=-0.5*data_max, vmax=0.5*data_max)
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (km)')
    # ax2.set_ylabel('Y (km)')
    ax2.set_yticklabels([])
    ax2.set_title('(b) Angle 2')
    ax3 = fig3.add_subplot(1, 3, 3)
    m3 = ax3.imshow(imgdiff3[:, :], extent=[0, ny * dy * 1e-3, 0, nx * dx * 1e-3], cmap='RdBu_r', vmin=-0.5*data_max, vmax=0.5*data_max)
    ax3.set_aspect('equal')
    ax3.set_xlabel('X (km)')
    # ax3.set_ylabel('Y (km)')
    ax3.set_yticklabels([])
    ax3.set_title('(c) Angle 3')
    fig3.colorbar(m1, ax=[ax1,ax2, ax3], orientation='horizontal', label='Radiance Difference', pad=0.2, shrink=0.4)
    fig3.suptitle('Radiance difference (first-guess - true)', fontsize=10)
    fig3.savefig(f'{work_dir}/01_radtoa_toa_detector_diff.png', dpi=300, bbox_inches='tight')

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
        np.max(radtracer_i1[:, :, :, :, :, 0]), 
        np.max(radtracer_i2[:, :, :, :, :, 0]),
        np.max(radtracer_i3[:, :, :, :, :, 0])
    )
    for ix in range(nx):
        for iy in range(ny):
            # data_plot = radtracer_i1[ix, iy, :, 0, :, 0] / radimg_i1[ix, iy, 0, np.newaxis, np.newaxis] * fmc_i.ksca[:, 0, :]
            # data_plot2 = radtracer_i2[ix, iy, :, 0, :, 0] / radimg_i2[ix, iy, 0, np.newaxis, np.newaxis] * fmc_i.ksca[:, 0, :]
            # data_plot3 = radtracer_i3[ix, iy, :, 0, :, 0] / radimg_i3[ix, iy, 0, np.newaxis, np.newaxis] * fmc_i.ksca[:, 0, :]
            data_plot = radtracer_i1[ix, iy, :, iy, :, 0]
            data_plot2 = radtracer_i2[ix, iy, :, iy, :, 0]
            data_plot3 = radtracer_i3[ix, iy, :, iy, :, 0]
            xx, zz = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, nz * 80.0, nz))

            fig = plt.figure(figsize=(9, 3.2))
            ax1 = fig.add_subplot(1, 3, 1)
            lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
            ax1.add_collection(lc)
            m = ax1.pcolormesh(xx*1e-3, zz*1e-3, data_plot.T, cmap='viridis', vmin=0., vmax=data_max)
            ax1.set_aspect('equal')
            ax1.set_xlabel('X (km)')
            ax1.set_ylabel('Z (km)')
            ax1.set_title('(a) Angle 1')

            ax2 = fig.add_subplot(1, 3, 2)
            lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
            ax2.add_collection(lc)
            m = ax2.pcolormesh(xx*1e-3, zz*1e-3, data_plot2.T, cmap='viridis', vmin=0., vmax=data_max)
            ax2.set_aspect('equal')
            ax2.set_xlabel('X (km)')
            ax2.set_yticklabels([])
            ax2.set_title('(b) Angle 2')

            ax3 = fig.add_subplot(1, 3, 3)
            lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
            ax3.add_collection(lc)
            m = ax3.pcolormesh(xx*1e-3, zz*1e-3, data_plot3.T, cmap='viridis', vmin=0., vmax=data_max)
            ax3.set_aspect('equal')
            ax3.set_xlabel('X (km)')
            ax3.set_yticklabels([])
            ax3.set_title('(c) Angle 3')

            fig.colorbar(m, ax=[ax1, ax2, ax3], orientation='horizontal', label=r'$k \partial_{k} I$', pad=0.2, shrink=0.4)
            fig.suptitle(f'Sensitivity (X={ix}, Y={iy})', fontsize=10)

            fig.savefig(f'{work_dir}/01_radtoa_tracer1_x{ix:05d}_y{iy:05d}.png', dpi=300, bbox_inches='tight')

            plt.close(fig)
    
    for ix in range(nx):
        for iy in range(ny):
            data_plot = update[ix, iy, 0, :, iy, :]
            data_plot2 = update[ix, iy, 1, :, iy, :]
            data_plot3 = update[ix, iy, 2, :, iy, :]
            xx, zz = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, nz * 80.0, nz))

            fig = plt.figure(figsize=(9, 3.2))
            ax1 = fig.add_subplot(1, 3, 1)
            lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
            ax1.add_collection(lc)
            m = ax1.pcolormesh(xx*1e-3, zz*1e-3, data_plot.T, cmap='coolwarm', vmin=-1., vmax=1.)
            ax1.set_aspect('equal')
            ax1.set_xlabel('X (km)')
            ax1.set_ylabel('Z (km)')
            ax1.set_title('(a) Angle 1')

            ax2 = fig.add_subplot(1, 3, 2)
            lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
            ax2.add_collection(lc)
            m = ax2.pcolormesh(xx*1e-3, zz*1e-3, data_plot2.T, cmap='coolwarm', vmin=-1., vmax=1.)
            ax2.set_aspect('equal')
            ax2.set_xlabel('X (km)')
            ax2.set_yticklabels([])
            ax2.set_title('(b) Angle 2')

            ax3 = fig.add_subplot(1, 3, 3)
            lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
            ax3.add_collection(lc)
            m = ax3.pcolormesh(xx*1e-3, zz*1e-3, data_plot3.T, cmap='coolwarm', vmin=-1., vmax=1.)
            ax3.set_aspect('equal')
            ax3.set_xlabel('X (km)')
            ax3.set_yticklabels([])
            ax3.set_title('(c) Angle 3')

            fig.colorbar(m, ax=[ax1, ax2, ax3], orientation='horizontal', label=r'$\Delta I / \partial_{k} I$ (1/m)', pad=0.22, shrink=0.45)
            fig.suptitle(f'Update (X={ix}, Y={iy})')

            fig.savefig(f'{work_dir}/01_radtoa_update_x{ix:05d}_y{iy:05d}.png', dpi=300, bbox_inches='tight')

            plt.close(fig)
    
    update_avg = np.nanmean(update, axis=(0, 1, 2))

    for iy in range(ny):

        data_plot1 = fmc.ksca[:, iy, :]
        data_plot2 = fmc_i.ksca[:, iy, :]
        data_plot3 = data_plot1 - data_plot2
        data_plot4 = update_avg[:, iy, :]
        maxval = max(np.nanmax(data_plot1), np.nanmax(data_plot2))
        maxabs = max(np.nanmax(np.abs(data_plot3)), np.nanmax(np.abs(data_plot4)))

        fig = plt.figure(figsize=(12, 4.8))
        ax1 = fig.add_subplot(1, 4, 1)
        lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
        ax1.add_collection(lc)
        pcm = ax1.pcolormesh(xx*1e-3, zz*1e-3, data_plot1.T, cmap='cividis', vmin=0., vmax=maxval)
        ax1.set_aspect('equal')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('Z (km)')
        ax1.set_title('(a) True scattering coefficient', fontsize=10)
        ax2 = fig.add_subplot(1, 4, 2)
        lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
        ax2.add_collection(lc)
        pcm = ax2.pcolormesh(xx*1e-3, zz*1e-3, data_plot2.T, cmap='cividis', vmin=0., vmax=maxval)
        ax2.set_aspect('equal')
        ax2.set_xlabel('X (km)')
        ax2.set_yticklabels([])
        ax2.set_title('(b) First-guess', fontsize=10)
        ax3 = fig.add_subplot(1, 4, 3)
        lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
        ax3.add_collection(lc)
        pcm2 = ax3.pcolormesh(xx*1e-3, zz*1e-3, data_plot3.T, cmap='coolwarm', vmin=-maxabs, vmax=maxabs)
        ax3.set_aspect('equal')
        ax3.set_xlabel('X (km)')
        ax3.set_yticklabels([])
        # ax3.set_ylabel('Z (km)')
        ax3.set_title('(c) Truth - First-guess', fontsize=10)
        ax4 = fig.add_subplot(1, 4, 4)
        lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
        ax4.add_collection(lc)
        pcm2 = ax4.pcolormesh(xx*1e-3, zz*1e-3, data_plot4.T, cmap='coolwarm', vmin=-maxabs, vmax=maxabs)
        ax4.set_aspect('equal')
        ax4.set_xlabel('X (km)')
        ax4.set_yticklabels([])
        ax4.set_title('(d) Mean Update', fontsize=10)
        ax4.set_aspect('equal')

        fig.colorbar(pcm, ax=[ax1, ax2], orientation='horizontal', label='(1/m)', pad=0.15, shrink=0.4)
        fig.colorbar(pcm2, ax=[ax3, ax4], orientation='horizontal', label='(1/m)', pad=0.15, shrink=0.4)

        fig.savefig(f'{work_dir}/01_radtoa_update_avg_{iy:05d}.png', dpi=300, bbox_inches='tight')


    # fig = plt.figure(figsize=(5, 3.2))
    # ax = fig.add_subplot(1, 1, 1)
    # lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
    # ax.add_collection(lc)
    # ax.set_aspect('equal')
    # ax.set_xlabel('X (km)')
    # ax.set_ylabel('Z (km)')
    # pcm = ax.pcolormesh(xx*1e-3, zz*1e-3, data_plot.T, cmap='coolwarm', vmin=-maxabs, vmax=maxabs)
    # ax.set_title('Average Update')
    # fig.colorbar(pcm, ax=ax, orientation='horizontal', label='Update', pad=0.22, shrink=0.45)
    # fig.savefig(f'{work_dir}/01_radtoa_update_avg.png', dpi=300, bbox_inches='tight')
    # plt.close(fig)

    # plot_data = fmc.ksca[:, 0, :]

    # fig = plt.figure(figsize=(5, 3.2))
    # ax = fig.add_subplot(1, 1, 1)
    # lc = LineCollection(cloud_segs, colors='k', linewidths=1.0)
    # ax.add_collection(lc)
    # ax.set_aspect('equal')
    # ax.set_xlabel('X (km)')
    # ax.set_ylabel('Z (km)')
    # pcm = ax.pcolormesh(xx*1e-3, zz*1e-3, plot_data.T, cmap='cividis', vmin=0., vmax=np.max(plot_data))
    # ax.set_title('Scattering Coefficient')
    # fig.colorbar(pcm, ax=ax, orientation='horizontal', label='(1/m)', pad=0.22, shrink=0.45)
    # fig.savefig(f'{work_dir}/01_radtoa_optical_depth.png', dpi=300, bbox_inches='tight')