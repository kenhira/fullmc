import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

exe_mc = 'exe_fullmc'

out_dir = 'test_01_radtoa'

debug = 0
# debug = 1

nx = 15
ny = 15
nz = 5
dx = 100.0
dy = 100.0
dz = 80.0

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
# nphoton = 2e1
# nphoton = 3e2
nphoton = 5e3
# nphoton = 2e4

# solmu = 1.0
solmu = np.sqrt(2.0) / 2.0
# solmu = 0.5
solphi = 0.5 * np.pi

viewmu = 1.0
viewphi = 0.0

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

kext = np.zeros((nx, ny, nz), dtype=np.float64)
kext[:, :, :] = 1e-5
kext[4:8, 4:8, 2:4] = 5e-1
kabs = np.zeros((nx, ny, nz), dtype=np.float64)
kabs[:, :, :] = 1e-6
kabs[4:8, 4:8, 2:4] = 1e-6
gparam = np.zeros((nx, ny, nz), dtype=np.float64)
gparam[:, :, :] = 0.0001
gparam[4:8, 4:8, 2:4] = 0.85
bplnk = np.zeros((nx, ny, nz), dtype=np.float64)
bgrnd = np.zeros((nx, ny), dtype=np.float64)
galb = np.zeros((nx, ny), dtype=np.float64)
galb[:, :] = 0.0

# nx = 10
# ny = 1
# nz = 10
# dx = 100.0
# dy = 100.0
# dz = 80.0
# source = 0
# solmu = 0.5
# solphi = -0.5 * np.pi

# nphoton = 1
# debug = 1

# kext = np.zeros((nx, ny, nz), dtype=np.float64)
# kext[:, :, :] = 1e-10
# # kext[0:2, 0:ny, -4:-2] = 1e-3
# kabs = np.zeros((nx, ny, nz), dtype=np.float64)
# kabs[:, :, :] = 1e-9
# # kabs[0:2, 0:ny, -4:-2] = 1e-4
# gparam = np.zeros((nx, ny, nz), dtype=np.float64)
# gparam[:, :, :] = 0.0001
# # gparam[0:2, 0:ny, -4:-2] = 0.85

os.system(f'mkdir -p out/{out_dir}')

with open(f'out/{out_dir}/config.txt', 'w') as fh:
    fh.write("%d %d %d\n" % (nx, ny, nz))
    fh.write("%g %g %g\n" % (dx, dy, dz))
    fh.write("%d %d\n" % (source, swlw))
    fh.write("%d\n" % transfermode)
    fh.write("%g %g\n" % (solmu, solphi))
    fh.write("%g %g\n" % (viewmu, viewphi))
    fh.write("%d\n" % int(nphoton))
    fh.write("%d\n" % debug)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                fh.write("%15.6e %15.6e %15.6e %15.6e\n" % (kext[ix,iy,iz], kabs[ix,iy,iz], gparam[ix,iy,iz], bplnk[ix,iy,iz]))
            fh.write("%15.6e %15.6e\n" % (galb[ix,iy], bgrnd[ix,iy]))
    

# os.system(f'./{exe_mc}')
# os.system(f'./{exe_mc} > out/{out_dir}/log.txt')
os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log.txt')

with open(f'out/{out_dir}/outradirr.txt', 'r') as fh:
    header = fh.readline()  # skip header
    dims_line = fh.readline()
    # nx, ny, nz, ncase, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
    nx, ny, nz, nd, na, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
    
    radirr = np.zeros((nx, ny, nz, nd, na, ncomp), dtype=np.float64)
    
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for idi in range(nd):
                    for ia in range(na):
                        line = fh.readline()
                        parts = line.strip().split()
                        radirr[ix, iy, iz, idi, ia, 0] = float(parts[0])
                        radirr[ix, iy, iz, idi, ia, 1] = float(parts[1])

with open(f'out/{out_dir}/outradconv.txt', 'r') as fh:
    header = fh.readline()  # skip header
    dims_line = fh.readline()
    # nx, ny, nz, ncase, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
    nx, ny, nz, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
    
    radconv = np.zeros((nx, ny, nz, ncomp), dtype=np.float64)
    
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                # for icase in range(ncase):
                line = fh.readline()
                parts = line.strip().split()
                radconv[ix, iy, iz, 0] = float(parts[0])
                radconv[ix, iy, iz, 1] = float(parts[1])

with open(f'out/{out_dir}/outradimg.txt', 'r') as fh:
    header = fh.readline()  # skip header
    dims_line = fh.readline()
    # nx, ny, ncase, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
    nx, ny, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
    
    radimg = np.zeros((nx, ny, ncomp), dtype=np.float64)
    
    for ix in range(nx):
        for iy in range(ny):
            # for icase in range(ncase):
            line = fh.readline()
            parts = line.strip().split()
            radimg[ix, iy, 0] = float(parts[0])
            radimg[ix, iy, 1] = float(parts[1])

if swlw == 0 and source == 2:
    indy = 0
    # data_plot = radirr[:, indy, :, 0, 0, 0]
    # data_plot2 = radirr[:, indy, :, 1, 0, 0]
    # data_plot3 = radirr[:, indy, :, 1, 1, 0]
    data_plot4 = radconv[:, indy, :, 0]
    norm_symlog = mcolors.SymLogNorm(
        linthresh=0.001*np.max(np.abs(data_plot4)),
        vmin=-np.max(np.abs(data_plot4)), 
        vmax=np.max(np.abs(data_plot4)),
        # clip=True
    )

    xx, zz = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, nz * 80.0, nz))

    fig = plt.figure(figsize=(8, 4))
    # fig = plt.figure(figsize=(16, 8))
    # ax = fig.add_subplot(2, 2, 1)
    # m1 = ax.pcolormesh(xx, zz, data_plot.T, cmap='cividis', vmin=0.)
    # ax.set_aspect('equal')
    # fig.colorbar(m1, ax=ax)
    # ax.set_xlabel('X (m)')
    # ax.set_ylabel('Z (m)')
    # ax.set_title(f'Y={indy * dy:.1f} m - Direct Irradiance')
    # ax2 = fig.add_subplot(2, 2, 2)
    # m2 = ax2.pcolormesh(xx, zz, (data_plot + data_plot2).T, cmap='cividis', vmin=0.)
    # ax2.set_aspect('equal')
    # fig.colorbar(m2, ax=ax2)
    # ax2.set_xlabel('X (m)')
    # ax2.set_ylabel('Z (m)')
    # ax2.set_title(f'Y={indy * dy:.1f} m - Global Downward Irradiance')
    # ax3 = fig.add_subplot(2, 2, 3)
    # m3 = ax3.pcolormesh(xx, zz, data_plot3.T, cmap='cividis', vmin=0.)
    # ax3.set_aspect('equal')
    # fig.colorbar(m3, ax=ax3)
    # ax3.set_xlabel('X (m)')
    # ax3.set_ylabel('Z (m)')
    # ax3.set_title(f'Y={indy * dy:.1f} m - Upward Irradiance')
    # ax4 = fig.add_subplot(2, 2, 4)
    ax4 = fig.add_subplot(1, 1, 1)
    m4 = ax4.pcolormesh(xx, zz, data_plot4.T, cmap='coolwarm', norm=norm_symlog)
    # m4 = ax4.pcolormesh(xx, zz, data_plot4.T, cmap='coolwarm', vmin=-np.abs(data_plot4).max(), vmax=np.abs(data_plot4).max())
    ax4.set_aspect('equal')
    fig.colorbar(m4, ax=ax4)
    ax4.set_xlabel('X (m)')
    ax4.set_ylabel('Z (m)')
    ax4.set_title(f'Y={indy * dy:.1f} m - Radiative Convergence')
    fig.tight_layout()

    plt.show()

if swlw == 1 and source in [0, 1]:
    indy = 0
    data_plot = radirr[:, indy, :, 0, 0, 0]
    data_plot2 = radirr[:, indy, :, 1, 0, 0]
    data_plot3 = radirr[:, indy, :, 1, 1, 0]
    data_plot4 = radconv[:, indy, :, 0]

    norm_log = mcolors.LogNorm(
        # vmin=np.max(np.abs(data_plot4)) * 1e-4,
        vmax=np.max(np.abs(data_plot4)),
        # clip=True
    )

    xx, zz = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, nz * 80.0, nz))

    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(2, 2, 1)
    m1 = ax.pcolormesh(xx, zz, data_plot.T, cmap='cividis', vmin=0.)
    ax.set_aspect('equal')
    fig.colorbar(m1, ax=ax)
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(f'Y={indy * dy:.1f} m - Direct Irradiance')
    ax2 = fig.add_subplot(2, 2, 2)
    m2 = ax2.pcolormesh(xx, zz, (data_plot + data_plot2).T, cmap='cividis', vmin=0.)
    ax2.set_aspect('equal')
    fig.colorbar(m2, ax=ax2)
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Z (m)')
    ax2.set_title(f'Y={indy * dy:.1f} m - Global Downward Irradiance')
    ax3 = fig.add_subplot(2, 2, 3)
    m3 = ax3.pcolormesh(xx, zz, data_plot3.T, cmap='cividis', vmin=0.)
    ax3.set_aspect('equal')
    fig.colorbar(m3, ax=ax3)
    ax3.set_xlabel('X (m)')
    ax3.set_ylabel('Z (m)')
    ax3.set_title(f'Y={indy * dy:.1f} m - Upward Irradiance')
    ax4 = fig.add_subplot(2, 2, 4)
    m4 = ax4.pcolormesh(xx, zz, data_plot4.T, cmap='inferno', norm=norm_log)
    ax4.set_aspect('equal')
    fig.colorbar(m4, ax=ax4)
    ax4.set_xlabel('X (m)')
    ax4.set_ylabel('Z (m)')
    ax4.set_title(f'Y={indy * dy:.1f} m - Radiative Convergence')
    fig.tight_layout()

    indz = 0
    data_plot = radirr[:, :, indz, 0, 0, 0] + radirr[:, :, indz, 1, 0, 0]
    hist, bin_edges = np.histogram(data_plot, bins=50, range=(0.0, data_plot.max()*1.1), density=False)

    fig2 = plt.figure(figsize=(7, 6))
    ax = fig2.add_subplot(111)
    ax.plot(0.5 * (bin_edges[1:] + bin_edges[:-1]), hist, drawstyle='steps-mid')
    ax.set_xlim(0.0, data_plot.max()*1.1)
    ax.set_xlabel('Irradiance at Z={:.1f} m (W/m²)'.format(indz * dz))
    ax.set_ylabel('Counts')
    ax.set_title('Histogram of Irradiance at Z={:.1f} m'.format(indz * dz))
    fig2.tight_layout()

    plt.show()

if swlw == 1 and source == 4:
    data_plot = np.pi * radimg[:, :, 0]
    data_plot2 = np.pi ** 2 * radimg[:, :, 1]
    xx, yy = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, ny * 100.0, ny))

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 2, 1)
    m = ax.pcolormesh(xx, yy, data_plot.T, cmap='viridis', vmin=0.)
    fig.colorbar(m, ax=ax)
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title('TOA Detector Irradiance')
    ax2 = fig.add_subplot(1, 2, 2)
    m2 = ax2.pcolormesh(xx, yy, data_plot2.T, cmap='viridis', vmin=0.)
    fig.colorbar(m2, ax=ax2)
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    ax2.set_title('TOA Detector Irradiance Variance')
    
    plt.show()