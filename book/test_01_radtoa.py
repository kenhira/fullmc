import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
matplotlib.rc('font', family='Hiragino Sans')


# matplotlib.use('Agg')
# # prefer Helvetica, fall back to other common sans-serifs
# matplotlib.rcParams['font.family'] = 'sans-serif'
# matplotlib.rcParams['font.sans-serif'] = ['Helvetica', 'DejaVu Sans', 'Arial']
# # embed fonts in vector outputs when possible
# matplotlib.rcParams['pdf.fonttype'] = 42
# matplotlib.rcParams['ps.fonttype'] = 42


exe_mc = 'exe_fullmc'

out_dir = 'test_01_radtoa'

debug = 0
# debug = 1

# nx = 30
# ny = 30
# nz = 6
nx = 28
ny = 28
nz = 28
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
nphoton = 1e3
# nphoton = 5e3
# nphoton = 2e4
# nphoton = 5e4

# solmu = 1.0
# solmu = np.sqrt(2.0) / 2.0
solmu = np.cos(np.radians(30.0))
# solmu = 0.5
solphi = 0.5 * np.pi

viewmu = 1.0
viewphi = 0.0

seedval = 1235

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
bplnk = np.zeros((nx, ny, nz), dtype=np.float64)
bgrnd = np.zeros((nx, ny), dtype=np.float64)
galb = np.zeros((nx, ny), dtype=np.float64)
galb[:, :] = 0.1

with open(f'dat/mod/les_mod_02.txt', 'r') as fh:
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

os.system(f'mkdir -p out/{out_dir}')

with open(f'out/{out_dir}/config.txt', 'w') as fh:
    seedval = seedval + 11
    fh.write("%d %d %d\n" % (nx, ny, nz))
    fh.write("%g %g %g\n" % (dx, dy, dz))
    fh.write("%d %d\n" % (source, swlw))
    fh.write("%d\n" % transfermode)
    fh.write("%g %g\n" % (solmu, solphi))
    fh.write("%g %g\n" % (viewmu, viewphi))
    fh.write("%d\n" % int(nphoton))
    fh.write("%d\n" % seedval)
    fh.write("%d\n" % debug)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                fh.write("%15.6e %15.6e %15.6e %15.6e\n" % (kext[ix,iy,iz], kabs[ix,iy,iz], gparam[ix,iy,iz], bplnk[ix,iy,iz]))
            fh.write("%15.6e %15.6e\n" % (galb[ix,iy], bgrnd[ix,iy]))
    

# os.system(f'./{exe_mc}')
# os.system(f'./{exe_mc} > out/{out_dir}/log.txt')
os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log.txt')

with open(f'out/{out_dir}/outradimg.txt', 'r') as fh:
    header = fh.readline()  # skip header
    dims_line = fh.readline()
    nx, ny, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
    
    radimg = np.zeros((nx, ny, ncomp), dtype=np.float64)
    
    for ix in range(nx):
        for iy in range(ny):
            # for icase in range(ncase):
            line = fh.readline()
            parts = line.strip().split()
            radimg[ix, iy, 0] = float(parts[0])
            radimg[ix, iy, 1] = float(parts[1])

radimg1 = radimg.copy()

transfermode = 0 # ICA

with open(f'out/{out_dir}/config.txt', 'w') as fh:
    seedval = seedval + 11
    fh.write("%d %d %d\n" % (nx, ny, nz))
    fh.write("%g %g %g\n" % (dx, dy, dz))
    fh.write("%d %d\n" % (source, swlw))
    fh.write("%d\n" % transfermode)
    fh.write("%g %g\n" % (solmu, solphi))
    fh.write("%g %g\n" % (viewmu, viewphi))
    fh.write("%d\n" % int(nphoton))
    fh.write("%d\n" % seedval)
    fh.write("%d\n" % debug)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                fh.write("%15.6e %15.6e %15.6e %15.6e\n" % (kext[ix,iy,iz], kabs[ix,iy,iz], gparam[ix,iy,iz], bplnk[ix,iy,iz]))
            fh.write("%15.6e %15.6e\n" % (galb[ix,iy], bgrnd[ix,iy]))
    

# os.system(f'./{exe_mc}')
# os.system(f'./{exe_mc} > out/{out_dir}/log.txt')
os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log.txt')

with open(f'out/{out_dir}/outradimg.txt', 'r') as fh:
    header = fh.readline()  # skip header
    dims_line = fh.readline()
    nx, ny, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
    
    radimg = np.zeros((nx, ny, ncomp), dtype=np.float64)
    
    for ix in range(nx):
        for iy in range(ny):
            # for icase in range(ncase):
            line = fh.readline()
            parts = line.strip().split()
            radimg[ix, iy, 0] = float(parts[0])
            radimg[ix, iy, 1] = float(parts[1])

radimg2 = radimg.copy()

data_plot = np.pi * radimg1[:, :, 0]
data_plot2 = np.pi * radimg2[:, :, 0]
data_max = max(np.max(data_plot), np.max(data_plot2))
diff = data_plot2 - data_plot
diff[data_plot2 < np.max(galb) * 1.1] = np.nan

xx, yy = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, ny * 100.0, ny))

fig = plt.figure(figsize=(5, 3.2))
# fig = plt.figure(figsize=(7, 3.2))
# ax2 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 2, 1)
m2 = ax2.pcolormesh(xx*1e-3, yy*1e-3, data_plot2.T, cmap='Greys_r', vmin=0., vmax=data_max)
ax2.set_aspect('equal')
ax2.set_xlabel('X (km)')
ax2.set_ylabel('Y (km)')
ax2.set_title('(a) IPA')
# ax = fig.add_subplot(1, 3, 2)
ax = fig.add_subplot(1, 2, 2)
m = ax.pcolormesh(xx*1e-3, yy*1e-3, data_plot.T, cmap='Greys_r', vmin=0., vmax=data_max)
ax.set_aspect('equal')
ax.set_xlabel('X (km)')
# ax.set_ylabel('Y (km)')
ax.set_yticklabels([])
ax.set_title('(b) 3D')
fig.colorbar(m, ax=[ax2, ax], orientation='horizontal', label='反射率', pad=0.22, shrink=0.45)
# ax3 = fig.add_subplot(1, 3, 3)
# ax3.pcolormesh(xx*1e-3, yy*1e-3, 0.75*np.isnan(diff).T.astype(float), cmap='Greys_r', vmin=0., vmax=1.)
# m3 = ax3.pcolormesh(xx*1e-3, yy*1e-3, diff.T, cmap='seismic', vmin=-np.nanmax(np.abs(diff)), vmax=np.nanmax(np.abs(diff)))
# ax3.set_aspect('equal')
# ax3.set_xlabel('X (km)')
# # ax3.set_ylabel('Y (km)')
# ax3.set_yticklabels([])
# ax3.set_title('(c) IPA-3D')
# fig.colorbar(m3, ax=ax3, orientation='horizontal', label='反射率の差', pad=0.22)

fig.savefig(f'out/{out_dir}/01_radtoa_toa_detector.png', dpi=300, bbox_inches='tight')