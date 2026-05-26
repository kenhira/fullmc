import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from matplotlib.collections import LineCollection
import glob
import matplotlib
matplotlib.rc('font', family='Hiragino Sans')

exe_mc = 'exe_fullmc'

out_dir = 'test_03_fluxxz'

debug = 0
# debug = 1

# nx = 15
# ny = 15
# nz = 5
# nx = 28
# ny = 1
# nz = 28
nx = 25
ny = 1
nz = 16
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
# nphoton = 3e2
# nphoton = 5e3
# nphoton = 2e4
nphoton = 4e4

# solmu = 1.0
# solmu = np.sqrt(2.0) / 2.0
# solmu = np.cos(np.radians(0.0))
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
# kext[4:8, 4:8, 2:4] = 5e-1
# kabs = np.zeros((nx, ny, nz), dtype=np.float64)
# kabs[:, :, :] = 1e-6
# kabs[4:8, 4:8, 2:4] = 1e-6
# gparam = np.zeros((nx, ny, nz), dtype=np.float64)
# gparam[:, :, :] = 0.0001
# gparam[4:8, 4:8, 2:4] = 0.85
bplnk = np.zeros((nx, ny, nz), dtype=np.float64)
bgrnd = np.zeros((nx, ny), dtype=np.float64)
galb = np.zeros((nx, ny), dtype=np.float64)
galb[:, :] = 0.0

xx, zz = np.meshgrid(np.linspace(0.0, nx * dx, nx) + dx/2., np.linspace(0.0, nz * dz, nz) + dz/2.)


# file_les = 'dat/mod/les_mod_01.txt'
# files_les = glob.glob('dat/mod/lescrop_01_24-01_*.txt')
# files_les = glob.glob('dat/mod/lescrop_01_24*.txt')
files_les = sorted(glob.glob('dat/mod/lescrop_01_*.txt'))

data_plot = np.zeros((nx, nz), dtype=np.float64)
data_plot2 = np.zeros((nx, nz), dtype=np.float64)
data_plotb = np.zeros((nx), dtype=np.float64)
data_plotb2 = np.zeros((nx), dtype=np.float64)

for file_les in files_les:
    bnd = int(os.path.basename(file_les).split('_')[2].split('-')[0])
    if bnd < 15:
        continue
    print(file_les)
    with open(file_les, 'r') as fh:
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
    absscale = np.exp(-np.mean(kabs[:, :, -1])*8.e3)
    print('Absorption scale factor: {:.6e}'.format(absscale))
    
    wgt = float(os.path.basename(file_les).split('_')[-2])
    fsol = float(os.path.basename(file_les).split('_')[-1].split('.txt')[0])

    os.system(f'mkdir -p out/{out_dir}')

    transfermode = 0 # ICA

    with open(f'out/{out_dir}/config.txt', 'w') as fh:
        seedval = seedval + 17
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

    radirr1 = radirr[...]
    radconv1 = radconv[...]

    transfermode = 1 # 3D

    with open(f'out/{out_dir}/config.txt', 'w') as fh:
        seedval = seedval + 17
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
    os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log2.txt')

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

    radirr2 = radirr[...]
    radconv2 = radconv[...]

    indy = 0
    # data_plot = radirr1[:, indy, :, 0, 0, 0] + radirr1[:, indy, :, 1, 0, 0]
    # data_plot2 = radirr2[:, indy, :, 0, 0, 0] + radirr2[:, indy, :, 1, 0, 0]
    # data_plot += (radirr1[:, indy, :, 0, 0, 0] + radirr1[:, indy, :, 1, 0, 0]) * fsol * wgt * absscale
    # data_plot2 += (radirr2[:, indy, :, 0, 0, 0] + radirr2[:, indy, :, 1, 0, 0]) * fsol * wgt * absscale
    heating_factor = 86400.0 / (1004.0 * 1.1)  # to convert W/m^3 to degC/day
    data_plot += radconv1[:, indy, :, 0] * fsol * wgt * absscale * heating_factor
    data_plot2 += radconv2[:, indy, :, 0] * fsol * wgt * absscale * heating_factor
    # heating_factor2 = 86400.0 / 4e6  # to convert W/m^2 to degC/day
    data_plotb += (radirr1[:, indy, 0, 0, 0, 0] + radirr1[:, indy, 0, 1, 0, 0]) * fsol * wgt * absscale
    data_plotb2 += (radirr2[:, indy, 0, 0, 0, 0] + radirr2[:, indy, 0, 1, 0, 0]) * fsol * wgt * absscale

print('Mean surface irradiance (ICA): {:.6e} W/m^2'.format(np.mean(data_plotb)))
print('Mean surface irradiance (3D RT): {:.6e} W/m^2'.format(np.mean(data_plotb2)))

# data_plot = data_plot[:-1, :]
# data_plot2 = data_plot2[:-1, :]
# nz -= 1

with open('dat/mod/lescrop_01_24-01_1.000e+00_4.409e+02.txt', 'r') as fh:
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

surf_max = max(np.max(data_plotb), np.max(data_plotb2))

# norm = mcolors.Normalize(vmin=0., vmax=data_max)
norm = mcolors.LogNorm(vmin=data_max*1e-4, vmax=data_max)

# fig = plt.figure(figsize=(7, 3))
fig = plt.figure(figsize=(3.5, 4))
# ax = fig.add_subplot(1, 2, 1)
# ax = fig.add_axes([0.08, 0.15, 0.35, 0.75])
ax = fig.add_axes([0.08, 0.50, 0.72, 0.4])
m = ax.pcolormesh(xx * 1e-3, zz * 1e-3, data_plot.T, cmap='cividis', norm=norm)
lc = LineCollection(segs, colors='k', linewidths=1.0)
ax.add_collection(lc)
# ax.set_aspect('equal')
# ax.set_xlabel('X (km)')
ax.set_ylabel('Z (km)')
ax.set_xticklabels([])
# ax.set_title('(a) ICA')
ax.text(0.03, 0.95, '(a) ICA', transform=ax.transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='square', facecolor='white', alpha=0.8, edgecolor='none'))

pos = ax.get_position()
left, bottom, width, height = pos.x0, pos.y0, pos.width, pos.height
twin_height = height * 0.15  # adjust fraction to make it "short"
ax_twin = fig.add_axes([left, bottom, width, twin_height], sharex=ax)
ax_twin.set_facecolor('none')
ax_twin.patch.set_visible(False)
x_centers_km = (np.linspace(0.0, nx * dx, nx) + dx/2.) * 1e-3
mb = ax_twin.plot(x_centers_km, data_plotb, color='red', marker='o', linestyle='-', markersize=2, lw=1, label='地表放射照度')
ax_twin.set_ylim(0, surf_max*1.1)
# ax_twin.yaxis.set_ticklabels([])
ax_twin.spines['top'].set_visible(False)
ax_twin.yaxis.tick_right()
ax_twin.set_zorder(ax.get_zorder() + 1)
ax_twin.set_ylabel('(W/m²)', ha='left', rotation=0, fontsize=9)
ax_twin.yaxis.set_label_coords(1.01, 2.)
ax_twin.tick_params(axis='y', color='red', labelcolor='red')
ax_twin.yaxis.label.set_color('red')

# fig.colorbar(m, ax=ax, label='Irradiance (W/m²)', orientation='horizontal', pad=0.2, shrink=0.65)
ax_cb = fig.add_axes([0.95, 0.2, 0.02, 0.5])  # [left, bottom, width, height] in figure coords
# cbar = fig.colorbar(m, cax=ax_cb, label='放射加熱率 (W/m³)', orientation='vertical')
cbar = fig.colorbar(m, cax=ax_cb, label='放射加熱率 (°C/day)', orientation='vertical')
# cbar.formatter.set_useMathText(True)
# ax2 = fig.add_subplot(1, 2, 2)
# ax2 = fig.add_axes([0.45, 0.15, 0.35, 0.75])
ax2 = fig.add_axes([0.08, 0.05, 0.72, 0.4])
m2 = ax2.pcolormesh(xx * 1e-3, zz * 1e-3, data_plot2.T, cmap='cividis', norm=norm)
lc = LineCollection(segs, colors='k', linewidths=1.0, label='雲')
ax2.add_collection(lc)
# ax2.set_aspect('equal')
ax2.set_xlabel('X (km)')
ax2.set_ylabel('Z (m)')
# ax2.set_yticklabels([])
# ax2.set_title('(b) 3D')
ax2.text(0.03, 0.95, '(b) 3D', transform=ax2.transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='square', facecolor='white', alpha=0.8, edgecolor='none'))

pos = ax2.get_position()
left, bottom, width, height = pos.x0, pos.y0, pos.width, pos.height
twin_height = height * 0.15  # adjust fraction to make it "short"
ax_twin2 = fig.add_axes([left, bottom, width, twin_height], sharex=ax2)
ax_twin2.set_facecolor('none')
ax_twin2.patch.set_visible(False)
x_centers_km = (np.linspace(0.0, nx * dx, nx) + dx/2.) * 1e-3
mb = ax_twin2.plot(x_centers_km, data_plotb2, color='red', marker='o', linestyle='-', markersize=2, lw=1, label='地表放射照度')
ax_twin2.set_ylim(0, surf_max*1.1)
ax_twin2.yaxis.set_label_position('right')
ax_twin2.spines['top'].set_visible(False)
ax_twin2.yaxis.tick_right()
ax_twin2.set_zorder(ax2.get_zorder() + 1)
ax_twin2.set_ylabel('(W/m²)', ha='left', rotation=0, fontsize=9)
ax_twin2.yaxis.set_label_coords(1.01, 2.1)
ax_twin2.tick_params(axis='y', color='red', labelcolor='red')
ax_twin2.yaxis.label.set_color('red')

# create a separate axis below the panels for the legend
handles, labels = ax2.get_legend_handles_labels()
ax_leg = fig.add_axes([0.82, 0.75, 0.12, 0.12])  # [left, bottom, width, height] in figure coords
ax_leg.axis('off')
patch = matplotlib.patches.Patch(facecolor='none', edgecolor='k', label=(labels[0] if labels else '雲'))
ax_leg.legend([patch], [patch.get_label()], loc='center left', ncol=1, frameon=False)

fig.text(0.92, -0.015, '太陽天頂角: 30°', ha='center', va='top', fontsize=10)

# fig.tight_layout()
# fig.subplots_adjust(bottom=0.18)
fig.savefig(f'out/{out_dir}/03_flux_xz_comparison.png', dpi=300, bbox_inches='tight')
