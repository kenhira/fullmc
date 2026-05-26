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

out_dir = 'test_00_mctestxz'

# debug = 0
# debug = 1
# nphoton = 1e4
debug = 2
nphoton = 1



# nx = 20
# ny = 1
# nz = 15
nx = 28
ny = 1
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

source = 0 # TOA Direct
# source = 1 # TOA Diffuse (Lambertian)
# source = 4 # TOA detector
swlw = 1 # SW

# solmu = 1.0
# solmu = np.sqrt(2.0) / 2.0
solmu = np.cos(np.radians(30.0))
# solmu = 0.5
solphi = 0.5 * np.pi

viewmu = 1.0
viewphi = 0.0

# seedval = 1242
seedval = 1284

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
# kext[:, :, :] = 1e-4
# kext[3:8, 0:ny, -10:-3] = 1.1e-1
# kabs = np.zeros((nx, ny, nz), dtype=np.float64)
# kabs[:, :, :] = 5e-5
# kabs[3:8, 0:ny, -10:-3] = 1e-4
# gparam = np.zeros((nx, ny, nz), dtype=np.float64)
# gparam[:, :, :] = 0.0001
# gparam[3:8, 0:ny, -10:-3] = 0.85
bplnk = np.zeros((nx, ny, nz), dtype=np.float64)
bgrnd = np.zeros((nx, ny), dtype=np.float64)
galb = np.zeros((nx, ny), dtype=np.float64)
galb[:, :] = 0.0

with open(f'dat/mod/les_mod_01.txt', 'r') as fh:
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

with open(f'out/{out_dir}/photon_trajectory.txt', 'r') as fh:
    header = fh.readline()  # skip header
    traj_data = [[]]
    iph = 0
    ixyzaph_prev = (0, 0, 0, 0, 0)
    for line in fh:
        parts = line.strip().split()
        ixyzaph = (int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
        if ixyzaph != ixyzaph_prev:
            if traj_data:
                traj_data.append([])
                iph += 1
            ixyzaph_prev = ixyzaph
        it = int(parts[5])
        xph = float(parts[8])
        yph = float(parts[9])
        zph = float(parts[10])
        xdir = float(parts[14])
        ydir = float(parts[15])
        zdir = float(parts[16])
        wgt = float(parts[17])
        if transfermode == 1:
            if xph >= nx*dx - 1e-10:
                traj_data[iph].append((0., yph, zph, wgt))
                traj_data[iph].append((np.nan, yph, zph, wgt))
            elif xph < 0.0 + 1e-10:
                traj_data[iph].append((nx*dx, yph, zph, wgt))
                traj_data[iph].append((np.nan, yph, zph, wgt))
            elif yph >= ny*dy - 1e-10:
                traj_data[iph].append((xph, 0., zph, wgt))
                traj_data[iph].append((xph, np.nan, zph, wgt))
            elif yph < 0.0 + 1e-10:
                traj_data[iph].append((xph, ny*dy, zph, wgt))
                traj_data[iph].append((xph, np.nan, zph, wgt))
        else: # ICA
            if np.mod(xph, dx) < 1e-10 or np.mod(xph, dx) >= dx - 1e-10:
                if xdir < 0.0:
                    traj_data[iph].append(((xph // dx)*dx, yph, zph, wgt))
                elif xdir >= 0.0:
                    traj_data[iph].append(((xph // dx + 1)*dx, yph, zph, wgt))
                traj_data[iph].append((np.nan, np.nan, np.nan, wgt))
            if np.mod(yph, dy) < 1e-10 or np.mod(yph, dy) >= dy - 1e-10:
                if ydir < 0.0:
                    traj_data[iph].append((xph, (yph // dy)*dy, zph, wgt))
                elif ydir >= 0.0:
                    traj_data[iph].append((xph, (yph // dy + 1)*dy, zph, wgt))
                traj_data[iph].append((np.nan, np.nan, np.nan, wgt))
        traj_data[iph].append((xph, yph, zph, wgt))
        # traj_data.append([float(x) if i >= 6 else int(x) for i, x in enumerate(parts)])

import warnings
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
def colored_line_between_pts(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified between (x, y) points by a third value.

    It does this by creating a collection of line segments between each pair of
    neighboring points. The color of each segment is determined by the
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should have a size one less than that of x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Check color array size (LineCollection still works, but values are unused)
    if len(c) != len(x) - 1:
        warnings.warn(
            "The c argument should have a length one less than the length of x and y. "
            "If it has the same length, use the colored_line function instead."
        )

    # Create a set of line segments so that we can color them individually
    # This creates the points as an N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, **lc_kwargs)

    # Set the values used for colormapping
    lc.set_array(c)

    return ax.add_collection(lc)

xx, zz = np.meshgrid(np.linspace(0, nx*dx*1e-3, nx+1), np.linspace(0, nz*dz*1e-3, nz+1), indexing='ij')

fig = plt.figure(figsize=(5,4))
ax = fig.add_axes([0.11, 0.18, 0.65, 0.7])
ax.set_aspect('equal')
kx1 = ax.pcolormesh(xx, zz, kext[:,0,:], norm=mcolors.LogNorm(vmin=kext[kext>0].min(), vmax=kext.max()), cmap='Blues_r')
# draw greyscale background and overlay hatch patterns by value
data = 1. - kabs[:,0,:]/kext[:,0,:]

n_hatches = 5
bins = np.linspace(np.nanmin(data), np.nanmax(data), n_hatches+1)
patterns = ['', '.', '..', '...','....']  # increasing hatch density

# draw hatch-patterned rectangles over cells according to binned values
# for ix in range(nx):
#     for iz in range(nz):
#         val = data[ix, iz]
#         if np.isnan(val):
#             continue
#         idx = np.digitize(val, bins) - 1
#         idx = max(0, min(idx, n_hatches-1))
#         pat = patterns[idx]
#         if pat:
#             rect = Rectangle((ix*dx*1e-3, iz*dz*1e-3), dx*1e-3, dz*1e-3, facecolor='none', edgecolor='k', hatch=pat, linewidth=0.0)
#             ax.add_patch(rect)

# plot photon trajectories
for iph in range(0, len(traj_data), 2):
    traj = np.array(traj_data[iph])
    lines = colored_line_between_pts(traj[:,0]*1e-3, traj[:,2]*1e-3, traj[:,3], ax, linewidth=1, cmap=plt.cm.rainbow)
    ax.scatter(traj[:,0]*1e-3, traj[:,2]*1e-3, c=traj[:,3], s=1.5, cmap=plt.cm.rainbow)

ax.set_yticks(np.arange(0, nz*dz*1e-3+1, 0.5))
ax.set_xlabel('X (km)')
ax.set_ylabel('Z (km)')
ax.set_xlim(0, nx*dx*1e-3)
ax.set_ylim(0, nz*dz*1e-3)
for ix in range(nx+1):
    ax.axvline(ix*dx*1e-3, color='gray', linestyle='--', linewidth=0.5)
for iz in range(nz+1):
    ax.axhline(iz*dz*1e-3, color='gray', linestyle='--', linewidth=0.5)
# colorbar for extinction
ax_cb1 = fig.add_axes([0.82, 0.56, 0.03, 0.23])  # small axes on the right
cbar1 = fig.colorbar(kx1, cax=ax_cb1)
ax_cb1.set_title(' '*10 +r'消散係数 ($\rm m^{-1}$)', fontsize=10)
 
# # hatch "colorbar" (legend) built from stacked rectangles showing bin ranges and hatch patterns
# ax_hatch = fig.add_axes([0.82, 0.20, 0.03, 0.23])  # small axes on the right
# ax_hatch.set_xlim(0, 1)
# ax_hatch.set_ylim(0, n_hatches)
# ax_hatch.axis('off')

# # draw stacked patches and tick labels
# labels = []
# for i in range(n_hatches):
#     y = i
#     pat = patterns[i]
#     # draw patch; use white face so hatch is visible
#     rect = Rectangle((0, y), 1, 1, facecolor='white', edgecolor='k', hatch=pat)
#     ax_hatch.add_patch(rect)
#     labels.append(f"{bins[i]:.2f}–{bins[i+1]:.2f}")

# # place tick labels centered on each patch
# ax_hatch_ticks = np.arange(0.5, n_hatches+0.5, 1.0)[:n_hatches]
# for y, lab in zip(ax_hatch_ticks, labels):
#     ax_hatch.text(1.2, y, lab, va='center', fontsize=8)
# ax_hatch.set_title('SSA', fontsize=10)

# optional colorbar for photon weight
# ax_cb3 = fig.add_axes([0.16, 0.13, 0.55, 0.03])  # small axes at the bottom
# cbar3 = fig.colorbar(lines, cax=ax_cb3, orientation='horizontal')
# cbar3.set_label('Photon Weight')
ax_cb3 = fig.add_axes([0.82, 0.20, 0.03, 0.23])  # small axes at the bottom
cbar3 = fig.colorbar(lines, cax=ax_cb3, orientation='vertical')
ax_cb3.set_title(' '*10 + '重み', fontsize=10, pad=10)

plt.savefig(f"out/{out_dir}/00_photon_trajectories_{'ICA' if transfermode == 0 else '3D'}.png", dpi=300)
plt.close()

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

indy = 0
data_plot = radirr[:, indy, :, 0, 0, 0]
data_plot2 = radirr[:, indy, :, 1, 0, 0]

xx, zz = np.meshgrid(np.linspace(0.0, nx * 100.0, nx), np.linspace(0.0, nz * 80.0, nz))

fig = plt.figure(figsize=(5, 4))

ax2 = fig.add_subplot(111)
m2 = ax2.pcolormesh(xx, zz, (data_plot + data_plot2).T, cmap='cividis', vmin=0.)
ax2.set_aspect('equal')
fig.colorbar(m2, ax=ax2, pad=0.05, shrink=0.65)
ax2.set_xlabel('X (m)')
ax2.set_ylabel('Z (m)')
ax2.set_title(f'Y={indy * dy:.1f} m - Global Downward Irradiance')
fig.tight_layout()

fig.savefig(f'out/{out_dir}/00_radirr_xz_downward_y{indy}.png', dpi=300)