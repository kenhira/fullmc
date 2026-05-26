import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from matplotlib.collections import LineCollection
import glob
import matplotlib
# matplotlib.rc('font', family='Hiragino Sans')

exe_mc = 'exe_fullmc'

out_dir_base = 'ex_05_fluxxzlw'

def execute_mc(iproc):
    out_dir = f'{out_dir_base}/proc{iproc:04d}'
    os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log.txt')
def execute_mc2(iproc):
    out_dir = f'{out_dir_base}/proc{iproc:04d}'
    os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log2.txt')

if __name__ == "__main__":

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
    source = 3 # Grid side sources
    swlw = 0 # LW
    # nphoton = 50
    nphoton = 500
    # nphoton = 4000
    # nphoton = 40000

    # source = 0 # TOA Direct
    # # source = 1 # TOA Diffuse (Lambertian)
    # # source = 4 # TOA detector
    # swlw = 1 # SW
    # nphoton = 2e1
    # # nphoton = 3e2
    # # nphoton = 5e3
    # # nphoton = 2e4

    # solmu = 1.0
    # solmu = np.sqrt(2.0) / 2.0
    solmu = np.cos(np.radians(30.0))
    # solmu = 0.5
    solphi = 0.5 * np.pi

    viewmu = 1.0
    viewphi = 0.0

    seedval = 1237

    ncpu = 1
    # ncpu = 10
    # ncpu = 32
    # ncpu = 128

    Ncpu = min(ncpu, mp.cpu_count())
    print(f'Using {Ncpu} CPU cores for MC simulation.')

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
    # bplnk = np.zeros((nx, ny, nz), dtype=np.float64)
    # bgrnd = np.zeros((nx, ny), dtype=np.float64)
    # bplnk[:, :, :] = 1.0
    # bgrnd[:, :] = 1.0
    galb = np.zeros((nx, ny), dtype=np.float64)
    galb[:, :] = 0.1

    xx, zz = np.meshgrid(np.linspace(0.0, nx * dx, nx) + dx/2., np.linspace(0.0, nz * dz, nz) + dz/2.)

    tair = np.zeros((nx, ny, nz), dtype=np.float64)
    tgrnd = np.zeros((nx, ny), dtype=np.float64)
    # tair[:, :, :] = 250.0
    tair[:, :, :] = 300.0 - 6.5e-3 * zz.T[:, np.newaxis, :]
    tgrnd[:, :] = 300.0


    # file_les = 'dat/mod/les_mod_01.txt'
    files_les = glob.glob('dat/mod/lescrop_01_10-*.txt')
    # files_les = glob.glob('dat/mod/lescrop_01_01-01_*.txt')
    # files_les = sorted(glob.glob('dat/mod/lescrop_01_*.txt'))

    data_plot = np.zeros((nx, nz), dtype=np.float64)
    data_plot2 = np.zeros((nx, nz), dtype=np.float64)
    data_plotb = np.zeros((nx, nz), dtype=np.float64)
    data_plotb2 = np.zeros((nx, nz), dtype=np.float64)
    data_plot_var = np.zeros((nx, nz), dtype=np.float64)
    data_plot2_var = np.zeros((nx, nz), dtype=np.float64)
    data_plotb_var = np.zeros((nx, nz), dtype=np.float64)
    data_plotb2_var = np.zeros((nx, nz), dtype=np.float64)

    for file_les in files_les:
        bnd = int(os.path.basename(file_les).split('_')[2].split('-')[0])
        if bnd >= 15:
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
            line = fh.readline()
            parts = line.strip().split()
            wls = [float(x) for x in parts]
            line = fh.readline()
            parts = line.strip().split()
            aplk = [float(x) for x in parts]
        
        wl = np.sqrt(wls[0] * wls[1])
        xxair = 1. / (wl * tair)
        bbair = np.polyval([aplk[4], aplk[3], aplk[2], aplk[1], aplk[0]], xxair)
        bplnk = 1./(np.exp(bbair)*wl**3.*xxair)
        xxgrnd = 1. / (wl * tgrnd)
        bbgrnd = np.polyval([aplk[4], aplk[3], aplk[2], aplk[1], aplk[0]], xxgrnd)
        bgrnd = 1./(np.exp(bbgrnd)*wl**3.*xxgrnd)
        
        absscale = 1.0 #np.exp(-np.mean(kabs[:, :, -1])*8.e3)
        
        wgt = float(os.path.basename(file_les).split('_')[-2])
        fsol = float(os.path.basename(file_les).split('_')[-1].split('.txt')[0])

        os.system(f'mkdir -p out/{out_dir_base}')
        for iproc in range(Ncpu):
            out_dir = f'{out_dir_base}/proc{iproc:04d}'
            os.system(f'mkdir -p out/{out_dir}')
        # os.system(f'mkdir -p out/{out_dir}')

        transfermode = 0 # ICA

        for iproc in range(Ncpu):
            out_dir = f'{out_dir_base}/proc{iproc:04d}'
            seedval = seedval + 17
            with open(f'out/{out_dir}/config.txt', 'w') as fh:
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
        # os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log.txt')
        # for iproc in range(Ncpu):
        #     out_dir = f'{out_dir_base}/proc{iproc:04d}'
        #     os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log.txt')
        # def execute_mc(iproc):
        #     out_dir = f'{out_dir_base}/proc{iproc:04d}'
        #     os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log.txt')
        
        with mp.Pool(processes=Ncpu) as pool:
            pool.map(execute_mc, range(Ncpu))

        for iproc in range(Ncpu):
            out_dir = f'{out_dir_base}/proc{iproc:04d}'
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

            radirr1prc = radirr[...]
            radconv1prc = radconv[...]

            radirr1 = radirr1prc / float(Ncpu) if iproc == 0 else radirr1 + radirr1prc / float(Ncpu)
            radconv1 = radconv1prc / float(Ncpu) if iproc == 0 else radconv1 + radconv1prc / float(Ncpu)

        transfermode = 1 # 3D

        for iproc in range(Ncpu):
            out_dir = f'{out_dir_base}/proc{iproc:04d}'
            seedval = seedval + 17
            with open(f'out/{out_dir}/config.txt', 'w') as fh:
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
        # os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log2.txt')
        # for iproc in range(Ncpu):
        #     out_dir = f'{out_dir_base}/proc{iproc:04d}'
        #     os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log2.txt')
        # def execute_mc2(iproc):
        #     out_dir = f'{out_dir_base}/proc{iproc:04d}'
        #     os.system(f'./{exe_mc} out/{out_dir} > out/{out_dir}/log2.txt')
        
        with mp.Pool(processes=Ncpu) as pool:
            pool.map(execute_mc2, range(Ncpu))

        for iproc in range(Ncpu):
            out_dir = f'{out_dir_base}/proc{iproc:04d}'
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

            radirr2prc = radirr[...]
            radconv2prc = radconv[...]

            radirr2 = radirr2prc / float(Ncpu) if iproc == 0 else radirr2 + radirr2prc / float(Ncpu)
            radconv2 = radconv2prc / float(Ncpu) if iproc == 0 else radconv2 + radconv2prc / float(Ncpu)

        indy = 0
        # data_plot = radirr1[:, indy, :, 0, 0, 0] + radirr1[:, indy, :, 1, 0, 0]
        # data_plot2 = radirr2[:, indy, :, 0, 0, 0] + radirr2[:, indy, :, 1, 0, 0]
        # data_plotb += (radirr1[:, indy, :, 0, 0, 0] + radirr1[:, indy, :, 1, 0, 0]) * wgt * absscale
        # data_plotb2 += (radirr2[:, indy, :, 0, 0, 0] + radirr2[:, indy, :, 1, 0, 0]) * wgt * absscale
        data_plotb += (radirr1[:, indy, :, 0, 0, 0]) * wgt * absscale
        data_plotb2 += (radirr2[:, indy, :, 0, 0, 0]) * wgt * absscale
        heating_factor = 86400.0 / (1004.0 * 1.1)  # to convert W/m^3 to degC/day
        data_plot += radconv1[:, indy, :, 0] * wgt * absscale * heating_factor
        data_plot2 += radconv2[:, indy, :, 0] * wgt * absscale * heating_factor

        data_plot_var += nphoton*(nphoton-1)*radconv1[:, indy, :, 1] * (wgt * absscale) **2
        data_plot2_var += nphoton*(nphoton-1)*radconv2[:, indy, :, 1] * (wgt * absscale) **2
        data_plotb_var += nphoton*(nphoton-1)*(radirr1[:, indy, :, 0, 0, 1]) * (wgt * absscale * heating_factor)**2
        data_plotb2_var += nphoton*(nphoton-1)*(radirr2[:, indy, :, 0, 0, 1]) * (wgt * absscale * heating_factor)**2
    
    data_plot_std = np.sqrt(data_plot_var / (Ncpu * nphoton * (Ncpu * nphoton - 1)))
    data_plot2_std = np.sqrt(data_plot2_var / (Ncpu * nphoton * (Ncpu * nphoton - 1)))
    data_plotb_std = np.sqrt(data_plotb_var / (Ncpu * nphoton * (Ncpu * nphoton - 1)))
    data_plotb2_std = np.sqrt(data_plotb2_var / (Ncpu * nphoton * (Ncpu * nphoton - 1)))

    print('Mean surface rad conv (ICA): {:.6e}'.format(np.mean(data_plot[:, 0])))
    print('Mean surface rad conv (3D RT): {:.6e}'.format(np.mean(data_plot2[:, 0])))

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

    data_max = max(np.max(data_plot), np.max(data_plot2))
    data_min = min(np.min(data_plot), np.min(data_plot2))
    data_abs = max(abs(data_min), abs(data_max))
    # norm = mcolors.Normalize(vmin=0., vmax=data_max)
    # norm = mcolors.Normalize(vmin=-data_abs, vmax=data_abs)
    norm = mcolors.SymLogNorm(linthresh=10**(np.floor(np.log10(data_abs) - 1.)), vmin=-data_abs, vmax=data_abs)
    # norm = mcolors.LogNorm(vmin=data_max*1e-4, vmax=data_max)

    fig = plt.figure(figsize=(7, 3))
    # ax = fig.add_subplot(1, 2, 1)
    ax = fig.add_axes([0.1, 0.15, 0.35, 0.75])
    m = ax.pcolormesh(xx * 1e-3, zz * 1e-3, data_plot.T, cmap='coolwarm', norm=norm)
    lc = LineCollection(segs, colors='k', linewidths=1.0)
    ax.add_collection(lc)
    ax.set_aspect('equal')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Z (km)')
    ax.set_title('(a) ICA')
    # fig.colorbar(m, ax=ax, label='Irradiance (W/m²)', orientation='horizontal', pad=0.2, shrink=0.65)
    ax_cb = fig.add_axes([0.87, 0.2, 0.02, 0.5])  # [left, bottom, width, height] in figure coords
    # cbar = fig.colorbar(m, cax=ax_cb, label='放射加熱率 (W/m³)', orientation='vertical')
    # cbar = fig.colorbar(m, cax=ax_cb, label='放射加熱率 (℃/day)', orientation='vertical')
    cbar = fig.colorbar(m, cax=ax_cb, label='Radiative Heating (℃/day)', orientation='vertical')
    # cbar.formatter.set_useMathText(True)
    # ax2 = fig.add_subplot(1, 2, 2)
    ax2 = fig.add_axes([0.5, 0.15, 0.35, 0.75])
    m2 = ax2.pcolormesh(xx * 1e-3, zz * 1e-3, data_plot2.T, cmap='coolwarm', norm=norm)
    # lc = LineCollection(segs, colors='k', linewidths=1.0, label='雲')
    lc = LineCollection(segs, colors='k', linewidths=1.0, label='cloud')
    ax2.add_collection(lc)
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (km)')
    # ax2.set_ylabel('Z (m)')
    ax2.set_yticklabels([])
    ax2.set_title('(b) 3D')
    # create a separate axis below the panels for the legend
    handles, labels = ax2.get_legend_handles_labels()
    ax_leg = fig.add_axes([0.86, 0.75, 0.12, 0.12])  # [left, bottom, width, height] in figure coords
    ax_leg.axis('off')
    # patch = matplotlib.patches.Patch(facecolor='none', edgecolor='k', label=(labels[0] if labels else '雲'))
    patch = matplotlib.patches.Patch(facecolor='none', edgecolor='k', label=(labels[0] if labels else 'cloud'))
    ax_leg.legend([patch], [patch.get_label()], loc='center left', ncol=1, frameon=False)

    # fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(f'out/{out_dir_base}/05_fconv_xz_comparison.png', dpi=300)


    data_max = max(np.max(data_plotb), np.max(data_plotb2))
    data_min = min(np.min(data_plotb), np.min(data_plotb2))
    data_abs = max(abs(data_min), abs(data_max))
    norm = mcolors.Normalize(vmin=0., vmax=data_max)
    
    fig = plt.figure(figsize=(7, 3))
    # ax = fig.add_subplot(1, 2, 1)
    ax = fig.add_axes([0.1, 0.15, 0.35, 0.75])
    m = ax.pcolormesh(xx * 1e-3, zz * 1e-3, data_plotb.T, cmap='cividis', norm=norm)
    lc = LineCollection(segs, colors='k', linewidths=1.0)
    ax.add_collection(lc)
    ax.set_aspect('equal')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Z (km)')
    ax.set_title('(a) ICA')
    # fig.colorbar(m, ax=ax, label='Irradiance (W/m²)', orientation='horizontal', pad=0.2, shrink=0.65)
    ax_cb = fig.add_axes([0.87, 0.2, 0.02, 0.5])  # [left, bottom, width, height] in figure coords
    # cbar = fig.colorbar(m, cax=ax_cb, label='放射フラックス (W/m²)', orientation='vertical')
    cbar = fig.colorbar(m, cax=ax_cb, label='Radiative Flux (W/m²)', orientation='vertical')
    # cbar.formatter.set_useMathText(True)
    # ax2 = fig.add_subplot(1, 2, 2)
    ax2 = fig.add_axes([0.5, 0.15, 0.35, 0.75])
    m2 = ax2.pcolormesh(xx * 1e-3, zz * 1e-3, data_plotb2.T, cmap='cividis', norm=norm)
    # lc = LineCollection(segs, colors='k', linewidths=1.0, label='雲')
    lc = LineCollection(segs, colors='k', linewidths=1.0, label='cloud')
    ax2.add_collection(lc)
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (km)')
    # ax2.set_ylabel('Z (m)')
    ax2.set_yticklabels([])
    ax2.set_title('(b) 3D')
    # create a separate axis below the panels for the legend
    handles, labels = ax2.get_legend_handles_labels()
    ax_leg = fig.add_axes([0.86, 0.75, 0.12, 0.12])  # [left, bottom, width, height] in figure coords
    ax_leg.axis('off')
    # patch = matplotlib.patches.Patch(facecolor='none', edgecolor='k', label=(labels[0] if labels else '雲'))
    patch = matplotlib.patches.Patch(facecolor='none', edgecolor='k', label=(labels[0] if labels else 'cloud'))
    ax_leg.legend([patch], [patch.get_label()], loc='center left', ncol=1, frameon=False)

    # fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(f'out/{out_dir_base}/05_flux_xz_comparison.png', dpi=300)
    

    data_max = max(np.max(data_plot_std), np.max(data_plot2_std))
    data_min = min(np.min(data_plot_std), np.min(data_plot2_std))
    data_abs = max(abs(data_min), abs(data_max))
    # norm = mcolors.Normalize(vmin=0., vmax=data_max)
    # norm = mcolors.Normalize(vmin=-data_abs, vmax=data_abs)
    # norm = mcolors.SymLogNorm(linthresh=10**(np.floor(np.log10(data_abs) - 1.)), vmin=-data_abs, vmax=data_abs)
    # norm = mcolors.LogNorm(vmin=data_max*1e-4, vmax=data_max)
    norm = mcolors.LogNorm(vmin=data_min, vmax=data_max)

    fig = plt.figure(figsize=(7, 3))
    # ax = fig.add_subplot(1, 2, 1)
    ax = fig.add_axes([0.1, 0.15, 0.35, 0.75])
    m = ax.pcolormesh(xx * 1e-3, zz * 1e-3, data_plot_std.T, cmap='viridis', norm=norm)
    lc = LineCollection(segs, colors='k', linewidths=1.0)
    ax.add_collection(lc)
    ax.set_aspect('equal')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Z (km)')
    ax.set_title('(a) ICA')
    # fig.colorbar(m, ax=ax, label='Irradiance (W/m²)', orientation='horizontal', pad=0.2, shrink=0.65)
    ax_cb = fig.add_axes([0.87, 0.2, 0.02, 0.5])  # [left, bottom, width, height] in figure coords
    # cbar = fig.colorbar(m, cax=ax_cb, label='放射加熱率 (W/m³)', orientation='vertical')
    # cbar = fig.colorbar(m, cax=ax_cb, label='放射加熱率 (℃/day)', orientation='vertical')
    cbar = fig.colorbar(m, cax=ax_cb, label='Radiative Heating (℃/day)', orientation='vertical')
    # cbar.formatter.set_useMathText(True)
    # ax2 = fig.add_subplot(1, 2, 2)
    ax2 = fig.add_axes([0.5, 0.15, 0.35, 0.75])
    m2 = ax2.pcolormesh(xx * 1e-3, zz * 1e-3, data_plot2_std.T, cmap='viridis', norm=norm)
    # lc = LineCollection(segs, colors='k', linewidths=1.0, label='雲')
    lc = LineCollection(segs, colors='k', linewidths=1.0, label='cloud')
    ax2.add_collection(lc)
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (km)')
    # ax2.set_ylabel('Z (m)')
    ax2.set_yticklabels([])
    ax2.set_title('(b) 3D')
    # create a separate axis below the panels for the legend
    handles, labels = ax2.get_legend_handles_labels()
    ax_leg = fig.add_axes([0.86, 0.75, 0.12, 0.12])  # [left, bottom, width, height] in figure coords
    ax_leg.axis('off')
    # patch = matplotlib.patches.Patch(facecolor='none', edgecolor='k', label=(labels[0] if labels else '雲'))
    patch = matplotlib.patches.Patch(facecolor='none', edgecolor='k', label=(labels[0] if labels else 'cloud'))
    ax_leg.legend([patch], [patch.get_label()], loc='center left', ncol=1, frameon=False)

    # fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(f'out/{out_dir_base}/05_fconv_xz_comparison2.png', dpi=300)

    
