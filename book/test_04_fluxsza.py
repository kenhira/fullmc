import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from matplotlib.collections import LineCollection
import glob
import matplotlib
matplotlib.rc('font', family='Hiragino Sans')

exe_mc = 'exe_fullmc'

out_dir = 'test_04_fluxsza'

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
nphoton = 5e3
# nphoton = 2e4

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

solmu_list = np.cos(np.radians(np.arange(0, 81, 10)))
irravg_ica = np.zeros_like(solmu_list)
irravg_3d = np.zeros_like(solmu_list)

for solmu_idx, solmu in enumerate(solmu_list):

    # file_les = 'dat/mod/les_mod_01.txt'
    # files_les = glob.glob('dat/mod/lescrop_01_24-01_*.txt')
    # files_les = glob.glob('dat/mod/lescrop_01_24*.txt')
    files_les = sorted(glob.glob('dat/mod/lescrop_01_*.txt'))

    data_plot = np.zeros((nx, nz), dtype=np.float64)
    data_plot2 = np.zeros((nx, nz), dtype=np.float64)
    data_plotb = np.zeros((nx, nz), dtype=np.float64)
    data_plotb2 = np.zeros((nx, nz), dtype=np.float64)

    for file_les in files_les:
        bnd = int(os.path.basename(file_les).split('_')[2].split('-')[0])
        if bnd < 15:
            continue
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
        data_plot += (radirr1[:, indy, :, 0, 0, 0] + radirr1[:, indy, :, 1, 0, 0]) * fsol * wgt * absscale
        data_plot2 += (radirr2[:, indy, :, 0, 0, 0] + radirr2[:, indy, :, 1, 0, 0]) * fsol * wgt * absscale
        # heating_factor = 86400.0 / (1004.0 * 1.1)  # to convert W/m^3 to degC/day
        # data_plot += radconv1[:, indy, :, 0] * fsol * wgt * absscale * heating_factor
        # data_plot2 += radconv2[:, indy, :, 0] * fsol * wgt * absscale * heating_factor
        # heating_factor2 = 86400.0 / 4e6  # to convert W/m^2 to degC/day
        data_plotb += (radirr1[:, indy, :, 1, 1, 0]) * fsol * wgt * absscale
        data_plotb2 += (radirr2[:, indy, :, 1, 1, 0]) * fsol * wgt * absscale

    # print('Mean surface irradiance (ICA): {:.6e} W/m^2'.format(np.mean(data_plot[0, :, :])))
    # print('Mean surface irradiance (3D RT): {:.6e} W/m^2'.format(np.mean(data_plot2[0, :, :])))
    irravg_ica[solmu_idx] = np.mean(data_plotb[:, -1]) / np.mean(data_plot[:, -1])
    irravg_3d[solmu_idx] = np.mean(data_plotb2[:, -1]) / np.mean(data_plot2[:, -1])

fig = plt.figure(figsize=(4,3))
plt.plot(np.arccos(solmu_list)*180.0/np.pi, irravg_ica, 'o-', label='ICA')
plt.plot(np.arccos(solmu_list)*180.0/np.pi, irravg_3d, 's--', label='3D RT')
plt.ylim(0.0, 1.0)
# plt.xlabel('Solar Zenith Angle (deg)')
# plt.ylabel('Mean Reflectance')
plt.xlabel('太陽天頂角 (deg)')
plt.ylabel('領域平均反射率')
# plt.title('TOA Reflectance vs Solar Zenith Angle')
plt.legend()
plt.tight_layout()
plt.savefig(f'out/{out_dir}/04_surface_irradiance_vs_sza.png', dpi=300)
plt.close()