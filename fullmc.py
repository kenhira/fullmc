import os
import numpy as np
import multiprocessing as mp
import time

exe_mc = 'exe_fullmc'

def run_mc_single(wrkdir, iproc):
    work_dir = f'{wrkdir}/proc{iproc:04d}'
    cmd = f'./{exe_mc} {work_dir} > {work_dir}/log.txt'
    os.system(cmd)

class FullMC:
    def __init__(self,
                 nx = None,
                 ny = None,
                 nz = None,
                 dx = None,
                 dy = None,
                 dz = None,
                 transfermode = None, 
                 source = None,
                 swlw = None,
                 solmu = None,
                 solphi = None,
                 viewmu = None,
                 viewphi = None,
                 nphoton = None,
                 seedval = 1,
                 wgttype = 1,
                 debug = 0,
                 kext = None,
                 kabs = None,
                 gparam = None,
                 bplnk = None,
                 galb = None,
                 bgrnd = None,
                 Ncpu = None,
                 wrkdir = '.'):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.transfermode = transfermode
        self.source = source
        self.swlw = swlw
        self.solmu = solmu
        self.solphi = solphi
        self.viewmu = viewmu
        self.viewphi = viewphi
        self.nphoton = nphoton
        self.seedval = seedval
        self.wgttype = wgttype
        self.debug = debug

        if kext is None:
            self.kext = np.full((nx, ny, nz), np.nan, dtype=np.float64)
        else:
            self.kext = kext
        if kabs is None:
            self.kabs = np.full((nx, ny, nz), np.nan, dtype=np.float64)
        else:
            self.kabs = kabs
        if gparam is None:
            self.gparam = np.full((nx, ny, nz), np.nan, dtype=np.float64)
        else:
            self.gparam = gparam
        if bplnk is None:
            self.bplnk = np.full((nx, ny, nz), np.nan, dtype=np.float64)
        else:
            self.bplnk = bplnk
        if galb is None:
            self.galb = np.full((nx, ny), np.nan, dtype=np.float64)
        else:
            self.galb = galb
        if bgrnd is None:
            self.bgrnd = np.full((nx, ny), np.nan, dtype=np.float64)
        else:
            self.bgrnd = bgrnd
        if Ncpu is None or Ncpu < 1:
            self.Ncpu = 0
        else:
            self.Ncpu = min(Ncpu, mp.cpu_count())
            print(f'Using {self.Ncpu} CPU cores.')
        self.wrkdir = wrkdir
    
    def read_atmtxt(self, infile, emission=False, tair=273.0, tgrnd=300.0):
        with open(infile, 'r') as fh:
            line = fh.readline()
            nx_file, ny_file, nz_file = [int(x) for x in line.strip().split()]
            assert nx_file == self.nx
            assert ny_file == self.ny
            assert nz_file == self.nz
            for ix in range(self.nx):
                for iy in range(self.ny):
                    for iz in range(self.nz):
                        line = fh.readline()
                        parts = line.strip().split()
                        self.kext[ix, iy, iz] = float(parts[0])
                        self.kabs[ix, iy, iz] = float(parts[1])
                        self.gparam[ix, iy, iz] = float(parts[2])
            if emission:
                line = fh.readline()
                parts = line.strip().split()
                wls = [float(x) for x in parts]
                line = fh.readline()
                parts = line.strip().split()
                aplk = [float(x) for x in parts]
                wl = np.sqrt(wls[0] * wls[1])
                xxair = 1. / (wl * tair)
                bbair = np.polyval([aplk[4], aplk[3], aplk[2], aplk[1], aplk[0]], xxair)
                self.bplnk = 1./(np.exp(bbair)*wl**3.*xxair)
                xxgrnd = 1. / (wl * tgrnd)
                bbgrnd = np.polyval([aplk[4], aplk[3], aplk[2], aplk[1], aplk[0]], xxgrnd)
                self.bgrnd = 1./(np.exp(bbgrnd)*wl**3.*xxgrnd)
    
    def write_config(self):
        os.system(f'mkdir -p {self.wrkdir}')
        if self.Ncpu < 1:
            self.write_config_single(self.wrkdir)
        else:
            for iproc in range(self.Ncpu):
                work_dir = f'{self.wrkdir}/proc{iproc:04d}'
                os.system(f'mkdir -p {work_dir}')
                self.write_config_single(work_dir)

    def write_config_single(self, workdir):
        with open(f'{workdir}/config.txt', 'w') as fh:
            self.seedval = self.seedval + 11
            fh.write("%d %d %d\n" % (self.nx, self.ny, self.nz))
            fh.write("%g %g %g\n" % (self.dx, self.dy, self.dz))
            fh.write("%d %d\n" % (self.source, self.swlw))
            fh.write("%d\n" % self.transfermode)
            fh.write("%g %g\n" % (self.solmu, self.solphi))
            fh.write("%g %g\n" % (self.viewmu, self.viewphi))
            fh.write("%d\n" % int(self.nphoton))
            fh.write("%d\n" % self.seedval)
            fh.write("%d\n" % self.wgttype)
            fh.write("%d\n" % self.debug)
            for ix in range(self.nx):
                for iy in range(self.ny):
                    for iz in range(self.nz):
                        if self.kext[ix,iy,iz] != self.kext[ix,iy,iz]:
                            raise ValueError("kext contains NaN values at index (%d,%d,%d)." % (ix,iy,iz))
                        if self.kabs[ix,iy,iz] != self.kabs[ix,iy,iz]:
                            raise ValueError("kabs contains NaN values at index (%d,%d,%d)." % (ix,iy,iz))
                        if self.gparam[ix,iy,iz] != self.gparam[ix,iy,iz]:
                            raise ValueError("gparam contains NaN values at index (%d,%d,%d)." % (ix,iy,iz))
                        if self.bplnk[ix,iy,iz] != self.bplnk[ix,iy,iz]:
                            raise ValueError("bplnk contains NaN values at index (%d,%d,%d)." % (ix,iy,iz))
                        fh.write("%15.6e %15.6e %15.6e %15.6e\n" % (self.kext[ix,iy,iz], self.kabs[ix,iy,iz], self.gparam[ix,iy,iz], self.bplnk[ix,iy,iz]))
                    if self.galb[ix,iy] != self.galb[ix,iy]:
                        raise ValueError("galb contains NaN values at index (%d,%d)." % (ix,iy))
                    if self.bgrnd[ix,iy] != self.bgrnd[ix,iy]:
                        raise ValueError("bgrnd contains NaN values at index (%d,%d)." % (ix,iy))
                    fh.write("%15.6e %15.6e\n" % (self.galb[ix,iy], self.bgrnd[ix,iy]))
        
    def run_mc(self):
        print("writing config file...")
        self.write_config()
        print("running FullMC...")
        start_time = time.time()
        if self.Ncpu < 1:
            cmd = f'./{exe_mc} {self.wrkdir} > {self.wrkdir}/log.txt'
            os.system(cmd)
        else:
            with mp.Pool(processes=self.Ncpu) as pool:
                pool.starmap(run_mc_single, [(self.wrkdir, iproc) for iproc in range(self.Ncpu)])
        end_time = time.time()
        print(f'done in: {end_time - start_time:.4f} sec')
        # with mp.Pool(processes=Ncpu) as pool:
        #     pool.map(execute_mc, range(Ncpu))

    def read_result(self, kind=None):
        if self.Ncpu < 1:
            if kind == 'img':
                outrad, nphoton = self.read_result_img_single(self.wrkdir)
            elif kind == 'irr':
                outrad, nphoton = self.read_result_irr_single(self.wrkdir)
            elif kind == 'conv':
                outrad, nphoton = self.read_result_conv_single(self.wrkdir)
        else:
            outrad_pts = []
            npho_pts = []
            for iproc in range(self.Ncpu):
                work_dir = f'{self.wrkdir}/proc{iproc:04d}'
                if kind == 'img':
                    outrad_part, nphoton = self.read_result_img_single(work_dir)
                elif kind == 'irr':
                    outrad_part, nphoton = self.read_result_irr_single(work_dir)
                elif kind == 'conv':
                    outrad_part, nphoton = self.read_result_conv_single(work_dir)
                outrad_pts.append(outrad_part)
                npho_pts.append(nphoton)
            npho_tot = sum(npho_pts)
            outrad = np.zeros_like(outrad_pts[0], dtype=np.float64)
            for iproc in range(self.Ncpu):
                outrad[..., 0] += outrad_pts[iproc][..., 0] * npho_pts[iproc] / npho_tot
                outrad[..., 1] += outrad_pts[iproc][..., 1] * npho_pts[iproc]*(npho_pts[iproc] - 1) / (npho_tot*(npho_tot - 1))
        return outrad
    
    def read_result_img_single(self, workdir):
        with open(f'{workdir}/outradimg.txt', 'r') as fh:
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
        return radimg, nphoton

    def read_result_irr_single(self, workdir):
        with open(f'{workdir}/outradirr.txt', 'r') as fh:
            header = fh.readline()  # skip header
            dims_line = fh.readline()
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
        return radirr, nphoton
    
    def read_result_conv_single(self, workdir):
        with open(f'{workdir}/outradconv.txt', 'r') as fh:
            header = fh.readline()  # skip header
            dims_line = fh.readline()
            nx, ny, nz, ncomp, nphoton = [int(x) for x in dims_line.strip().split()]
            
            radconv = np.zeros((nx, ny, nz, ncomp), dtype=np.float64)
            for ix in range(nx):
                for iy in range(ny):
                    for iz in range(nz):
                        line = fh.readline()
                        parts = line.strip().split()
                        radconv[ix, iy, iz, 0] = float(parts[0])
                        radconv[ix, iy, iz, 1] = float(parts[1])
        return radconv, nphoton