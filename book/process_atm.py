import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import glob

def read_ctl_file(ctl_path):
    """Parse the .ctl file to extract metadata."""
    metadata = {}
    with open(ctl_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("DSET"):
                metadata['dset'] = line.split('^')[-1].strip()
            elif line.startswith("UNDEF"):
                metadata['undef'] = float(line.split()[-1])
            elif line.startswith("XDEF"):
                _, count, _, start, step = line.split()
                metadata['xdef'] = (int(count), float(start), float(step))
            elif line.startswith("YDEF"):
                _, count, _, start, step = line.split()
                metadata['ydef'] = (int(count), float(start), float(step))
            elif line.startswith("ZDEF"):
                _, count, _, start, step = line.split()
                metadata['zdef'] = (int(count), float(start), float(step))
            elif line.startswith("VARS"):
                metadata['vars'] = []
            elif line.startswith("ENDVARS"):
                break
            elif 'vars' in metadata:
                parts = line.split()
                metadata['vars'].append((parts[0], int(parts[1]), parts[2], ' '.join(parts[3:])))
    return metadata

def read_grd_file(grd_path, metadata):
    """Read the binary .grd file based on metadata from the .ctl file."""
    x_count, _, _ = metadata['xdef']
    y_count, _, _ = metadata['ydef']
    z_count, _, _ = metadata['zdef']
    num_vars = len(metadata['vars'])
    
    data = np.zeros((num_vars, z_count, y_count, x_count), dtype=np.float32)
    with open(grd_path, 'rb') as f:
        for var_idx in range(num_vars):
            for z in range(metadata['vars'][var_idx][1]):
                raw_data = f.read(x_count * y_count * 4)  # 4 bytes per float
                layer_data = np.frombuffer(raw_data, dtype='>f4')  # Big-endian float32
                layer_data = layer_data.reshape((y_count, x_count))
                data[var_idx, z, :, :] = layer_data
    return data

def load_data(ctl_path):
    """Load metadata and data from .ctl and .grd files."""
    dirpath = os.path.dirname(ctl_path)
    metadata = read_ctl_file(ctl_path)
    grd_path = os.path.join(dirpath, metadata['dset'])  # Path to the .grd file
    data = read_grd_file(grd_path, metadata)
    return (data, metadata)


def data_slice(data, metadata, varname, dx, dy, dz, axi, level, func=None):
    var_index = next(i for i, var in enumerate(metadata['vars']) if var[0] == varname)
    x_count, x_start, x_step = metadata['xdef']
    y_count, y_start, y_step = metadata['ydef']
    _, z_start, z_step = metadata['zdef']
    z_count = metadata['vars'][var_index][1]

    x_arr = np.linspace(x_start, x_start + (x_count - 1) * x_step * dx, x_count)
    y_arr = np.linspace(y_start, y_start + (y_count - 1) * y_step * dy, y_count)
    z_arr = np.linspace(z_start, z_start + (z_count - 2) * z_step * dz, z_count - 1)

    if axi == 'z':
        xx, yy = np.meshgrid(x_arr, y_arr)
        if func is not None:
            return xx, yy, func(data[var_index, level, :, :], axis=0)
        else:
            return xx, yy, data[var_index, level, :, :]
    elif axi == 'y':
        xx, zz = np.meshgrid(x_arr, z_arr)
        if func is not None:
            return xx, zz, func(data[var_index, :z_count-1, :, :][:, level, :], axis=1)
        else:
            return xx, zz, data[var_index, :z_count-1, :, :][:, level, :]
    elif axi == 'x':
        yy, zz = np.meshgrid(y_arr, z_arr)
        if func is not None:
            return yy, zz, func(data[var_index, :z_count-1, :, :][:, :, level], axis=2)
        else:
            return yy, zz, data[var_index, :z_count-1, :, :][:, :, level]


# xslice = slice(0, None)
# yslice = slice(0, None)
# zslice = slice(0, None)
xslice = slice(34, 62)
# yslice = slice(18, 46)
yslice = slice(31, 32)
zslice = slice( 0, 28)

# ctl_file = 'dat/orig/les_24-01.atm.ctl'
# ctl_files = glob.glob('dat/orig/les_*.atm.ctl')
ctl_files = glob.glob('dat/orig/shcu/les_*.atm.ctl')
# ctl_files = glob.glob('dat/orig/les_24-01.atm.ctl')

# mod_id = '01'
mod_id = '02'

with open('dat/orig/PARAG.29', 'r') as fp:
    lines = fp.readlines()
    line_len = len(lines)
    il = 0
    nsub_list = []
    weights_list = []
    while il < line_len:
        if 'number of subinterval' in lines[il].lower():
            # read number of subintervals
            il += 1
            while il < line_len and lines[il].strip() == '':
                il += 1
            try:
                n_sub = int(lines[il].strip().split()[0])
                nsub_list.append(n_sub)
            except Exception:
                il += 1
                continue
            il += 1
            # advance to "weight of subintervals"
            while il < line_len and 'weight of subinterval' not in lines[il].lower():
                il += 1
            if il >= line_len:
                break
            # collect numeric tokens until we have n_sub floats
            tokens = []
            # consider numbers possibly on the same "weight ..." line after the text
            after = lines[il].split('weight of subintervals')[-1].strip()
            if after:
                tokens.extend(after.split())
            il += 1
            while len(tokens) < n_sub and il < line_len:
                line_tokens = lines[il].strip().split()
                # stop if a clearly non-numeric header is encountered
                # if line_tokens and any(c.isalpha() for c in line_tokens[0]) and len(line_tokens) == 1:
                if line_tokens and 'Line absorption' in line_tokens[0]:
                    break
                tokens.extend(line_tokens)
                il += 1
            # parse floats
            vals = []
            for t in tokens:
                try:
                    vals.append(float(t))
                except Exception:
                    continue
            weights_list.append(vals[:n_sub])
        else:
            il += 1

with open('dat/orig/PARAPC.29', 'r') as fp:
    lines = fp.readlines()
    line_len = len(lines)
    il = 0
    fsol_list = []
    aplk_list = []
    wv_list = []
    while il < line_len:
        if lines[il].startswith(' Band boundary wave number (cm-1)'):
            il += 1
            # skip band boundary numbers
            while il < line_len and 'Planck' not in lines[il]:
                wv_list.extend([float(x) for x in lines[il].strip().split()])
                il += 1
        if lines[il].startswith('  Incident solar flux (W/m2)'):
            il += 1
            try:
                fsol = float(lines[il].strip().split()[0])
                fsol_list.append(fsol)
            except Exception:
                continue
        if lines[il].startswith('  Planck function fitting parameters'):
            il += 1
            try:
                aplk = [float(x) for x in lines[il].strip().split()]
                aplk_list.append(aplk)
            except Exception:
                continue
        il += 1

print(nsub_list)
print(weights_list)
print(fsol_list)
print(aplk_list)
print(wv_list)

# mod_file = 'dat/mod/les_mod_01.txt'
# mod_file = 'dat/mod/les_mod_02.txt'
# mod_files = ['dat/mod/les_crop-%s_%.3e_%.3e.txt' % ('%02d' % (i+1), 1.e-3, 1.e-3) for i in range(len(ctl_files))]

for ctl_file in ctl_files:
    basename = os.path.basename(ctl_file)
    bndstr = basename.split('_')[1].split('.')[0].split('-')[0]
    substr = basename.split('_')[1].split('.')[0].split('-')[1]
    weight = weights_list[int(bndstr)-1][int(substr)-1]
    fsol = fsol_list[int(bndstr)-1]
    aplk = aplk_list[int(bndstr)-1]
    wl1 = 1.e4 / wv_list[int(bndstr)-1]
    wl2 = 1.e4 / wv_list[int(bndstr)]
    mod_file = 'dat/mod/lescrop_%s_%s-%s_%.3e_%.3e.txt' % (mod_id, bndstr, substr, weight, fsol)

    data, metadata = load_data(ctl_file)
    data_mod = data[:, zslice, yslice, xslice]
    metadata_mod = metadata.copy()
    metadata_mod['xdef'] = (data_mod.shape[3], metadata['xdef'][1], metadata['xdef'][2])
    metadata_mod['ydef'] = (data_mod.shape[2], metadata['ydef'][1], metadata['ydef'][2])
    for var_index in range(len(metadata['vars'])):
        metadata_mod['vars'][var_index] = (metadata['vars'][var_index][0], data_mod.shape[1], metadata['vars'][var_index][2:])

    # print(metadata_mod)

    with open(mod_file, 'w') as fh:
        nx = data_mod.shape[3]
        ny = data_mod.shape[2]
        nz = data_mod.shape[1]
        fh.write(f"{nx} {ny} {nz}\n")
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    varname = 'extp3d'
                    var_index = next(i for i, var in enumerate(metadata['vars']) if var[0] == varname)
                    kext = data_mod[var_index, iz, iy, ix]
                    varname = 'omgp3d'
                    var_index = next(i for i, var in enumerate(metadata['vars']) if var[0] == varname)
                    omg = data_mod[var_index, iz, iy, ix]
                    kabs = kext * (1.0 - omg)
                    varname = 'apfp3d'
                    var_index = next(i for i, var in enumerate(metadata['vars']) if var[0] == varname)
                    gparam = data_mod[var_index, iz, iy, ix]
                    fh.write(f"{kext:.6e} {kabs:.6e} {gparam:.6e}\n")
        fh.write(f"{wl1:.2e} {wl2:.2e}\n")
        fh.write(f"{aplk[0]:.6e} {aplk[1]:.6e} {aplk[2]:.6e} {aplk[3]:.6e} {aplk[4]:.6e}\n")

    dx = 100.
    dy = 100.
    dz = 80.

    # axi = 'y'
    # level = 62
    # func = None

    # axi = 'z'
    # level = slice(0, None)
    # func = np.max
    axi = 'y'
    level = slice(0, None)
    func = np.max

    # xx, zz, ext = data_slice(data, metadata, 'extp3d', dx, dy, dz, axi, level, func)
    # xx, zz, omg = data_slice(data, metadata, 'omgp3d', dx, dy, dz, axi, level, func)
    # xx, zz, gparam = data_slice(data, metadata, 'apfp3d', dx, dy, dz, axi, level, func)
    xx, zz, ext = data_slice(data_mod, metadata_mod, 'extp3d', dx, dy, dz, axi, level, func)
    xx, zz, omg = data_slice(data_mod, metadata_mod, 'omgp3d', dx, dy, dz, axi, level, func)
    xx, zz, gparam = data_slice(data_mod, metadata_mod, 'apfp3d', dx, dy, dz, axi, level, func)

    if True:
        fig = plt.figure(figsize=(12, 3.5))
        ax1 = fig.add_subplot(1, 3, 1)
        m1 = ax1.pcolormesh(xx*1e-3, zz*1e-3, ext, cmap='Greys_r', norm=mcolors.LogNorm(vmin=np.min(ext[ext > 0]), vmax=np.max(ext)))
        ax1.set_aspect('equal')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('Z (km)')
        ax1.set_title('Extinction Coefficient')
        fig.colorbar(m1, ax=ax1, label='Extinction Coefficient (1/m)', orientation='horizontal', pad=0.2, shrink=0.65)
        ax2 = fig.add_subplot(1, 3, 2)
        m2 = ax2.pcolormesh(xx*1e-3, zz*1e-3, omg, cmap='cividis')
        ax2.set_aspect('equal')
        ax2.set_xlabel('X (km)')
        # ax2.set_ylabel('Z (km)')
        ax2.set_yticklabels([])
        ax2.set_title('Single Scattering Albedo')
        fig.colorbar(m2, ax=ax2, label='Single Scattering Albedo', orientation='horizontal', pad=0.2, shrink=0.65)
        ax3 = fig.add_subplot(1, 3, 3)
        m3 = ax3.pcolormesh(xx*1e-3, zz*1e-3, gparam, cmap='viridis', vmin=0., vmax=1.)
        ax3.set_aspect('equal')
        ax3.set_xlabel('X (km)')
        # ax3.set_ylabel('Z (km)')
        ax3.set_yticklabels([])
        ax3.set_title('Asymmetry Parameter')
        fig.colorbar(m3, ax=ax3, label='Asymmetry Parameter', orientation='horizontal', pad=0.2, shrink=0.65)
        plt.tight_layout()
        plt.savefig('out/process_atm/les_crop_%s-%s.png' % (bndstr, substr), dpi=300)