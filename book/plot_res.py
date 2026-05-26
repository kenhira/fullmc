import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('font', family='Hiragino Sans')

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

if __name__ == "__main__":

    dirname = 'out/parasol'
    # dirname = 'case_test_32'
    # dirname = 'case_test_33'

    dirname_out = 'out/parasol'

    ctl_list = [
        '%s/res_b_mstrn.ctl' % dirname,
        '%s/res_b.ctl' % dirname,
        '%s/res_mca.ctl' % dirname,
                ]
    
    title_list = [
        "(a) ICA",
        "(b) 3D近似解法",
        "(c) 3Dモンテカルロ法",
                 ]
    # ctl_list = ['%s/res_b.ctl' % dirname,
    #             # '%s/res_c.ctl' % dirname,
    #             '%s/res_mca.ctl' % dirname,
    #             '%s/res_b_mstrn.ctl' % dirname,
    #             '%s/res_mcs.ctl' % dirname]
    
    # title_list = ["Parasol-10",
    #             #   "Parasol-24",
    #               "MCARaTS",
    #               "MSTRN",
    #               "MCstar"]

    typ_list = ['a1', 'a2', 'a3', 'b1']

    dx = 100.
    dy = 100.
    dz = 80.

    # axis_level_pairs = [('y', 2), ('y', 61), ('z', 0), ('z', 0), ('all', 0)]
    # plot_types = ['cross', 'cross', 'cross', 'hist', 'hist']
    axis_level_pairs = [('y', 61)]
    plot_types = ['cross']


    data_meta_list = [load_data(ctl) for ctl in ctl_list]

    data_slices = [[] for _ in data_meta_list]
    for idata, (data, metadata) in enumerate(data_meta_list):
        for ivar, varname in enumerate(typ_list):
            data_list = []
            var_index = next(i for i, var in enumerate(metadata['vars']) if var[0] == varname)
            x_count, x_start, x_step = metadata['xdef']
            y_count, y_start, y_step = metadata['ydef']
            _, z_start, z_step = metadata['zdef']
            z_count = metadata['vars'][var_index][1]
            for ial, (axis, level) in enumerate(axis_level_pairs):
                if axis == 'z':
                    x = np.linspace(x_start, x_start + (x_count - 1) * x_step * dx, x_count)
                    y = np.linspace(y_start, y_start + (y_count - 1) * y_step * dy, y_count)
                    X, Y = np.meshgrid(x, y)
                    data_slice = data[var_index, level, :, :]
                    xlabel, ylabel = "X", "Y"
                elif axis == 'y':
                    x = np.linspace(x_start, x_start + (x_count - 1) * x_step * dx, x_count)
                    z = np.linspace(z_start, z_start + (z_count - 2) * z_step * dz, z_count - 1)
                    X, Y = np.meshgrid(x, z)
                    data_slice = data[var_index, :z_count - 1, level, :]
                    xlabel, ylabel = "X", "Z"
                elif axis == 'x':
                    y = np.linspace(y_start, y_start + (y_count - 1) * y_step * dy, y_count)
                    z = np.linspace(z_start, z_start + (z_count - 2) * z_step * dz, z_count - 1)
                    X, Y = np.meshgrid(y, z)
                    data_slice = data[var_index, :z_count - 1, :, level]
                    xlabel, ylabel = "Y", "Z"
                elif axis == 'all':
                    x = np.linspace(x_start, x_start + (x_count - 1) * x_step * dx, x_count)
                    y = np.linspace(y_start, y_start + (y_count - 1) * y_step * dy, y_count)
                    z = np.linspace(z_start, z_start + (z_count - 2) * z_step * dz, z_count - 1)
                    X, Y, Z = np.meshgrid(x, y, z)
                    data_slice = data[var_index, :, :, :]
                    xlabel, ylabel = "X/Y/Z", "Value"
                
                if plot_types[ial] == 'cross':
                    data_list.append(((X, Y), data_slice))
                elif plot_types[ial] == 'hist':
                    data_list.append(data_slice)
            data_slices[idata].append(data_list)
                
    
    for ivar, varname in enumerate(typ_list):
        for ial, (axis, level) in enumerate(axis_level_pairs):
            if plot_types[ial] == 'cross':
                n = len(data_slices)
                # ncols = int(np.ceil(np.sqrt(n)))
                # nrows = int(np.ceil(n / ncols))
                ncols = 1
                nrows = n
                fig = plt.figure(figsize=(4.8 * ncols, 2 * nrows))
                # vmin = min([np.nanmin(data_slices[idata][ivar][ial][1]) for idata in range(len(data_slices)) if np.isfinite(data_slices[idata][ivar][ial][1]).any()])
                vmin = 0.0
                vmax = max([np.nanmax(data_slices[idata][ivar][ial][1]) for idata in range(len(data_slices)) if np.isfinite(data_slices[idata][ivar][ial][1]).any()])
                cmap = 'Reds' if varname == 'b1' else 'Blues_r'
                if axis == 'z':
                    xlabel, ylabel = "X (km)", "Y (km)"
                elif axis == 'y':
                    xlabel, ylabel = "X (km)", "Z (km)"
                elif axis == 'x':
                    xlabel, ylabel = "Y (km)", "Z (km)"
                for idata, data_slice in enumerate(data_slices):
                    X, Y = data_slice[ivar][ial][0]
                    dat = data_slice[ivar][ial][1]
                    ax = fig.add_subplot(nrows, ncols, idata + 1)
                    title = title_list[idata] if idata < len(title_list) else os.path.basename(ctl_list[idata])
                    contour = ax.pcolormesh(X*1e-3, Y*1e-3, dat, cmap=cmap, vmin=vmin, vmax=vmax)
                    ax.text(0.02, 0.95, f"{title}", transform=ax.transAxes, 
                            va='top', ha='left', fontsize=10,)
                    ax.set_aspect('equal')
                    if idata // ncols == nrows - 1:
                        ax.set_xlabel(xlabel)
                    else:
                        ax.set_xticklabels([])
                    ax.set_ylabel(ylabel)
                cb = fig.colorbar(contour, ax=fig.axes, orientation='horizontal', pad=0.12, shrink=0.6)
                # cb.set_label("Irradiance (W/m^2)" if varname.startswith('a') else "Flux convergence (W/m^3)")
                cb.set_label("放射フラックス (W/m²)" if varname.startswith('a') else "(W/m³)")
                # plt.tight_layout()
                fig.savefig("%s/cross_%s_%s_%03d.png" % (dirname_out, varname, axis, level), bbox_inches='tight', dpi=450)
            elif plot_types[ial] == 'hist':
                vmin = min([np.nanmin(data_slices[idata][ivar][ial]) for idata in range(len(data_slices)) if np.isfinite(data_slices[idata][ivar][ial]).any()])
                vmax = max([np.nanmax(data_slices[idata][ivar][ial]) for idata in range(len(data_slices)) if np.isfinite(data_slices[idata][ivar][ial]).any()])
                bins = np.linspace(vmin, vmax, 50)
                
                col_list = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']
                fig = plt.figure(figsize=(8, 7))
                ax1 = fig.add_subplot(211)
                for idata, data_slice in enumerate(data_slices[:2]):
                    dat = data_slice[ivar][ial]
                    ax1.hist(dat.flatten(), bins=bins, density=True, histtype='stepfilled',
                             linewidth=1.5, label=title_list[idata] if idata < len(title_list) else os.path.basename(ctl_list[idata]),
                             color=col_list[idata % len(col_list)], alpha=0.5)
                ax1.set_xlim(vmin, vmax)
                ax1.set_xticklabels([])
                ax1.set_title(f"Density Histogram of {varname} ({axis}={level})")
                ax1.set_ylabel("Density")
                ax1.legend()
                ax2 = fig.add_subplot(212)
                for idata, data_slice in enumerate(data_slices[2:]):
                    dat = data_slice[ivar][ial]
                    ax2.hist(dat.flatten(), bins=bins, density=True, histtype='stepfilled',
                             linewidth=1.5, label=title_list[idata + 2] if (idata + 2) < len(title_list) else os.path.basename(ctl_list[idata + 2]),
                             color=col_list[(idata + 2) % len(col_list)], alpha=0.5)
                ax2.set_xlim(vmin, vmax)
                ax2.set_xlabel("Irradiance (W/m^2)" if varname.startswith('a') else "Flux convergence (W/m^3)")
                ax2.set_ylabel("Density")
                ax2.legend()
                plt.tight_layout()
                fig.savefig("%s/hist_%s_%s_%03d.png" % (dirname_out, varname, axis, level), dpi=300)

