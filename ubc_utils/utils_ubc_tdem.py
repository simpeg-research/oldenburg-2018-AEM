import numpy as np
from SimPEG import Utils, Mesh
from scipy.interpolate import NearestNDInterpolator


def generate_trapezoidal_waveform(
    ramp_off_time=3e-6,
    on_time=0.02,
    n_sample_to_peak_time=5,
    n_sample_in_ramp_off=3,
    n_segment = 5,
    n_samples_per_segment = 5,
    t_min = 1e-5,
    t_max = 1e-2,
    time_steps_off_time = None,
    fname_wave='wave.txt'
):
    dt_on = (0.02-ramp_off_time) / (n_sample_to_peak_time)
    dt_ramp_off = ramp_off_time / n_sample_in_ramp_off

    # TODO: make it bit more general, so use n_segment, n_samples_per_segment
    # t_min and t_max
    if time_steps_off_time is None:
        raise Exception()

    ht = [(dt_on, n_sample_to_peak_time), (dt_ramp_off, n_sample_in_ramp_off)] + time_steps_off_time
    mesh_time = Mesh.TensorMesh([ht], x0=[-on_time])
    time_current = mesh_time.vectorNx
    time_current[abs(time_current)<1e-10] = 0.
    current = trapezoidal_waveform(time_current, on_time, ramp_off_time)
    np.savetxt(fname_wave, np.c_[time_current, current])
    return time_current, current


def trapezoidal_waveform(time_current, on_time, ramp_off_time):
    current = np.zeros_like(time_current)
    # immediate ramp on
    current[time_current<=-ramp_off_time] = 1.
    current[0] = 0.
    t_inds = np.logical_and(time_current>-ramp_off_time, time_current<=0.)
    current[t_inds] = -1./ramp_off_time * time_current[t_inds]
    return current


def generate_airborne_tdem_data(
    nx, ny, dx, dy,
    time,
    x0=0., y0=0.,
    topo=None,
    loop_radius = 13.,
    src_height = 30.,
    rx_height = 30.,
    src_rx_offset = 0.,
    src_type = 'TRX_LOOP',
    fname_data="data_locations.txt",
    fname_src_locations="src_points.txt",
    fname_rx_locations="rx_points.txt",
    fname_topo = 'topography.txt',
    data_ignore = '-0',
    alpha=0.,
    theta=0.
):
    x = x0 + np.arange(nx)*dx
    y = y0 + np.arange(ny)*dy
    n_src = nx*ny
    n_time = time.size

    # Update this for x and z
    offset = np.ones((nx, ny)) * src_rx_offset
    offset[::2, :] = -src_rx_offset
    offset = Utils.mkvc(offset)

    if topo is None:
        z = 0.
        xyz = Utils.ndgrid(x, y, np.r_[z+src_height])
        xyz_rx = np.c_[xyz[:,0]-offset, xyz[:,1], xyz[:,2]]
        xyz_rx = np.c_[xyz[:,0]-offset, xyz[:,1], xyz[:,2]*z+rx_height]
    else:
        F = NearestNDInterpolator(topo[:,:2], topo[:,2])
        xyz = Utils.ndgrid(x, y, np.r_[0.])
        z = F(xyz[:,:2])
        xyz[:,2] = z + src_height
        xyz_rx = np.c_[xyz[:,0]-offset, xyz[:,1], z+rx_height]

    with open(fname_data, 'w') as file:
        file.write("IGNORE " + data_ignore + " \n")
        file.write("N_TRX " + str(n_src) + " \n")
        file.write("\n")
        dummy = " -0 " * 16 + " 1 1 "
        for i_src in range(n_src):
            file.write(src_type+"\n")
            file.write(
              ("   %.5e %.5e %.5e %.5e %.5e %.5e \n") %
              (xyz[i_src,0], xyz[i_src,1], xyz[i_src,2], loop_radius, alpha, theta)
            )
            file.write("N_RECV 1\n")
            file.write("N_TIME "+str(n_time)+" \n")
            for t in time:
                file.write(
                  ("   %.5e %.5e %.5e %.5e " + dummy+ " \n") %
                  (xyz[i_src,0]-offset[i_src], xyz[i_src,1], xyz[i_src,2], t)
                )
            file.write("\n")
            file.write("\n")
            file.write("\n")
        np.savetxt(fname_src_locations, xyz)
        np.savetxt(fname_rx_locations, xyz_rx)
        if topo is not None:
            with open(fname_topo, 'w') as file:
                file.write(("%i \n") % (int(topo.shape[0])))
                for i in range(topo.shape[0]):
                    file.write(
                        (" %.5e %.5e %.5e \n") %
                        (topo[i, 0], topo[i, 1], topo[i, 2])
                    )
    output = {
        "src_locs": xyz,
        "rx_locs": xyz_rx,
        "topography": topo,
        "time": time
    }
    return output


def read_tdoctree_dpred(fname):
    with open(fname) as file:

        lines = file.readlines()

        src = lines[1].split()
        if src[0] == 'N_TRX':
            n_src = int(src[1])
        else:
            print (src)
            raise Exception("Wrong format")

        src_type = lines[3].split()[0]
        i_count = 3
        data_type = []
        data = []
        time = []
        data_header = [
                        "Ex",
                        "Ey",
                        "Ez",
                        "Hx",
                        "Hy",
                        "Hz",
                        "dBxdt",
                        "dBydt",
                        "dBzdt"
        ]

        for line in lines[2:]:
            if len(line.split()) == 13:
                data.append(np.array(line.split(), dtype=float))
        data = np.vstack(data)
        output = {
            'data_header': data_header,
            'data': data[:, 4:],
            'xyz_data': data[:, :3],
        }
    return output


def read_tdoctree_dpred_inversion(fname):
    data = np.loadtxt(fname, comments='%')
    data_header = [
                    "Ex",
                    "Ey",
                    "Ez",
                    "Hx",
                    "Hy",
                    "Hz",
                    "dBxdt",
                    "dBydt",
                    "dBzdt"
    ]
    output = {
        'data_header': data_header,
        'data': data[:, 3:],
        'xyz_data': data[:, :3],
    }
    return output


def read_tdoctree_dobs(fname):
    with open(fname) as file:

        lines = file.readlines()

        src = lines[1].split()
        if src[0] == 'N_TRX':
            n_src = int(src[1])
        else:
            print (src)
            raise Exception("Wrong format")

        src_type = lines[3].split()[0]
        i_count = 3
        data_type = []
        data = []
        time = []
        data_header = [
                        "Ex",
                        "Ey",
                        "Ez",
                        "Hx",
                        "Hy",
                        "Hz",
                        "dBxdt",
                        "dBydt",
                        "dBzdt"
        ]
        uncertainty_header = [
                        "Ex_uncertainty",
                        "Ey_uncertainty",
                        "Ez_uncertainty",
                        "Hx_uncertainty",
                        "Hy_uncertainty",
                        "Hz_uncertainty",
                        "dBxdt_uncertainty",
                        "dBydt_uncertainty",
                        "dBzdt_uncertainty"
        ]
        for line in lines[2:]:
            if len(line.split()) == 22:
                data.append(np.array(line.split(), dtype=float))
        data = np.vstack(data)
        data_inds = 4 + np.arange(9)*2
        uncertainty_inds = 5 + np.arange(9)*2
        output = {
            'data_header': data_header,
            'uncertainty_header': uncertainty_header,
            'data': data[:, data_inds],
            'uncenrtainty': data[:, uncertainty_inds],
            'xyz_data': data[:, :3],
            'time_data': data[:, 3]
        }
    return output


def write_tdoctree_dobs(
    src_locations,
    rx_locations,
    time,
    data,
    uncenrtainty,
    src_type='TRX_LOOP',
    fname_dobs="dobs.txt",
    data_ignore='-0',
    loop_radius=13.,
    alpha=0.,
    theta=0.,
    dtype='dBzdt'
):

    n_src = src_locations.shape[0]
    n_time = time.size

    with open(fname_dobs, 'w') as file:
        file.write("IGNORE " + data_ignore + " \n")
        file.write("N_TRX " + str(n_src) + " \n")
        file.write("\n")
        if dtype == 'dBzdt':
            dummy = " -0 " * 16
            print_string = "   %.5e %.5e %.5e %.5e " + dummy + " %.5e %.5e \n"
        elif dtype =='Hz':
            dummy_1 = " -0 " * 10
            dummy_2 = " -0 " * 6
            print_string = "   %.5e %.5e %.5e %.5e " + dummy_1 + " %.5e %.5e" + dummy_2 + "\n"
        
        i_count = 0
        for i_src in range(n_src):
            file.write(src_type+"\n")
            file.write(
              ("   %.5e %.5e %.5e %.5e %.5e %.5e \n") %
              (
                src_locations[i_src, 0],
                src_locations[i_src, 1],
                src_locations[i_src, 2],
                loop_radius, alpha, theta
              )
            )
            file.write("N_RECV 1\n")
            file.write("N_TIME "+str(n_time)+" \n")
            for t in time:
                file.write(
                  (print_string) %
                  (
                    rx_locations[i_src, 0],
                    rx_locations[i_src, 1],
                    rx_locations[i_src, 2],
                    t, data[i_count], uncenrtainty[i_count]
                  )
                )
                i_count += 1
            file.write("\n")
            file.write("\n")
            file.write("\n")
