import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backend_bases import MouseButton
import matlab.engine as me
from time import perf_counter_ns
from scipy.interpolate import griddata
import os

data_dir = "2d_cache"
query_matlab = False
point_size = 40
selected_point_size = 100
default_font_size = 18 #14
project3d = False #True
mpl.rcParams['font.size'] = default_font_size
mpl.rcParams['axes.labelsize'] = default_font_size
mpl.rcParams['axes.titlesize'] = default_font_size
mpl.rcParams['xtick.labelsize'] = default_font_size - 4
mpl.rcParams['ytick.labelsize'] = default_font_size

def prior_sampling_2d():
    engine = me.start_matlab()
    engine.cd("~/Documents/sabl/examples/Adv_Diff_2D_Discrepancy_Calibration", nargout=0)
    result = engine.Visualization_Prior_Sampling(nargout=7)
    engine.quit()
    results = map(np.array, result)
    return results

if query_matlab:
    os.makedirs(f"{data_dir}", exist_ok=True)
    t_query_start = perf_counter_ns()
    x,y, delta_samples_z_opt, delta_samples_z_pert, z_pert, corlengths, magnitudes = prior_sampling_2d()
    t_query_end = perf_counter_ns()
    query_time = np.round((t_query_end - t_query_start) / 1e9,3)
    print(f"Matlab query time: {query_time//60} minutes {query_time%60} seconds")
    np.save(f"{data_dir}/x.npy", x)
    np.save(f"{data_dir}/y.npy", y)
    np.save(f"{data_dir}/delta_samples_z_opt.npy", delta_samples_z_opt)
    np.save(f"{data_dir}/delta_samples_z_pert.npy", delta_samples_z_pert)
    np.save(f"{data_dir}/z_pert.npy", z_pert)
    np.save(f"{data_dir}/corlengths.npy", corlengths)
    np.save(f"{data_dir}/magnitudes.npy", magnitudes)
else:
    x,y,delta_samples_z_opt,delta_samples_z_pert,z_pert,corlengths,magnitudes = [np.load(f"{data_dir}/{name}.npy") for name in ["x","y","delta_samples_z_opt","delta_samples_z_pert","z_pert","corlengths","magnitudes"]]

# set the resolution of 2D domain evaluation
x_mesh,y_mesh = np.meshgrid(np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100))
# grid coordinates of the scalar field
grid_points = np.concat([x,y], axis=1)
n_perts, n_samples = delta_samples_z_pert.shape[:2]
print(f'{n_perts=}, {n_samples=}')
# use linear interpolation to resample the data to the 2D grid
def resample_to_image(data):
    return griddata(grid_points, data, (x_mesh, y_mesh), method='linear')

z_pert_indices = np.arange(n_perts) # tracking the perturbation indices
corlengths = corlengths.squeeze()
sorted_indices = np.argsort(corlengths)

sorted_corlength = corlengths[sorted_indices]
unique_corlength, unique_counts = np.unique(sorted_corlength, return_counts=True)
sorted_delta_samples_z_pert = delta_samples_z_pert[sorted_indices]
sorted_magnitudes = magnitudes[:,sorted_indices]
sorted_z_pert = z_pert[:,sorted_indices]

sorted_points = np.stack([np.broadcast_to(sorted_corlength,sorted_magnitudes.shape), sorted_magnitudes], axis=2)

rescaled_points = sorted_points.copy()
x_range = sorted_points[:,:,0].max() - sorted_points[:,:,0].min()
y_range = sorted_points[:,:,1].max() - sorted_points[:,:,1].min()
rescaled_points[:,:,0] /= x_range
rescaled_points[:,:,1] /= y_range
print(sorted_points.shape)

selected_index = (0,0)
clicked_point = None

fig = plt.figure()
gs_root = fig.add_gridspec(6,28)

scatter_axis = fig.add_subplot(gs_root[:3,:])
delta_z_opt_axis = fig.add_subplot(gs_root[4:,0:6],projection='3d') if project3d else fig.add_subplot(gs_root[4:,0:6])
delta_z_opt_cb_axis = fig.add_subplot(gs_root[4:,6:7])
delta_z_pert_axis = fig.add_subplot(gs_root[4:,10:16],projection='3d') if project3d else fig.add_subplot(gs_root[4:,10:16])
delta_z_pert_cb_axis = fig.add_subplot(gs_root[4:,16:17])
z_pert_axis = fig.add_subplot(gs_root[4:,20:26],projection='3d') if project3d else fig.add_subplot(gs_root[4:,20:26])
z_pert_cb_axis = fig.add_subplot(gs_root[4:,26:27])

vmin = magnitudes.min()
vmax = magnitudes.max()

corlength_diff = corlengths[1:] - corlengths[:-1]

sc_colors = []

def plot_images():
    theta_idx = selected_index[0]
    pert_idx = selected_index[1]

    delta_z_opt_image = resample_to_image(delta_samples_z_opt[:,theta_idx])
    actor = delta_z_opt_axis.plot_surface(x_mesh, y_mesh, delta_z_opt_image, cmap='viridis', vmin=delta_samples_z_opt.min(), vmax=delta_samples_z_opt.max(), edgecolor='none') if project3d else delta_z_opt_axis.imshow(delta_z_opt_image, cmap='viridis', vmin=delta_samples_z_opt.min(), vmax=delta_samples_z_opt.max(), origin='lower')
    delta_z_opt_axis.set_title(f"Discrepancy sample {theta_idx} at z_opt")
    fig.colorbar(actor, cax=delta_z_opt_cb_axis)
    delta_z_opt_axis.set_xticks([])
    delta_z_opt_axis.set_yticks([])
    
    
    delta_z_pert_image = resample_to_image(sorted_delta_samples_z_pert[pert_idx,:,theta_idx])
    #actor = delta_z_pert_axis.plot_surface(x_mesh, y_mesh, delta_z_pert_image, cmap='viridis', vmin=delta_samples_z_pert.min(), vmax=delta_samples_z_pert.max(), edgecolor='none') if project3d else delta_z_pert_axis.imshow(delta_z_pert_image, cmap='viridis', vmin=delta_samples_z_pert.min(), vmax=delta_samples_z_pert.max(), origin='lower')
    actor = delta_z_pert_axis.plot_surface(x_mesh, y_mesh, delta_z_pert_image, cmap='viridis', edgecolor='none') if project3d else delta_z_pert_axis.imshow(delta_z_pert_image, cmap='viridis', origin='lower')
    delta_z_pert_axis.set_title(f"Discrepancy sample {theta_idx} at z_{z_pert.shape[1]-pert_idx}")
    fig.colorbar(actor, cax=delta_z_pert_cb_axis)
    delta_z_pert_axis.set_xticks([])
    delta_z_pert_axis.set_yticks([])
    
    z_pert_image = resample_to_image(sorted_z_pert[:,pert_idx])
    actor = z_pert_axis.plot_surface(x_mesh, y_mesh, z_pert_image, cmap='viridis', vmin=z_pert.min(), vmax=z_pert.max(), edgecolor='none') if project3d else z_pert_axis.imshow(z_pert_image, cmap='viridis', vmin=z_pert.min(), vmax=z_pert.max(), origin='lower')
    z_pert_axis.set_title(f"Perturbed z {z_pert.shape[1]-pert_idx}")
    fig.colorbar(actor, cax=z_pert_cb_axis)
    z_pert_axis.set_xticks([])
    z_pert_axis.set_yticks([])

    if project3d:
        delta_z_opt_axis.set_zlim(delta_samples_z_opt.min(),delta_samples_z_opt.max())
        delta_z_pert_axis.set_zlim(delta_samples_z_pert.min(),delta_samples_z_pert.max())
        z_pert_axis.set_zlim(z_pert.min(),z_pert.max())

def plot_scatter(plot_selection=False):
    scatter_axis.clear()
    delta_z_opt_axis.clear()
    delta_z_opt_cb_axis.clear()
    delta_z_pert_axis.clear()
    delta_z_pert_cb_axis.clear()
    z_pert_axis.clear()
    z_pert_cb_axis.clear()
    prev_corlength = sorted_corlength[0]
    _violin_data = []
    _group_data = []
    for i in range(n_perts):
        if sorted_corlength[i] != prev_corlength:
            _violin_data.append(np.array(_group_data).flatten())
            _group_data = []
        _group_data.append(sorted_magnitudes[:,i])
        prev_corlength = sorted_corlength[i]
    _violin_data.append(np.array(_group_data).flatten())
    parts = scatter_axis.violinplot(_violin_data, positions=unique_corlength, showmeans=False, showmedians=False, showextrema=False,widths=np.abs(corlength_diff[corlength_diff.nonzero()]).min()/2)

    
    for pc in parts['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('black')
        pc.set_alpha(1)
        pc.set_linewidth(2)

    for i in range(n_perts):
        _scatter_x = np.zeros_like(sorted_magnitudes[:,i]) + sorted_corlength[i]
        sc = scatter_axis.scatter(_scatter_x, sorted_magnitudes[:,i], s=point_size, alpha=1, marker='o')
        if len(sc_colors) <= n_perts:
            sc_colors.append(sc.get_facecolor())

    scatter_axis.set_xticks(unique_corlength)
    labels = [f'{unique_corlength[i]:.2f}({unique_counts[i]})' for i in range(len(np.unique(sorted_corlength)))]
    scatter_axis.set_xticklabels(labels)
    
    
    scatter_axis.set_xlabel("correlation lengths")
    scatter_axis.set_ylabel("magnitude")
    if plot_selection:
        _x = selected_index[0]
        _y = selected_index[1]

        scatter_axis.scatter(sorted_points[_x,_y,0], sorted_points[_x,_y,1], s=selected_point_size, alpha=1, marker='o', color=sc_colors[_y], edgecolors='black', linewidths=2)
        #scatter_axis.axhline(y=sorted_points[_x,_y,1], color='black', linestyle='--')
        plot_images()
    fig.canvas.draw()


def on_click(event):
    global selected_index
    if event.inaxes == scatter_axis:
        clicked_point = (event.xdata/x_range, event.ydata/y_range)
        distances = np.linalg.norm(rescaled_points - clicked_point, axis=2)
        selected_index = np.unravel_index(np.argmin(distances), distances.shape)
        plot_scatter(True)
    if event.button == MouseButton.RIGHT:
        selected_index = (0,0)
        plot_scatter(False)




def on_move(event):
    if not project3d:
        return
    if event.inaxes == delta_z_opt_axis:
        in_ax = delta_z_opt_axis
        other_ax = (delta_z_pert_axis,z_pert_axis)
    elif event.inaxes == delta_z_pert_axis:
        in_ax = delta_z_pert_axis
        other_ax = (delta_z_opt_axis,z_pert_axis)
    elif event.inaxes == z_pert_axis:
        in_ax = z_pert_axis
        other_ax = (delta_z_opt_axis,delta_z_pert_axis)
    else:
        return
    for ax in other_ax:
        ax.view_init(in_ax.elev, in_ax.azim)
    fig.canvas.draw_idle()

plot_scatter(False)
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('motion_notify_event', on_move)
plt.show()