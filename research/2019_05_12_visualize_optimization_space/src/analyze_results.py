import os, glob
import pdb
from vqf.visualization import plot_energy_landscape
import numpy as np


def recreate_plots_for_given_m(current_dir):
    data_file_paths = glob.glob(os.path.join(current_dir, "*.csv"))
    list_of_data = []
    for data_file_path in data_file_paths:
        list_of_data.append(np.genfromtxt(data_file_path, delimiter=","))

    all_data = np.hstack(list_of_data)
    max_value = np.max(all_data[:,2])
    min_value = np.min(all_data[:,2])
    m = current_dir.split("_")[0]
    plots_path = m + "_plots"
    log_plots_path = m + "_plots_log"

    os.mkdir(plots_path)
    os.mkdir(log_plots_path)
    for data, file_path in zip(list_of_data, data_file_paths):
        print("Plotting", file_path, end="\r")
        grid_size = int(np.sqrt(len(data)))
        x = data[:, 0][::grid_size]
        y = data[:, 1][:grid_size]
        values = data[:, 2]
        file_name = file_path.split("/")[1][:-4]
        plot_energy_landscape(x, y, values, title=file_name, dir_path=plots_path, legend_min=min_value, legend_max=max_value)
        plot_energy_landscape(x, y, values, title=file_name+"_log", dir_path=log_plots_path, log_legend=True, legend_min=min_value, legend_max=max_value)


def main():
    result_dirs = glob.glob('*_results')
    for current_dir in result_dirs:
        recreate_plots_for_given_m(current_dir)
    

if __name__ == '__main__':
    main()