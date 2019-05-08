import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import pdb


def plot_energy_landscape(x, y, values, title=None, legend_title=None, legend_min=0, legend_max=None):
    x, y = preprocess(x, y)
    fig, ax = plt.subplots()
    XX, YY = np.meshgrid(x, y)
    z = values.reshape(len(x)-1, len(y)-1)
    # ax.grid(True, which='minor', axis='both', linestyle='-', color='k')
    # ax.set_xticks(x, minor=True)
    # ax.set_yticks(y, minor=True)

    ## This is for having ticks in the plot as multiples of pi
    ax.xaxis.set_major_formatter(tck.FuncFormatter(
       lambda val,pos: '{:.2f}$\pi$'.format(val/np.pi) if val !=0 else '0'
    ))
    ax.xaxis.set_major_locator(tck.MultipleLocator(base=np.pi/4))

    ax.yaxis.set_major_formatter(tck.FuncFormatter(
       lambda val,pos: '{:.2f}$\pi$'.format(val/np.pi) if val !=0 else '0'
    ))
    ax.yaxis.set_major_locator(tck.MultipleLocator(base=np.pi/4))


    mesh_plot = ax.pcolormesh(XX, YY, z, cmap='RdBu', vmin=legend_min, vmax=legend_max)

    ax.set_xlabel("beta")
    ax.set_ylabel("gamma")
    if title is None:
        title = "QAOA energy landscape"
    ax.set_title(title)
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    cbar = fig.colorbar(mesh_plot, ax=ax)
    if legend_title is None:
        legend_title = "energy"
    cbar.set_label(legend_title)

    plt.savefig(title)
    plt.clf()

def plot_variance_landscape(betas, gammas, values):
    steps = betas.shape[1]
    grid_size = int(betas.shape[0]**(1/steps))
    for p in range(steps):
        x = betas[:, p]
        y = gammas[:, p]
        z = values.reshape([grid_size**betas.shape[1]]*2)
        mean_values = []
        var_values = []
        for beta_0 in np.unique(x):
            for gamma_0 in np.unique(y):
                betas_mask = betas[:, p] == beta_0
                gammas_mask = gammas[:, p] == gamma_0
                z_subset = z[np.ix_(betas_mask, gammas_mask)]
                mean_z = np.mean(z_subset)
                var_z = np.var(z_subset)
                mean_values.append(mean_z)
                var_values.append(var_z)


        plot_energy_landscape(np.unique(x), np.unique(y), np.array(mean_values), 
            title='steps_'+ str(steps) + ' Mean energy layer '+str(p), 
            legend_title='energy', 
            legend_min=0, 
            legend_max=np.max(values))
        plot_energy_landscape(np.unique(x), np.unique(y), np.array(var_values),
            title='steps_'+ str(steps) + ' Energy variance layer '+str(p), legend_title='variance')


def preprocess(x, y):
    x_diff = x[1] - x[0]
    y_diff = y[1] - y[0]
    x = np.append(x, x[-1] + x_diff)
    y = np.append(y, y[-1] + y_diff)
    x = x - x_diff/2
    y = y - y_diff/2
    return x, y
