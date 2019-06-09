import pandas as pd
import matplotlib.pyplot as plt
import os
import pdb

def main():
    data_path = os.path.join("..", "results")
    bfgs_basic_df = pd.read_csv(os.path.join(data_path, "results_bfgs.csv"), delimiter=',')
    bfgs_interp_df = pd.read_csv(os.path.join(data_path, "results_interp_bfgs.csv"), delimiter=',')
    lbfgsb_basic_df = pd.read_csv(os.path.join(data_path, "results_lbfgsb.csv"), delimiter=',')
    lbfgsb_interp_df = pd.read_csv(os.path.join(data_path, "results_interp_lbfgsb.csv"), delimiter=',')
    nm_basic_df = pd.read_csv(os.path.join(data_path, "results_nm.csv"), delimiter=',')
    nm_interp_df = pd.read_csv(os.path.join(data_path, "results_interp_nm.csv"), delimiter=',')

    all_data = [(bfgs_basic_df, bfgs_interp_df, 'bfgs'), (lbfgsb_basic_df, lbfgsb_interp_df, 'lbfgsb'), (nm_basic_df, nm_interp_df, 'nm')]

    for basic_df, interp_df, name in all_data:
        fig_0, ax_0 = plt.subplots()
        fig_1, ax_1 = plt.subplots()
        linestyle_basic = '-'
        linestyle_interp = '--'
        ax_0, ax_1 = create_plots(ax_0, ax_1, basic_df, linestyle_basic)
        ax_0, ax_1 = create_plots(ax_0, ax_1, interp_df, linestyle_interp)

        ax_0.set_xlabel("Number of circuit layers")
        ax_0.set_ylabel("Squared overlap")
        ax_0.set_ylim([0, 1.05])
        ax_0.legend()

        ax_1.set_xlabel("Number of circuit layers")
        ax_1.set_ylabel("BFGS function evaluations")
        ax_1.legend()

        fig_0.savefig("squared_overlap_" + name)
        fig_1.savefig("evaluations_" + name)

def create_plots(ax_0, ax_1, results_df, linestyle):
    color_dict = {1981: 'tab:blue', 319: 'tab:orange', 69169: 'tab:green', 2893: 'tab:red', 56153: 'tab:purple', 291311: 'tab:brown'}
    for m in results_df.m.unique():
        color = color_dict[m]
        m_subset = results_df[results_df.m == m]
        m_means = []
        m_vars = []
        m_evals = []
        m_evals_var = []
        for steps in m_subset.steps.unique():
            steps_subset = m_subset[m_subset.steps == steps]
            mean_overlap = steps_subset.squared_overlap.mean()
            var_overlap = steps_subset.squared_overlap.var()
            m_means.append(mean_overlap)
            m_vars.append(var_overlap)

            mean_evals = steps_subset.bfgs_evaluations.mean()
            var_evals = steps_subset.bfgs_evaluations.var()
            m_evals.append(mean_evals)
            m_evals_var.append(var_evals)
        if linestyle == '-':
            label = str(int(m))
        else:
            label = None
        ax_0.errorbar(m_subset.steps.unique(), m_means, yerr=m_vars, label=label, marker="o", linestyle=linestyle, capsize=2, color=color)
        ax_1.errorbar(m_subset.steps.unique(), m_evals, yerr=m_evals_var, label=label, marker="o", linestyle=linestyle, capsize=2, color=color)
    return ax_0, ax_1

if __name__ == '__main__':
    main()