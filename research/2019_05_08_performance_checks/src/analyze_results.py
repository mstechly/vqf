import pandas as pd
import matplotlib.pyplot as plt
import pdb

def main():
    results_df = pd.read_csv("results_nm.csv", delimiter=",")
    fig_0, ax_0 = plt.subplots()
    fig_1, ax_1 = plt.subplots()
    for m in results_df.m.unique():
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
        ax_0.errorbar(m_subset.steps.unique(), m_means, yerr=m_vars, label=str(int(m)), marker="o", capsize=2)
        ax_1.errorbar(m_subset.steps.unique(), m_evals, yerr=m_evals_var, label=str(int(m)), marker="o", capsize=2)

    ax_0.set_xlabel("Number of circuit layers")
    ax_0.set_ylabel("Squared overlap")
    ax_0.set_ylim([0, 1.05])
    ax_0.legend()

    ax_1.set_xlabel("Number of circuit layers")
    ax_1.set_ylabel("BFGS function evaluations")
    ax_1.legend()


    fig_0.savefig("squared_overlap")
    fig_1.savefig("bfgs_evaluations")

if __name__ == '__main__':
    main()