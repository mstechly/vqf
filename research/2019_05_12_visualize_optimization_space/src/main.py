from vqf.preprocessing import create_clauses, calculate_number_of_unknowns
from vqf.preprocessing import factor_56153, factor_291311
from vqf.optimization import OptimizationEngine
from sympy import Add, Mul, Symbol
import numpy as np
import pdb
import gc

def main():
    p_q_m_list = [[283, 7, 1981], [263, 11, 2893], [241, 233, 56153], [557, 523, 291311], [29, 11, 319], [263, 263, 69169]]
    results = []
    unknowns_list = [[2, 0], [3, 1], [4, 0], [6, 0], [6, 3], [8, 5]]
    for p_q_m, unknowns in zip(p_q_m_list, unknowns_list):
        true_p = p_q_m[0]
        true_q = p_q_m[1]
        m = p_q_m[2]
        apply_preprocessing = True
        preprocessing_verbose = False
        optimization_verbose = True
        number_of_unknowns = 0
        carry_bits = 0
        counter = 0
        if m == 56153:
            p_dict, q_dict, z_dict, clauses = factor_56153()
        elif m == 291311:
            p_dict, q_dict, z_dict, clauses = factor_291311()
        else:
            p_dict, q_dict, z_dict, clauses = create_clauses(m, true_p, true_q, apply_preprocessing, preprocessing_verbose)
        number_of_unknowns, carry_bits = calculate_number_of_unknowns(p_dict, q_dict, z_dict)

        if number_of_unknowns != unknowns[0] or carry_bits != unknowns[1]:
            print("Got wrong number of unknowns!")
            continue

        p_q_info = [true_p, true_q, p_dict, q_dict]
        step_by_step_results = None
        for steps in range(1, 9):
            for i in range(3):
                print(m, steps, i)                
                optimization_engine = OptimizationEngine(clauses, steps=steps, grid_size=grid_size, tol=1e-10, gate_noise=1e-3, verbose=optimization_verbose, visualize=False)
                optimization_engine.step_by_step_results = step_by_step_results
                squared_overlap, bfgs_evaluations, step_by_step_results = run_single_case(p_q_info, optimization_engine)
                
                print(squared_overlap)
                results.append([m, steps, squared_overlap, bfgs_evaluations])

                # optimization_history = optimization_engine.optimization_history
                # history_file_name = "_".join([str(m), str(steps), str(i), "history"]) + ".csv"
                # np.savetxt(history_file_name, optimization_history, delimiter=",")
            np.savetxt("results.csv", results, delimiter=",", header="m,steps,squared_overlap,bfgs_evaluations", fmt='%.4f', comments='')

if __name__ == '__main__':
    main()