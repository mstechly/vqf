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
        optimization_verbose = False
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
            pdb.set_trace()
            continue

        p_q_info = [true_p, true_q, p_dict, q_dict]
        steps = 1
        grid_size = 200
        for i in range(1):
            for gate_noise in [None]:#[1e-4, 1e-3, 1e-2]:
                for samples in [None]:#[1e2, 1e3]:
                    if samples is not None:
                        samples = int(samples)
                    print(m, gate_noise, samples, i, "     ")
                    optimization_engine = OptimizationEngine(clauses, steps=steps, grid_size=grid_size, gate_noise=gate_noise, samples=samples, verbose=optimization_verbose, visualize=True)
                    if gate_noise is None:
                        gate_noise_str = '0'
                    else:
                        gate_noise_str = str(gate_noise).split('.')[1]
                    label = "_".join([str(m), "noise", gate_noise_str, "samples", str(samples), "i", str(i)])
                    optimization_engine.simple_grid_search_angles(label=label, save_data=True)
                    del optimization_engine
                    gc.collect()
                    

if __name__ == '__main__':
    main()