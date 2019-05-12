from preprocessing import create_clauses, calculate_number_of_unknowns
from preprocessing import factor_56153, factor_291311
from optimization import OptimizationEngine
from sympy import Add, Mul, Symbol
import pdb


def factor_number(m, true_p, true_q, use_true_values=False):
    apply_preprocessing = True
    preprocessing_verbose = False
    optimization_verbose = False
    if m == 56153:
        p_dict, q_dict, z_dict, clauses = factor_56153()
    elif m == 291311:
        p_dict, q_dict, z_dict, clauses = factor_291311()
    elif use_true_values:
        p_dict, q_dict, z_dict, clauses = create_clauses(m, true_p, true_q, apply_preprocessing, preprocessing_verbose)
    else:
        p_dict, q_dict, z_dict, clauses = create_clauses(m, None, None, apply_preprocessing, preprocessing_verbose)

    number_of_uknowns, number_of_carry_bits = calculate_number_of_unknowns(p_dict, q_dict, z_dict)
    print("Number of unknowns:", number_of_uknowns)
    print("Number of carry bits:", number_of_carry_bits)
    if clauses[0] == 0 and len(set(clauses)) == 1:
        if number_of_uknowns == 0:
            p, q = decode_solution(p_dict, q_dict)
            return p, q, None

    optimization_engine = OptimizationEngine(clauses, m, steps=1, grid_size=20, gate_noise=None, verbose=optimization_verbose, visualize=True)
    sampling_results, mapping = optimization_engine.perform_qaoa()
    most_frequent_bit_string = max(sampling_results, key=lambda x: sampling_results[x])
    
    squared_overlap = calculate_squared_overlap(mapping, sampling_results, true_p, true_q, p_dict, q_dict)
    p_dict = update_dictionary(most_frequent_bit_string, mapping, p_dict)
    q_dict = update_dictionary(most_frequent_bit_string, mapping, q_dict)
    z_dict = update_dictionary(most_frequent_bit_string, mapping, z_dict)

    p, q = decode_solution(p_dict, q_dict)
    return p, q, squared_overlap


def calculate_squared_overlap(mapping, sampling_results, true_p, true_q, p_dict, q_dict):
    p_binary_string = bin(true_p)[2:][::-1]
    q_binary_string = bin(true_q)[2:][::-1]

    p_binary = [int(char) for char in p_binary_string]
    q_binary = [int(char) for char in q_binary_string]
    if len(p_binary) < len(p_dict):
        trailing_zeros = len(p_dict) - len(p_binary)
        for zero in range(trailing_zeros):
            p_binary.append(0)

    if len(q_binary) < len(q_dict):
        trailing_zeros = len(q_dict) - len(q_binary)
        for zero in range(trailing_zeros):
            q_binary.append(0)

    all_correct_assignments = []
    correct_assignment = {}
    for q_id, q_val in q_dict.items():
        if type(q_val) is Symbol:
            bit_id = mapping[str(q_val)]
            correct_value = q_binary[q_id]
            if bit_id not in correct_assignment.keys():
                correct_assignment[bit_id] = correct_value

    for p_id, p_val in p_dict.items():
        if type(p_val) is Symbol:
            bit_id = mapping[str(p_val)]
            correct_value = p_binary[p_id]
            if bit_id not in correct_assignment.keys():
                correct_assignment[bit_id] = correct_value

    all_correct_assignments.append(correct_assignment)
    # TODO: 
    # This is just a hack for 56153 and 291311 to work properly.
    # It should be generalized to work for any symmetric case.

    if (true_p == 241 and true_q == 233):
        assignment_1 = {mapping['p_3']: 0, mapping['p_4']: 1, mapping['q_3']: 1, mapping['q_4']: 0}
        assignment_2 = {mapping['p_3']: 1, mapping['p_4']: 0, mapping['q_3']: 0, mapping['q_4']: 1}
        all_correct_assignments = [assignment_1, assignment_2]

    if (true_p == 557 and true_q == 523):
        assignment_1 = {mapping['p_1']: 0, mapping['p_2']: 1, mapping['p_5']: 1, mapping['q_1']: 1, mapping['q_2']: 0, mapping['q_5']: 0}
        assignment_2 = {mapping['p_1']: 1, mapping['p_2']: 0, mapping['p_5']: 0, mapping['q_1']: 0, mapping['q_2']: 1, mapping['q_5']: 1}
        all_correct_assignments = [assignment_1, assignment_2]


    total_overlap = 0
    total_count = 0
    print(all_correct_assignments)
    print(mapping)
    squared_overlap = 0
    for correct_assignment in all_correct_assignments:
        for bit_string, count in sampling_results.most_common():
            correct_count = 0
            for bit_id, bit_value in enumerate(bit_string):
                # This accounts for the fact some of the bits of the sampling results
                # are irrelevant to the result - namely, carry bits.
                if bit_id not in correct_assignment.keys():
                    continue
                if bit_value == correct_assignment[bit_id]:
                    correct_count += 1
            overlap = (correct_count / len(correct_assignment))**2 * count
            total_count += count
            print(bit_string, count, correct_count, overlap)
            total_overlap += overlap
        print("_"*10)
        total_overlap = total_overlap / total_count
        squared_overlap += total_overlap
    return squared_overlap


def update_dictionary(qaoa_solution, mapping, x_dict):
    # values_dict = {symbol_str: qaoa_solution[index] for symbol_str, index in mapping.items()}
    symbols_dict = {Symbol(symbol_str): qaoa_solution[index] for symbol_str, index in mapping.items()}
    for key, value in x_dict.items():
        if value in symbols_dict.keys():
            x_dict[key] = symbols_dict[value]
        if type(value) == Add or type(value) == Mul:
            x_dict[key] = value.subs(symbols_dict)
    return x_dict


def decode_solution(p_dict, q_dict):
    p = 0
    for key, value in p_dict.items():
        p += value * 2**key

    q = 0
    for key, value in q_dict.items():
        q += value * 2**key

    return p, q


def main():
    p_q_m_list = [[283, 7, 1981], [29, 11, 319], [263, 263, 69169], [263, 11, 2893], [241, 233, 56153], [557, 523, 291311]]
    for p_q_m in p_q_m_list:
        true_p = p_q_m[0]
        true_q = p_q_m[1]
        m = p_q_m[2]
        use_true_values = True

        print("M:", m)
        if m % 2 == 0:
            p = 2
            q = int(m / 2)
            print("The primes are:", p, "and", q)
        else:
            p, q, squared_overlap = factor_number(m, true_p, true_q, use_true_values=use_true_values)
            print("Calculated primes of ",m, "are:", p, "and", q)
            print("      True primes of ",m, "are:", true_p, "and", true_q)
            print("Squared overlap:", squared_overlap)



if __name__ == '__main__':
    main()