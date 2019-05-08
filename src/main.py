from preprocessing import create_clauses, calculate_number_of_unknowns
from optimization import OptimizationEngine
from sympy import Add, Symbol
import pdb


def factor_number(m, true_p, true_q, use_true_values=False):
    apply_preprocessing = True
    verbose = True
    if use_true_values:
        p_dict, q_dict, z_dict, clauses = create_clauses(m, true_p, true_q, apply_preprocessing, verbose)
    else:
        p_dict, q_dict, z_dict, clauses = create_clauses(m, None, None, apply_preprocessing, verbose)

    number_of_uknowns, number_of_carry_bits = calculate_number_of_unknowns(p_dict, q_dict, z_dict)
    print("Number of unknowns:", number_of_uknowns)
    print("Number of carry bits:", number_of_carry_bits)
    if clauses[0] == 0 and len(set(clauses)) == 1:
        if number_of_uknowns == 0:
            return decode_solution(p_dict, q_dict)

    optimization_engine = OptimizationEngine(clauses, steps=1, grid_size=10, verbose=True, visualize=True)
    sampling_results, mapping = optimization_engine.perform_qaoa()
    most_frequent_bit_string = max(sampling_results, key=lambda x: sampling_results[x])
    
    squared_overlap = calculate_squared_overlap(mapping, sampling_results, true_p, true_q, p_dict, q_dict)
    p_dict = update_dictionary(most_frequent_bit_string, mapping, p_dict)
    q_dict = update_dictionary(most_frequent_bit_string, mapping, q_dict)
    z_dict = update_dictionary(most_frequent_bit_string, mapping, z_dict)

    p, q = decode_solution(p_dict, q_dict)
    return p, q


def calculate_squared_overlap(mapping, sampling_results, true_p, true_q, p_dict, q_dict):
    p_binary_string = bin(true_p)[2:][::-1]
    q_binary_string = bin(true_q)[2:][::-1]

    p_binary = [int(char) for char in p_binary_string]
    q_binary = [int(char) for char in q_binary_string]
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

    total_overlap = 0
    total_count = 0
    print(correct_assignment)
    for bit_string, count in sampling_results.most_common():
        correct_count = 0
        for bit_id, bit_value in enumerate(bit_string):
            if bit_value == correct_assignment[bit_id]:
                correct_count += 1
        overlap = correct_count / len(bit_string) * count
        total_count += count
        print(bit_string, count, overlap)
        total_overlap += overlap
    total_overlap = total_overlap / total_count
    return total_overlap * total_overlap


def update_dictionary(qaoa_solution, mapping, x_dict):
    values_dict = {symbol_str: qaoa_solution[index] for symbol_str, index in mapping.items()}

    for key, value in x_dict.items():
        if str(value) in values_dict.keys():
            x_dict[key] = values_dict[str(value)]
        if type(value) == Add:
            x_dict[key] = value.subs(values_dict)
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
    p_q_m_list = [[7, 5, 35], [283, 11, 2893], [67, 37, 2479], [263, 263, 69169], [17, 3, 51]]
    for p_q_m in p_q_m_list:
        true_p = p_q_m[0]
        true_q = p_q_m[1]
        m = p_q_m[2]
        use_true_values = True
        if m == 35 or m == 51:
            use_true_values = False

        print("M:", m)
        if m % 2 == 0:
            p = 2
            q = int(m / 2)
            print("The primes are:", p, "and", q)
        else:
            p, q = factor_number(m, true_p, true_q, use_true_values=use_true_values)
            print("Calculated primes of ",m, "are:", p, "and", q)
            print("      True primes of ",m, "are:", true_p, "and", true_q)

if __name__ == '__main__':
    main()