from preprocessing import create_clauses, assess_number_of_unknowns
from vqf_quantum import perform_qaoa
from sympy import Add
import pdb

def factor_number(m, true_p=None, true_q=None):
    p_dict, q_dict, z_dict, clauses = create_clauses(m, true_p, true_q, apply_preprocessing=True, verbose=False)
    number_of_uknowns, number_of_carry_bits = assess_number_of_unknowns(p_dict, q_dict, z_dict)
    
    if clauses[0] == 0 and len(set(clauses)) == 1:
        if number_of_uknowns == 0:
            return decode_solution(p_dict, q_dict)

    qaoa_solution, mapping = perform_qaoa(clauses, steps=2, grid_size=5, visualize=True)
    p_dict, q_dict, z_dict = update_dictionaries(qaoa_solution, mapping, p_dict, q_dict, z_dict)
    p, q = decode_solution(p_dict, q_dict)
    return p, q


def update_dictionaries(qaoa_solution, mapping, p_dict, q_dict, z_dict):
    values_dict = {symbol_str: qaoa_solution[index] for symbol_str, index in mapping.items()}

    for x_dict in [p_dict, q_dict, z_dict]:
        for key, value in x_dict.items():
            if str(value) in values_dict.keys():
                x_dict[key] = values_dict[str(value)]
            if type(value) == Add:
                x_dict[key] = value.subs(values_dict)
    return p_dict, q_dict, z_dict


def decode_solution(p_dict, q_dict):
    p = 0
    for key, value in p_dict.items():
        p += value * 2**key

    q = 0
    for key, value in q_dict.items():
        q += value * 2**key

    return p, q



def main():
    p_q_m_list = [[7, 5, 35], [11, 7, 77], [71, 17, 1207], [257, 131, 33667], [241, 233, 56153], [557,523, 291311]]

    for p_q_m in p_q_m_list:
        true_p = p_q_m[0]
        true_q = p_q_m[1]
        m = p_q_m[2]

        print("M:", m)
        if m % 2 == 0:
            p = 2
            q = int(m / 2)
            print("The primes are:", p, "and", q)
        else:
            p, q = factor_number(m, true_p, true_q)
            print("Calculated primes of ",m, "are:", p, "and", q)
            print("      True primes of ",m, "are:", true_p, "and", true_q)

if __name__ == '__main__':
    main()