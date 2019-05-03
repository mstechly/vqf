from preprocessing import create_clauses
from vqf_quantum import perform_qaoa
from sympy import Add
import pdb

def factor_number(m, true_p=None, true_q=None):
    p_dict, q_dict, z_dict, clauses = create_clauses(m, true_p, true_q, apply_preprocessing=True, verbose=True)
    assess_number_of_unknowns(p_dict, q_dict, z_dict)

    if clauses[0] == 0 and len(set(clauses)) == 1:
        return decode_solution(p_dict, q_dict)
    qaoa_solution, mapping = perform_qaoa(clauses, steps=1, grid_size=20, visualize=True)
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

def extract_unknowns(x_dict):
    all_values = list(x_dict.values())
    list_of_variables = []
    for x in all_values:
        if type(x) != int and len(x.free_symbols) != 0:
            list_of_variables += (list(x.free_symbols))

    unknowns = list(set(list_of_variables))
    return unknowns

def assess_number_of_unknowns(p_dict, q_dict, z_dict):
    p_unknowns = extract_unknowns(p_dict)
    q_unknowns = extract_unknowns(q_dict)
    z_unknowns = extract_unknowns(z_dict)
    all_unknowns = list(set(p_unknowns + q_unknowns + z_unknowns))
    carry_bits = [value for value in z_unknowns if 'z' in str(value)]
    print("Number of unknowns:", len(all_unknowns))
    print("Number of carry bits:", len(carry_bits))



def main():
    # for m in [35, 77, 1207, 33667, 56153, 291311]:
    p_q_m_list = [[7, 5, 35], [11, 7, 77], [71, 17, 1207], [241, 233, 56153], [257, 131, 33667], [557,523, 291311]]
    # for m in [35]:
    # p_q_m_list = [[7, 3, 21]]
    # p_q_m_list = [[71, 17, 1207]]

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