from preprocessing import create_clauses
from vqf_quantum import perform_qaoa
from sympy import Add
import pdb

def factor_number(m):
    p_dict, q_dict, z_dict, clauses = create_clauses(m, apply_preprocessing=False)
    if clauses[0] == 0 and len(set(clauses)) == 1:
        return decode_solution(p_dict, q_dict)
    qaoa_solution, mapping = perform_qaoa(clauses)
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
    m = 15
    # for m in [15, 21, 25, 33, 35, 39]:
    for m in [15]:
        print("M:", m)
        if m % 2 == 0:
            p = 2
            q = int(m / 2)
            print("The primes are:", p, "and", q)
        else:
            p, q = factor_number(m)
            print("The primes of ",m, "are:", p, "and", q)

if __name__ == '__main__':
    main()