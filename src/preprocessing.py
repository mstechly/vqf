import numpy as np
import pdb
import sympy
from sympy import symbols, Eq, FiniteSet, solveset


def create_problem_representation(m_int, apply_preprocessing=True):
    """
    Creates clauses VQF algorithm.
    """

    m_dict, p_dict, q_dict, z_dict = create_initial_dicts(m_int)

    
    clauses = create_clauses(m_dict, p_dict, q_dict, z_dict)
    if apply_preprocessing:
        known_symbols = {}
        simplified_clauses = clauses
        counter = 0
        while len(known_symbols) > 0 or counter == 0:
            print("Preprocessing iteration:", counter)
            counter += 1
            simplified_clauses, known_symbols = apply_preprocessing_rules(simplified_clauses)

    p_dict, q_dict, z_dict = update_dictionaries(known_symbols, p_dict, q_dict, z_dict)
    pdb.set_trace()
    return p_dict, q_dict, z_dict, simplified_clauses


def create_initial_dicts(m_int):
    m_binary = bin(m_int)[2:][::-1]

    m_dict = {}
    n_m = len(m_binary)
    for i, item in enumerate(m_binary):
        m_dict[i] = int(item)

    p_dict = {}
    n_p = len(m_dict)
    for i in range(n_p):
        p_dict[i] = symbols('p_'+str(i))

    q_dict = {}
    n_q = int(np.ceil(len(m_dict)/2))
    for i in range(n_q):
        q_dict[i] = symbols('q_'+str(i))

    n_c = n_p + n_q - 1

    z_dict = {}
    for i in range(n_c):
        for j in range(i+1):
            if i!=j:
                if i >= n_m:
                    z_dict[(j, i)] = 0
                elif j==0:
                    z_dict[(j, i)] = 0
                else:
                    z_dict[(j, i)] = symbols('z_'+str(j)+'_'+str(i))
            else:
                z_dict[(j, i)] = 0

    return m_dict, p_dict, q_dict, z_dict


def create_clauses(m_dict, p_dict, q_dict, z_dict):
    clauses = []
    n_c = len(p_dict) + len(q_dict) - 1
    for i in range(n_c):
        clause = 0
        for j in range(i+1):
            clause += q_dict.get(j,0)*p_dict.get(i-j, 0)

        for j in range(i+1):
            clause += z_dict.get((j, i), 0)

        clause += -m_dict.get(i, 0)

        for j in range(1, n_c):
            clause += - 2**j * z_dict.get((i, i+j), 0)
        clauses.append(clause)

    return clauses


def apply_preprocessing_rules(clauses):
    solution_set = FiniteSet(0, 1)
    simple_rules = []
    x, y, a, b = symbols('x y a b')
    simple_rules.append([x * y - 1, {'x': 1, 'y': 1}])
    simple_rules.append([x + y - 1, {'xy': 0}])
    simple_rules.append([a - b * x, {'x': 1}])
    known_symbols = {}

    for clause in clauses:
        clause.subs(known_symbols)
        clause_variables = list(clause.free_symbols)
        
        ## Rule 1:
        rule = x * y - 1
        if len(clause_variables) == 2:
            substitution = clause.subs({clause_variables[0]: x, clause_variables[1]: y})
            if substitution - rule == 0:
                print("Rule 1 applied!")
                known_symbols[clause_variables[0]] = 1
                known_symbols[clause_variables[1]] = 1
                continue

        ## Rule 2:
        rule = x + y - 1
        if len(clause_variables) == 2:
            substitution = clause.subs({clause_variables[0]: x, clause_variables[1]: y})
            if substitution - rule == 0:
                print("Rule 2 applied!")
                known_symbols[clause_variables[0] * clause_variables[1]] = 0
                continue

    simplified_clauses = []
    for clause in clauses:
        simplified_clause = clause.subs(known_symbols)
        if simplified_clause != 0:
            simplified_clauses.append(simplified_clause)

    return simplified_clauses, known_symbols


def update_dictionaries(known_symbols, p_dict, q_dict, z_dict):
    for symbol in known_symbols:
        str_symbol = str(symbol)
        symbol_type = str_symbol[0]
        if symbol_type == 'p':
            symbol_number = int(str_symbol.split('_')[1])
            p_dict[symbol_number] = known_symbols[symbol]
        if symbol_type == 'q':
            symbol_number = int(str_symbol.split('_')[1])
            q_dict[symbol_number] = known_symbols[symbol]
        if symbol_type == 'z':
            symbol_number_0 = int(str_symbol.split('_')[1])
            symbol_number_1 = int(str_symbol.split('_')[2])
            z_dict[(symbol_number_0, symbol_number_1)] = known_symbols[symbol]

    return p_dict, q_dict, z_dict


if __name__ == '__main__':
    main()