import numpy as np
import pdb
import sympy
from sympy import symbols, Eq, FiniteSet, solveset

def create_problem_representation(m_int, apply_preprocessing=True):
    """
    Creates clauses VQF algorithm.
    """

    m_dict, p_dict, q_dict, z_dict = create_initial_dicts(m_int)

    # simplify_clauses(clauses)
    if apply_preprocessing:
        p_dict, q_dict, z_dict = apply_preprocessing_rules(m_dict, p_dict, q_dict, z_dict)

    clauses = create_clauses(m_dict, p_dict, q_dict, z_dict)

    return p_dict, q_dict, z_dict, clauses


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


def apply_preprocessing_rules(m_dict, p_dict, q_dict, z_dict):
    #Comes from the assumption that p and q are odd. If they are even, the problem is trivial.
    p_dict[0] = 1
    q_dict[0] = 1
    z_dict[(0, 1)] = 0
    return p_dict, q_dict, z_dict


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

        # for j in range(1, n_c):
        #     clause += - 2**j * z_dict.get((i, i+j), 0)
        clauses.append(clause)
    return clauses


def simplify_clauses(clauses):
    solution_set = FiniteSet(0, 1)
    for clause in clauses:
        variables = list(clause.free_symbols)
        solveset(Eq(clause, 0), variables, domain=solution_set)
        #TODO


if __name__ == '__main__':
    main()