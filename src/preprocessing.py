import numpy as np
import pdb


def create_problem_representation(m_int, preprocessing=True):
    """
    Creates clauses VQF algorithm.
    """
    m_dict, p_dict, q_dict, z_dict = create_initial_dicts(m_int)
    p_dict, q_dict, z_dict = apply_preprocessing_rules(m_dict, p_dict, q_dict, z_dict)


def create_initial_dicts(m_int):
    m_binary = bin(m_int)

    m_dict = {}
    for i, item in enumerate(m_binary[2:][::-1]):
        m_dict[i] = int(item)

    p_dict = {}
    n_p = len(m_dict)
    for i in range(n_p):
        p_dict[i] = np.nan

    q_dict = {}
    n_q = int(np.ceil(len(m_dict)/2))
    for i in range(n_q):
        q_dict[i] = np.nan

    n_c = n_p + n_q - 1

    z_dict = {}
    for i in range(n_c):
        for j in range(i+1):
            if j+1 == i:
                z_dict[(j, i)] = np.nan
            else:
                z_dict[(j, i)] = 0

    return m_dict, p_dict, q_dict, z_dict


def apply_preprocessing_rules(m_dict, p_dict, q_dict, z_dict):
    #Comes from the assumption that p and q are odd. If they are even, the problem is trivial.
    p_dict[0] = 1
    q_dict[0] = 1
    z_dict[(0, 1)] = 1
    return p_dict, q_dict, z_dict


def main():
    m = 15
    if m % 2 == 0:
        p = 2
        q = int(m / 2)
        print("The primes are:", p,"and", q)
    else:
        create_clauses(m)


if __name__ == '__main__':
    main()