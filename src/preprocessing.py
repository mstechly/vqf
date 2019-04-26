import numpy as np
from sympy import Symbol, Add, Mul, Pow, factor
import pdb

def create_clauses(m_int, apply_preprocessing=True, verbose=True):
    """
    Creates clauses for the VQF algorithm.
    """

    m_dict, p_dict, q_dict, z_dict = create_initial_dicts(m_int)
    
    clauses = create_basic_clauses(m_dict, p_dict, q_dict, z_dict, apply_preprocessing)
    known_symbols = {}
    simplified_clauses = clauses


    if apply_preprocessing:
        counter = 0
        should_continue = True
        while should_continue:
            if verbose:
                print("Preprocessing iteration:", counter)
            simplified_clauses, new_known_symbols = apply_preprocessing_rules(simplified_clauses, verbose)
            if len(new_known_symbols) > len(known_symbols):
                should_continue = True
                known_symbols = new_known_symbols
            else:
                should_continue = False

            if counter == 0:
                should_continue = True
            counter += 1

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
        p_dict[i] = Symbol('p_'+str(i))

    q_dict = {}
    n_q = int(np.ceil(len(m_dict)/2))
    for i in range(n_q):
        q_dict[i] = Symbol('q_'+str(i))

    n_c = n_p + n_q - 1

    z_dict = {}
    for i in range(n_c):
        for j in range(i+1):
            if i!=j:
                if i >= n_m:
                    pass
                elif j==0:
                    pass
                else:
                    z_dict[(j, i)] = Symbol('z_'+str(j)+'_'+str(i))
            else:
                pass

    return m_dict, p_dict, q_dict, z_dict


def create_basic_clauses(m_dict, p_dict, q_dict, z_dict, apply_preprocessing=True):
    clauses = []
    n_c = len(p_dict) + len(q_dict) - 1
    for i in range(n_c):
        clause = 0
        for j in range(i+1):
            clause += q_dict.get(j, 0) * p_dict.get(i-j, 0)

        clause += -m_dict.get(i, 0)

        for j in range(i+1):
            clause += z_dict.get((j, i), 0)
        if apply_preprocessing:
            # This part exists in order to limit the number of z terms.
            if clause.func == Mul:
                max_sum = 0
            elif clause.func == Add:
                max_sum = len(clause.args) - m_dict.get(i, 0)

            if max_sum != 0:
                max_carry = int(np.floor(np.log2(max_sum)))
            else:
                max_carry = 0

            for j in range(i + max_carry + 1, n_c):
                if z_dict.get((i, j), 0) != 0:
                    z_dict[(i, j)] = 0
        
        for j in range(1, n_c):
            clause += - 2**j * z_dict.get((i, i+j), 0)
        clauses.append(clause)

    return clauses


def apply_preprocessing_rules(clauses, verbose=True):
    known_symbols = {}

    for clause in clauses:
        if clause == 0:
            continue
        clause = simplify_clause(clause, known_symbols)
        if verbose:
            print("Current clause:", clause)

        known_symbols = apply_z_rule_1(clause, known_symbols, verbose)
        clause = simplify_clause(clause, known_symbols)

        known_symbols = apply_z_rule_2(clause, known_symbols, verbose)
        clause = simplify_clause(clause, known_symbols)

        known_symbols = apply_rule_1(clause, known_symbols, verbose)
        clause = simplify_clause(clause, known_symbols)

        known_symbols = apply_rule_2(clause, known_symbols, verbose)
        clause = simplify_clause(clause, known_symbols)

        known_symbols = apply_rule_3(clause, known_symbols, verbose)
        clause = simplify_clause(clause, known_symbols)

        known_symbols = apply_rules_4_and_5(clause, known_symbols, verbose)
        clause = simplify_clause(clause, known_symbols)

        known_symbols = apply_rule_of_equality(clause, known_symbols, verbose)
        clause = simplify_clause(clause, known_symbols)


    simplified_clauses = []
    for clause in clauses:
        simplified_clause = simplify_clause(clause, known_symbols)
        # if simplified_clause != 0:
        simplified_clauses.append(simplified_clause)

    return simplified_clauses, known_symbols

def simplify_clause(clause, known_symbols):
    simplified_clause = clause.subs(known_symbols).expand()
    if simplified_clause.func == Add:
        factored_clause = factor(simplified_clause)
        if factored_clause.func == Mul:
            simplified_clause = simplified_clause / factored_clause.args[0]

        # Simplifies x**2 -> x, since the variables we use are binary.
        for term in simplified_clause.args:
            if term.func == Pow:
                simplified_clause = simplified_clause - term + term.args[0]
    return simplified_clause


def apply_z_rule_1(clause, known_symbols, verbose):
    # Example: p_1 + q_1 - 1 - 2*z_1_2 = 0
    # z12 must be equal to 0, otherwise the equation can't be satisfied
    max_non_z_sum = 0
    z_variables = {}
    for term in clause.args:
        term_variables = list(term.free_symbols)

        if len(term_variables) == 0:
            max_non_z_sum += term
        elif len(term_variables) == 1 and 'z' in str(term_variables[0]):
            # We care only for z-terms with coefficient other than 1
            if term.func == Mul:
                z_variables[term_variables[0]] = term.args[0]
        elif len(term_variables) == 1 and 'z' not in str(term_variables[0]):
            if term.func == Symbol:
                max_non_z_sum += 1
            elif term.func == Mul:
                max_non_z_sum += term.args[0]
        else:
            if len(term_variables) == 2 and len(term.args) == 2:
                max_non_z_sum += 1
            else:
                max_non_z_sum += term.args[0]

    if len(z_variables) > 0:
        for variable, coefficient in z_variables.items():
            if -coefficient > max_non_z_sum:
                if verbose:
                    print("Z rule 1 applied!", variable, "= 0")
                known_symbols[variable] = 0

    return known_symbols


def apply_z_rule_2(clause, known_symbols, verbose):
    # Example: p_1 + q_1 - 2*z_1_2 = 0
    # p1 must be equal to q1, otherwise the equation can't be satisfied
    if clause.func == Add and len(clause.args) == 3:
        z_term = None
        z_mul = None
        non_z_terms = []
        for term in clause.args:
            term_variables = list(term.free_symbols)
            if len(term_variables) == 1 and 'z' in str(term_variables[0]):
                if term.func == Mul:
                    z_mul = term.args[0]
                    z_term = term_variables[0]
            else:
                non_z_terms.append(term)

        if z_term is not None and z_mul==-2:
            if verbose:
                print("Z rule 2 applied!", non_z_terms[0], "=", non_z_terms[1])

            known_symbols[non_z_terms[0]] = non_z_terms[1]

    return known_symbols

def apply_rule_of_equality(clause, known_symbols, verbose):
    if clause.func == Add and len(clause.args) == 2:
        ## Basic rule of equality
        # Example: x - 1 = 0
        clause_variables = list(clause.free_symbols)
        # if len(clause_variables) == 1:
        if verbose:
            print("Rule of equality applied!", clause)
        known_symbols[clause.args[1]] = -clause.args[0]
    return known_symbols


def apply_rule_1(clause, known_symbols, verbose):
    clause_variables = list(clause.free_symbols)
    if clause.func == Add and len(clause.args)==2:
        if len(clause_variables) == 2:
            x = Symbol('x')
            y = Symbol('y')
            rule = x * y - 1
            substitution = clause.subs({clause_variables[0]: x, clause_variables[1]: y})
            if substitution - rule == 0:
                if verbose:
                    print("Rule 1 applied!", clause)
                known_symbols[clause_variables[0]] = 1
                known_symbols[clause_variables[1]] = 1
    return known_symbols


def apply_rule_2(clause, known_symbols, verbose):
    ## Rule 2:
    x = Symbol('x')
    y = Symbol('y')

    rule = x + y - 1
    clause_variables = list(clause.free_symbols)
    if clause.func == Add and len(clause.args) == 3 and len(clause_variables)==2:
        substitution = clause.subs({clause_variables[0]: x, clause_variables[1]: y})
        if substitution - rule == 0:
            if verbose:
                print("Rule 2 applied!", clause)
            known_symbols[clause_variables[0] * clause_variables[1]] = 0
            known_symbols[clause_variables[1]] = 1 - clause_variables[0]

    return known_symbols


def apply_rule_3(clause, known_symbols, verbose):
    ## Rule 3:
    if clause.func == Add and len(clause.args) == 2:
        if len(clause.args[0].free_symbols) == 0:
            constant_a = clause.args[0]
            if clause.args[1].func == Mul:
                constant_b = clause.args[1].args[0]
                symbol = clause.args[1] / constant_b
                if constant_a > 0 or constant_b < 0:
                    if verbose:
                        print("Rule 3 applied", clause)
                    known_symbols[symbol] = 1
    return known_symbols


def apply_rules_4_and_5(clause, known_symbols, verbose):
    ## Rule 4 & 5:
    constant = 0
    if clause.func == Mul:
        if verbose:
            print("Basic rule of x=0 applied!", clause)
        known_symbols[clause] = 0
    elif clause.func == Add:
        for part in clause.args:
            variables = list(part.free_symbols)

            if len(variables) == 0:
                constant += part

            elif len(variables) == 1:
                # This means, that the coefficient is equal to 1
                if part.func == Symbol:
                    continue
                if part.args[0] == variables[0] and part.args[1] != 0:
                    break
                elif part.args[1] == variables[0] and part.args[0] != 0:
                    break

            elif len(variables) == 2:
                # This means there is a coefficient other than 1
                if len(part.args) != 2:
                    break

        else:
            if constant == 0:
                if verbose:
                    print("Rule 4 applied!", clause)
                for part in clause.args:
                    known_symbols[part] = 0
            elif constant == len(clause.args) - 1:
                if verbose:
                    print("Rule 5 applied!", clause)
                for part in clause.args:
                    known_symbols[part] = 1
    return known_symbols


def update_dictionaries(known_symbols, p_dict, q_dict, z_dict):
    for symbol in known_symbols:
        str_symbol = str(symbol)
        symbol_type = str_symbol[0]
        if '*' in str_symbol:
            continue
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
    
    z_dict = {key:value for key, value in z_dict.items() if value != 0}

    return p_dict, q_dict, z_dict


if __name__ == '__main__':
    main()