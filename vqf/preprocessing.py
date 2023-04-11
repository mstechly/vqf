import numpy as np
from sympy import Symbol, Add, Mul, Pow, Number
from sympy import factor, srepr, sympify
import pdb

"""
This is script for preprocessing part of the VQF algorithm.
Notation used here refers to section IIA and IIB from the VQF article (see readme).
I specific equation is refered, it's equation from the paper.

Since the same variables are being used throughout the whole script, 
I've decided to write their description here to avoid redundancy.

** Arguments **
verbose 
Boolean flag.
If True, information about the execution will be printed to the console.

m_dict, p_dict, q_dict
Dictionaries representing numbers m, p and q.
Keys are indices of the bits.
Values are integers (0 or 1) or sympy expressions representing the variables (see eq. 1).

z_dict
Dictionary representing carry bits.
Keys are tuples of integers, where the first one is a starting bit and the second one is a target bit.
Values are the same as in case of p_dict and q_dict.

clauses
A list of clauses corresponding to the equations (2) and (3). 
They are represented as sympy expressions.

known_expressions
A dictionary of expressions which has been deducted from the clauses.
Keys are sympy expressions, either simple ones (e.g.: q_0), which represent single unknowns,
or more complex (e.g.: p_1*q_1). 
Values are either sympy expressions or integers (0 or 1).


"""


def create_clauses(m_int, true_p_int=None, true_q_int=None, apply_preprocessing=True, verbose=True):
    """
    Creates clauses for the VQF algorithm.

    If true_p_int or true_q_int are provided, algorithm treats their length as known.
    This also means that it sets the leading bits to 1.
    It returns dictionaries, which represent p, q and carry bits, and a list of clauses. 
    For p and q keys represent the bit index, for carry bits, between which bits carrying occurs.
    Clauses are sympy expressions which represent the optimization problem.

    Args:
        m_int (int): number to be factored, as an integer.
        true_p_int (int, optional): p - first factor, as an integer. Default: None.
        true_p_int (int, optional): q - second factor, as an integer. Default: None.
        apply_preprocessing (bool, optional): If True, the preprocessing will be applied. Default: True
        verbose (bool, optional): See module documentation at the top.

    Returns:
        p_dict, q_dict, z_dict: See module documentation at the top.
        final_clauses (list): See 'clauses' in module documentation at the top.
    """

    m_dict, p_dict, q_dict, z_dict = create_initial_dicts(m_int, true_p_int, true_q_int)
    if apply_preprocessing:
        q_dict[0] = 1
        p_dict[0] = 1
        if len(q_dict) == 2:
            q_dict[1] = 1
    clauses = create_basic_clauses(m_dict, p_dict, q_dict, z_dict, apply_preprocessing)


    if apply_preprocessing:
        simplified_clauses, known_expressions = simplify_clauses(clauses, verbose)
        p_dict, q_dict, z_dict = update_dictionaries(known_expressions, p_dict, q_dict, z_dict)
    
    known_symbols = create_known_symbols_dict(p_dict, q_dict, z_dict)

    final_clauses = []
    for clause in clauses:
        final_clauses.append(simplify_clause(clause, known_symbols))

    z_dict = {key:value for key, value in z_dict.items() if value != 0}

    # TODO: In principle this should be recurrent / in while loop.
    if apply_preprocessing and final_clauses[0] == 0 and len(set(final_clauses)) == 1:
        number_of_unknowns, _ = calculate_number_of_unknowns(p_dict, q_dict, z_dict)
        if number_of_unknowns != 0:
            p_dict, q_dict = simplify_symmetric_case(p_dict, q_dict)
            number_of_unknowns, _ = calculate_number_of_unknowns(p_dict, q_dict, z_dict)
            if number_of_unknowns != 0:
                final_clauses = create_basic_clauses(m_dict, p_dict, q_dict, z_dict, apply_preprocessing)
                final_clauses, known_expressions = simplify_clauses(final_clauses, verbose)
                p_dict, q_dict, z_dict = update_dictionaries(known_expressions, p_dict, q_dict, z_dict)

    if final_clauses[0] == 0 and len(set(final_clauses)) == 1:
        number_of_unknowns, _ = calculate_number_of_unknowns(p_dict, q_dict, z_dict)
        if number_of_unknowns != 0:
            # This should be fixed after implementation of the TODO above.
            raise Exception("All clauses equal to 0, but unknowns still exist.")

    for clause in final_clauses:
        if isinstance(clause, Number) and clause != 0:
            raise Exception("Clause is a number and is not equal to 0!")

    if verbose:
        print("Final clauses:")
        for clause in final_clauses:
            print(clause)

    return p_dict, q_dict, z_dict, final_clauses


def simplify_clauses(clauses, verbose=True):
    """
    Performs simplification of clauses.

    Args:
        clauses (list): See module documentation at the top.
        verbose (bool, optional): See module documentation at the top.

    Returns:
        simplified_clauses (list): See 'clauses' in module documentation at the top.
        known_expressions (dict): See module documentation at the top.
    """

    known_expressions = {}
    counter = 0
    should_continue = True
    simplified_clauses = clauses
    while should_continue:
        if verbose:
            print("Preprocessing iteration:", counter)
        new_simplified_clauses, new_known_expressions = apply_preprocessing_rules(simplified_clauses, verbose)
        for new_clause, old_clause in zip(new_simplified_clauses, simplified_clauses):
            if new_clause != old_clause:
                break
        else:
            should_continue = False

        simplified_clauses = new_simplified_clauses
        known_expressions = {**known_expressions, **new_known_expressions}

        if counter == 0:
            should_continue = True
        counter += 1
        if verbose:
            print("\n")

    return simplified_clauses, known_expressions


def create_initial_dicts(m_int, true_p_int=None, true_q_int=None):
    """
    Creates dictionaries representing m, p and q.

    If true_p_int or true_q_int are provided, algorithm treats their length as known 
    and it sets the leading bits to 1.

    Args:
        m_int (int): number to be factored, as an integer.
        true_p_int (int, optional): p - first factor, as an integer. Default: None.
        true_p_int (int, optional): q - second factor, as an integer. Default: None.

    Returns:
        m_dict, p_dict, q_dict, z_dict: See module documentation at the top.
    """

    m_binary = bin(m_int)[2:][::-1]

    m_dict = {}
    n_m = len(m_binary)
    for i, item in enumerate(m_binary):
        m_dict[i] = int(item)

    p_dict = {}
    if true_p_int is None:
        n_p = len(m_dict)
    else:
        true_p_binary = bin(true_p_int)[2:][::-1]
        n_p = len(true_p_binary)

    if true_q_int is None:
        n_q = int(np.ceil(len(m_dict)/2))
    else:
        true_q_binary = bin(true_q_int)[2:][::-1]
        n_q = len(true_q_binary)


    for i in range(n_p):
        p_dict[i] = Symbol('p_'+str(i))

    if true_p_int is not None:
        p_dict[n_p-1] = 1

    q_dict = {}
    for i in range(n_q):
        q_dict[i] = Symbol('q_'+str(i))

    if true_q_int is not None:
        q_dict[n_q-1] = 1


    n_c = len(m_dict) + int(np.ceil(len(m_dict)/2)) - 1

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
    """
    Creates clauses based on dictionaries representing m, p, q and carry bits.
    
    Clauses created here are "basic", which means they are bases solely on the 
    information about the number of bits.
    Preprocessing in this case is limited to calculating what carry bits are possible
    for p and q of given length.

    Args:
        m_dict, p_dict, q_dict, z_dict: See module documentation at the top.
        apply_preprocessing (bool, optional): If True, the preprocessing will be applied. Default: True

    Returns:
        clauses (list): See module documentation at the top.
    """
    clauses = []
    n_c = len(m_dict) + int(np.ceil(len(m_dict)/2)) - 1
    for i in range(n_c):
        clause = 0
        for j in range(i+1):
            clause += q_dict.get(j, 0) * p_dict.get(i-j, 0)
        clause += -m_dict.get(i, 0)

        for j in range(i+1):
            clause += z_dict.get((j, i), 0)
        
        if type(clause) == int:
            clause = sympify(clause)

        if apply_preprocessing and clause != 0:
            # This part exists in order to limit the number of z terms.
            max_sum = get_max_sum_from_clause(clause)


            if max_sum != 0:
                max_carry = int(np.floor(np.log2(max_sum)))
            else:
                max_carry = 0
            for j in range(i + max_carry + 1, n_c):
                if z_dict.get((i, j), 0) != 0:
                    z_dict[(i, j)] = 0
        
        for j in range(1, n_c):
            clause += - 2**j * z_dict.get((i, i+j), 0)
        
        if clause == 0:
            clause = sympify(clause)

        clauses.append(clause)

    return clauses


def get_max_sum_from_clause(clause):
    """
    Calculates maximum sum that can be achieved for given clause.
    
    Args:
        clause: sympy expression representing a clause.

    Returns:
        max_sum (int): Maximum sum that can be achieved by given clause.
    """

    max_sum = 0
    if clause.func == Mul:
        if isinstance(clause.args[0], Number) and clause.args[0] > 0:
            max_sum += int(clause.args[0])
        else:
            max_sum += 1
    elif clause.func == Add:
        for term in clause.args:
            if isinstance(term, Number):
                max_sum += int(term)
            elif term.func == Symbol:
                max_sum += 1
            elif term.func == Mul:
                if isinstance(term.args[0], Number) and term.args[0] > 0:
                    max_sum += int(term.args[0])
                elif isinstance(term.args[0], Number) and term.args[0] < 0:
                    pass
                else:
                    max_sum += 1

    elif clause.func == Symbol:
        max_sum = 1
    elif isinstance(clause, Number):
        max_sum += int(clause)
    return max_sum


def simplify_symmetric_case(p_dict, q_dict):
    """
    Simplifies expressions if p_dict and q_dict are "symmetric".
    
    Example: numbers p = [1, x, 1] and q = [1, 1-x, 1] are symmetric, since 
    if multiplied, they will produce the same result regardless of the value of x.
    This function handles such cases and sets the bits in such way, 
    that p is greater than q.

    Args:
        p_dict, q_dict: See module documentation at the top.

    Returns:
        p_dict, q_dict: See module documentation at the top.
    """

    if len(p_dict) != len(q_dict):
        return p_dict, q_dict

    for key in sorted(p_dict.keys())[::-1]:
        if type(p_dict[key]) != int or type(q_dict[key]) != int:
            if p_dict[key] + q_dict[key] == 1:
                known_expressions = {p_dict[key]: 1, q_dict[key]: 0}
                p_dict, q_dict, _ = update_dictionaries(known_expressions, p_dict, q_dict, {})
    return p_dict, q_dict


def apply_preprocessing_rules(clauses, verbose=True):
    """
    Applies preprocessing rules to given set of clauses.

    Preprocessing rules are those mentioned in section IIB of the paper,
    especially equation (5), extended by some additional rules (see readme).

    Args:
        clauses (list): See module documentation at the top.
        verbose (bool, optional): See module documentation at the top.

    Returns:
        simplified_clauses (list): See 'clauses' in module documentation at the top.
        known_expressions (dict): See module documentation at the top.
    """

    known_expressions = {}
    counter = 0 

    for clause in clauses:
        clause = simplify_clause(clause, known_expressions)
        if verbose and clause != 0:
            print("Current clause", counter, ":", clause)
        counter += 1
        if clause == 0:
            continue


        known_expressions = apply_z_rule(clause, known_expressions, verbose)
        clause = simplify_clause(clause, known_expressions)

        known_expressions = apply_rule_1(clause, known_expressions, verbose)
        clause = simplify_clause(clause, known_expressions)

        known_expressions = apply_rule_2(clause, known_expressions, verbose)
        clause = simplify_clause(clause, known_expressions)

        known_expressions = apply_rule_3(clause, known_expressions, verbose)
        clause = simplify_clause(clause, known_expressions)

        known_expressions = apply_rules_4_and_5(clause, known_expressions, verbose)
        clause = simplify_clause(clause, known_expressions)

        known_expressions = apply_rule_of_equality(clause, known_expressions, verbose)
        clause = simplify_clause(clause, known_expressions)

        known_expressions = apply_parity_rule(clause, known_expressions, verbose)
        clause = simplify_clause(clause, known_expressions)


    simplified_clauses = []
    for clause in clauses:
        simplified_clause = simplify_clause(clause, known_expressions)
        simplified_clauses.append(simplified_clause)
    return simplified_clauses, known_expressions


def simplify_clause(clause, known_expressions, iterations=2):
    """
    Simplifies clauses based on known_expressions and algebraic rules.
    
    Substitutes some variables in given claus with known expressions (if possible).
    Also performs addition simplification, like dividing clause by a constant (if it makes sense)
    and substituting x**2 -> x, since the variables are binary.

    TODO: instead of using iterations, it should use some form of recursion.
    Args:
        clause: sympy expression representing a clause.
        known_expressions (dict): See module documentation at the top.

    Returns:
        simplified_clause: sympy expression representing a simplified clause.
    """
    simplified_clause = clause
    for i in range(iterations):
        simplified_clause = simplified_clause.subs(known_expressions).expand()
        if simplified_clause.func == Add:
            # Simplifies x**2 -> x, since the variables we use are binary.
            for term in simplified_clause.args:
                if term.func == Mul and 'Pow' in srepr(term):
                    for subterm in term.args:
                        if subterm.func == Pow:
                            simplified_clause = simplified_clause.subs({subterm: subterm.args[0]})

                if term.func == Pow:
                    simplified_clause = simplified_clause - term + term.args[0]

            # factor() is very resource-heavy - this intends to limit its usage.
            # It gives even 20x speedup for large numbers!
            for term in simplified_clause.args:
                if term.func == Mul or isinstance(term, Number):
                    continue
                else:
                    break
            else:
                factored_clause = factor(simplified_clause)
                if factored_clause.func == Mul:
                    if isinstance(factored_clause.args[0], Number):
                        simplified_clause = simplified_clause / factored_clause.args[0]

    return simplified_clause


def apply_z_rule(clause, known_expressions, verbose=False):
    """
    Extends known_expressions by applying "Z-rule" (see example below).
    
    This rule has been described in the section IIB of the article.
    Example: p_1 + q_1 - 1 - 2*z_1_2 = 0
    z12 must be equal to 0, otherwise the equation can't be satisfied
    
    Args:
        clause: sympy expression representing a clause.
        known_expressions (dict): See module documentation at the top.
        verbose (bool): See module documentation at the top.

    Returns:
        known_expressions (dict): See module documentation at the top.
    """

    # TODO: The following equations should add the following rule z_2_3*z_1_3 = 0
    # TODO: p_1 + p_2 + p_3 + p_4 - 2*z_2_3 - 4*z_1_3 = 0
    max_sum = get_max_sum_from_clause(clause)

    negative_terms = []
    for term in clause.args:
        if term.func == Mul and isinstance(term.args[0], Number) and term.args[0] < 0:
            negative_terms.append(term)

    if len(negative_terms) > 0:
        for term in negative_terms:
            if -term.args[0] > max_sum:
                variable = term / term.args[0]
                if verbose:
                    print("Z rule 1 applied!", variable, "= 0")
                known_expressions[variable] = 0

    return known_expressions


def apply_parity_rule(clause, known_expressions, verbose=False):
    """
    Extends known_expressions by applying parity rule (see example below).
    
    Example: p_1 + q_1 - 2*z_1_2 = 0
    p1 and z_1_2 must be equal to q1, otherwise the equation can't be satisfied.
    For more examples please refer to the tests.
    This rule turned out to be the most problematic one. 
    There are known cases it doesn't solve (see comments) and it might require
    refactoring/reworking.
    
    Args:
        clause: sympy expression representing a clause.
        known_expressions (dict): See module documentation at the top.
        verbose (bool): See module documentation at the top.

    Returns:
        known_expressions (dict): See module documentation at the top.
    """

    new_known_expressions = {}
    even_positive_terms = []
    even_negative_terms = []
    odd_terms = []
    if clause.func == Add:
        for term in clause.args:
            if term.func == Symbol:
                odd_terms.append(term)
            if isinstance(term, Number):
                if term % 2 == 0 and term > 0:
                    even_positive_terms.append(term)
                elif term % 2 == 0 and term < 0:
                    even_negative_terms.append(term)
                else:
                    odd_terms.append(term)
            if term.func == Mul:
                first_argument = term.args[0]
                if isinstance(first_argument, Number):
                    if first_argument % 2 == 0 and first_argument > 0:
                        even_positive_terms.append(term)
                    elif first_argument % 2 == 0 and first_argument < 0:
                        even_negative_terms.append(term)
                    else:
                        odd_terms.append(term)
                else:
                    odd_terms.append(term)

    if len(odd_terms) == 1:
        if type(odd_terms[0]) == Symbol:
            new_known_expressions[odd_terms[0]] = 0
        elif type(odd_terms[0]) == Mul:
            term = odd_terms[0]
            if isinstance(term.args[0], Number):
                term = term / term.args[0]
            new_known_expressions[term] = 0
        else:
            print("TODO: Z rule 2: don't know this type!")
            pdb.set_trace()

    if len(odd_terms) == 2:
        non_number_index = None
        if isinstance(odd_terms[0], Number):
            non_number_index = 1
        elif isinstance(odd_terms[1], Number):
            non_number_index = 0
        if non_number_index is not None:
            term = odd_terms[non_number_index]
            if type(term) == Symbol:
                new_known_expressions[term] = 1
            elif type(term) == Mul:
                for arg in term.args:
                    if not isinstance(arg, Number):
                        new_known_expressions[arg] = 1
            else:
                # TODO: Example of clause which results in this case:
                # 2*q_2 + z_4_6 + z_5_6 - 4
                # (p=23, q=23, m=529)
                # This should be handled by rule 5
                # print("TODO: Z rule 2: don't know this type!")
                # pdb.set_trace()
                pass

        else:
            if 'q' in str(odd_terms[0]):
                non_q_index = 1
            else:
                non_q_index = 0
            variable_0 = odd_terms[1 - non_q_index]
            variable_1 = odd_terms[non_q_index]

            if type(variable_0) == Mul:
                if isinstance(variable_0.args[0], Number):
                    variable_0 = variable_0/variable_0.args[0]

            if type(variable_1) == Mul:
                if isinstance(variable_1.args[0], Number):
                    variable_1 = variable_1/variable_1.args[0]


            new_known_expressions[variable_1] = variable_0

            if len(even_negative_terms) == 1:
                term = even_negative_terms[0]
                if isinstance(term, Number):
                    # TODO: Example of clause which results in this case:
                    # q_2 + z_5_6 + 2*z_7_8 - 2
                    # (p=29, q=23, m=667)
                    # pdb.set_trace()
                    pass
                elif type(term) == Mul:
                    if len(even_positive_terms) == 0:
                        term = term / term.args[0]
                        new_known_expressions[term] = variable_0
                    else:
                        pass
                        # pdb.set_trace()
                else:
                    print("TODO: Z rule 2: don't know this type!")
                    pdb.set_trace()

    if len(odd_terms) == 3:
        number_index = None
        if isinstance(odd_terms[0], Number):
            number_index = 0
        elif isinstance(odd_terms[1], Number):
            number_index = 1
        elif isinstance(odd_terms[2], Number):
            number_index = 2
        if number_index is not None:
            indices = [0, 1, 2]
            indices.remove(number_index)
            new_term = odd_terms[indices[0]] * odd_terms[indices[1]]
            if isinstance(new_term.args[0], Number):
                new_term = new_term / new_term.args[0]
            if 'Pow' in srepr(new_term):
                new_term = simplify_clause(new_term, {})
            new_known_expressions[new_term] = 0
            

 

    if len(new_known_expressions) != 0:
        known_expressions = {**known_expressions, **new_known_expressions}
        if verbose:
            print("Z rule 2 applied:", new_known_expressions)


    return known_expressions


def apply_rule_of_equality(clause, known_expressions, verbose=False):
    """
    Extends known_expressions by leveraging that clause equals to 0.
    
    Example:  x - 1 = 0
    
    Args:
        clause: sympy expression representing a clause.
        known_expressions (dict): See module documentation at the top.
        verbose (bool): See module documentation at the top.

    Returns:
        known_expressions (dict): See module documentation at the top.
    """
    if clause.func == Symbol:
        known_expressions[clause] = 0
    elif clause.func == Add and len(clause.args) == 2:
        if isinstance(clause.args[0], Number):
            known_expressions[clause.args[1]] = -clause.args[0]
        elif isinstance(clause.args[1], Number):
            known_expressions[clause.args[0]] = -clause.args[1]
        else:
            if 'q' in str(clause.args[0]):
                non_q_index = 1
            else:
                non_q_index = 0
            if '-' in str(clause.args[1 - non_q_index]):
                known_expressions[clause.args[non_q_index]] = -clause.args[1 - non_q_index]
            else:
                known_expressions[-clause.args[non_q_index]] = clause.args[1 - non_q_index]

    elif clause.func == Mul:
        if len(clause.free_symbols) == 1:
            known_expressions[list(clause.free_symbols)[0]] = 0
        else:
            known_expressions[clause] = 0
    else:
        return known_expressions
    if verbose:
        print("Rule of equality applied!", clause)
    return known_expressions


def apply_rule_1(clause, known_expressions, verbose=False):
    """
    Extends known_expressions by applying rule 1 from eq. (5).
        
    Args:
        clause: sympy expression representing a clause.
        known_expressions (dict): See module documentation at the top.
        verbose (bool): See module documentation at the top.

    Returns:
        known_expressions (dict): See module documentation at the top.
    """

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
                known_expressions[clause_variables[0]] = 1
                known_expressions[clause_variables[1]] = 1
    return known_expressions


def apply_rule_2(clause, known_expressions, verbose=False):
    """
    Extends known_expressions by applying rule 2 from eq. (5).
        
    Args:
        clause: sympy expression representing a clause.
        known_expressions (dict): See module documentation at the top.
        verbose (bool): See module documentation at the top.

    Returns:
        known_expressions (dict): See module documentation at the top.
    """
    x = Symbol('x')
    y = Symbol('y')

    rule = x + y - 1
    clause_variables = list(clause.free_symbols)
    if clause.func == Add and len(clause.args) == 3 and len(clause_variables)==2:
        substitution = clause.subs({clause_variables[0]: x, clause_variables[1]: y})
        if substitution - rule == 0:
            if verbose:
                print("Rule 2 applied!", clause_variables[0], "=", 1 - clause_variables[1])
            known_expressions[clause_variables[0] * clause_variables[1]] = 0
            if 'q' in str(clause_variables[1]):
                known_expressions[clause_variables[0]] = 1 - clause_variables[1]
            else:
                known_expressions[clause_variables[1]] = 1 - clause_variables[0]

    return known_expressions


def apply_rule_3(clause, known_expressions, verbose=False):
    """
    Extends known_expressions by applying rule 3 from eq. (5).
        
    Args:
        clause: sympy expression representing a clause.
        known_expressions (dict): See module documentation at the top.
        verbose (bool): See module documentation at the top.

    Returns:
        known_expressions (dict): See module documentation at the top.
    """
    if clause.func == Add and len(clause.args) == 2:
        if len(clause.args[0].free_symbols) == 0:
            constant_a = clause.args[0]
            if clause.args[1].func == Mul:
                constant_b = clause.args[1].args[0]
                symbol = clause.args[1] / constant_b
                if constant_a > 0 or constant_b < 0:
                    if verbose:
                        print("Rule 3 applied", clause)
                    known_expressions[symbol] = 1
    return known_expressions


def apply_rules_4_and_5(clause, known_expressions, verbose=False):
    """
    Extends known_expressions by applying rules 4 and 5 from eq. (5).
        
    Args:
        clause: sympy expression representing a clause.
        known_expressions (dict): See module documentation at the top.
        verbose (bool): See module documentation at the top.

    Returns:
        known_expressions (dict): See module documentation at the top.
    """
    constant = 0
    if clause.func == Add:
        for term in clause.args:
            variables = list(term.free_symbols)

            if len(variables) == 0:
                constant += term

            elif len(variables) == 1:
                # This means, that the coefficient is equal to 1
                if term.func == Symbol:
                    continue
                if term.args[0] == variables[0] and term.args[1] != 0:
                    break
                elif term.args[1] == variables[0] and term.args[0] != 0:
                    break

            elif len(variables) == 2:
                # This means there is a coefficient other than 1
                if len(term.args) != 2:
                    break

        else:
            if constant == 0:
                if verbose:
                    print("Rule 4 applied!", clause)
                for term in clause.args:
                    known_expressions[term] = 0
            elif constant == -(len(clause.args) - 1):
                if verbose:
                    print("Rule 5 applied!", clause)
                for term in clause.args:
                    if term != constant:
                        known_expressions[term] = 1
    return known_expressions


def update_dictionaries(known_expressions, p_dict, q_dict, z_dict):
    """
    Updates values of p_dict, q_dict and z_dict with the values stored in known_expressions.
        
    Args:
        known_expressions (dict): See module documentation at the top.
        p_dict, q_dict, z_dict: See module documentation at the top.

    Returns:
        known_expressions (dict): See module documentation at the top.
    """

    all_known_expressions = {**known_expressions}

    for symbol in known_expressions:
        str_symbol = str(symbol)
        symbol_type = str_symbol[0]
        if '*' in str_symbol:
            continue
        if symbol_type == 'p':
            symbol_number = int(str_symbol.split('_')[1])
            p_dict[symbol_number] = known_expressions[symbol]
        if symbol_type == 'q':
            symbol_number = int(str_symbol.split('_')[1])
            q_dict[symbol_number] = known_expressions[symbol]

        if symbol_type == 'z':
            symbol_number_0 = int(str_symbol.split('_')[1])
            symbol_number_1 = int(str_symbol.split('_')[2])
            z_dict[(symbol_number_0, symbol_number_1)] = known_expressions[symbol]


    known_symbols = create_known_symbols_dict(p_dict, q_dict, z_dict)
    all_known_expressions = {**all_known_expressions, **known_symbols}

    for x_dict in [p_dict, q_dict, z_dict]:
        for index, value in x_dict.items():
            if type(value) in [Symbol, Add, Mul]:
                x_dict[index] = x_dict[index].subs(all_known_expressions)


    return p_dict, q_dict, z_dict


def create_known_symbols_dict(p_dict, q_dict, z_dict):
    """
    Creates a dictionary of known simple symbols.
        
    Args:
        p_dict, q_dict, z_dict: See module documentation at the top.

    Returns:
        known_symbols (dict): Dictionary of known symbols.
    """

    known_symbols = {}
    for index, value in p_dict.items():
        known_symbols[Symbol('p_' + str(index))] = value

    for index, value in q_dict.items():
        known_symbols[Symbol('q_' + str(index))] = value

    for index, value in z_dict.items():
        known_symbols[Symbol('z_' + str(index[0]) + "_" + str(index[1]))] = value
    return known_symbols


def calculate_number_of_unknowns(p_dict, q_dict, z_dict):
    """
    Calculates how many unknowns and unknown carry bits are in given dictionaries.

    Args:
        p_dict, q_dict, z_dict: See module documentation at the top.

    Returns:
        number_of_unknowns (int): Number of unknowns.
        number_of_carry_bits (int): Number of unknown carry bits.
    """

    p_unknowns = extract_unknowns(p_dict)
    q_unknowns = extract_unknowns(q_dict)
    z_unknowns = extract_unknowns(z_dict)
    all_unknowns = list(set(p_unknowns + q_unknowns + z_unknowns))
    non_carry_unknowns = p_unknowns + q_unknowns
    carry_bits = [value for value in z_unknowns if 'z' in str(value) and value not in non_carry_unknowns]
    number_of_unknowns = len(all_unknowns)
    number_of_carry_bits = len(carry_bits)
    return number_of_unknowns, number_of_carry_bits


def extract_unknowns(x_dict):
    """
    Extracts unknown variable (sympy expressions) from given dictionary
        
    Args:
        x_dict: One of q_dict, p_dict, z_dict, see module documentation at the top.

    Returns:
        unknowns: list of unknown variables in given dictionary
    """

    all_values = list(x_dict.values())
    list_of_variables = []
    for x in all_values:
        if type(x) != int and len(x.free_symbols) != 0:
            list_of_variables += (list(x.free_symbols))

    unknowns = list(set(list_of_variables))
    return unknowns


def factor_56153():
    clauses = []
    p_3 = Symbol('p_3')
    q_3 = Symbol('q_3')
    p_4 = Symbol('p_4')
    q_4 = Symbol('q_4')

    clauses.append(p_3 + q_3 - 1)
    clauses.append(p_4 + q_4 - 1)
    clauses.append(p_4*q_3 + p_3*q_4 - 1)
    p_dict = {0: 1, 1: 0, 2: 0, 3: p_3, 4: p_4, 5: 1, 6: 1, 7: 1}
    q_dict = {0: 1, 1: 0, 2: 0, 3: q_3, 4: q_4, 5: 1, 6: 1, 7: 1}
    z_dict = {}
    return p_dict, q_dict, z_dict, clauses


def factor_291311():
    clauses = []
    p_1 = Symbol('p_1')
    p_2 = Symbol('p_2')
    p_5 = Symbol('p_5')
    q_1 = Symbol('q_1')
    q_2 = Symbol('q_2')
    q_5 = Symbol('q_5')

    clauses.append(p_1 + q_1 - 1)
    clauses.append(p_2 + q_2 - 1)
    clauses.append(p_5 + q_5 - 1)
    clauses.append(p_1*q_2 + p_2*q_1 - 1)
    clauses.append(p_2*q_5 + p_5*q_2 - 1)
    clauses.append(p_5*q_1 + q_5*p_1 - 1)

    p_dict = {0: 1, 1: p_1, 2: p_2, 3: 1, 4: 0, 5: p_5, 6: 0, 7: 0, 8: 0, 9: 1}
    q_dict = {0: 1, 1: q_1, 2: q_2, 3: 1, 4: 0, 5: q_5, 6: 0, 7: 0, 8: 0, 9: 1}
    z_dict = {}
    return p_dict, q_dict, z_dict, clauses