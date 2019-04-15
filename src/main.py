from preprocessing import create_problem_representation
import pdb
from grove.pyqaoa.qaoa import QAOA
from pyquil.paulis import PauliTerm, PauliSum
import pyquil.quil as pq
from pyquil.gates import X, I


def factor_number(m):
    p_dict, q_dict, z_dict, clauses = create_problem_representation(m)
    cost_operators, mapping = create_operators_from_clauses(clauses)
    driver_operators = create_driver_operators(mapping)


def create_operators_from_clauses(clauses):
    operators = []
    mapping = {}
    variable_counter = 0
    for clause in clauses:
        if clause == 0:
            continue
        variables = list(clause.free_symbols)
        for variable in variables:
            if str(variable) not in mapping.keys():
                mapping[str(variable)] = variable_counter
                variable_counter += 1
        pauli_terms = []
        quadratic_pauli_terms = []
        for single_term in clause.args:
            if len(single_term.free_symbols) == 0:
                print("Constant term", single_term)
                pauli_terms.append(PauliTerm("I", 0, int(single_term) / 2))
            elif len(single_term.free_symbols) == 1:
                print("Single term", single_term)
                symbol = list(single_term.free_symbols)[0]
                symbol_id = mapping[str(symbol)]
                pauli_terms.append(PauliTerm("I", symbol_id, 1/2))
                pauli_terms.append(PauliTerm("Z", symbol_id, -1/2))
            elif len(single_term.free_symbols) == 2:
                print("Double term", single_term)
                symbol_1 = list(single_term.free_symbols)[0]
                symbol_2 = list(single_term.free_symbols)[1]
                symbol_id_1 = mapping[str(symbol_1)]
                symbol_id_2 = mapping[str(symbol_2)]
                pauli_term_1 = PauliTerm("I", symbol_id_1, 1/2) - PauliTerm("Z", symbol_id_1, 1/2)
                pauli_term_2 = PauliTerm("I", symbol_id_2, 1/2) - PauliTerm("Z", symbol_id_2, 1/2)
                quadratic_pauli_terms.append(pauli_term_1 * pauli_term_2)
            else:
                print("Terms of orders higher than quadratic are not handled.")
        clause_operator = PauliSum(pauli_terms)
        for quadratic_term in quadratic_pauli_terms:
            clause_operator += quadratic_term
        squared_clause_operator = clause_operator*clause_operator
        print("C:", clause_operator)
        print("C**2:", squared_clause_operator)
        operators.append(squared_clause_operator)


    return operators, mapping

def create_driver_operators(mapping):
    driver_operators = []
    
    for key, value in mapping.items():
        driver_operators.append(PauliSum([PauliTerm("X", value, -1.0)]))

    return driver_operators


def main():
    m = 15
    if m % 2 == 0:
        p = 2
        q = int(m / 2)
        print("The primes are:", p,"and", q)
    else:
        factor_number(m)

if __name__ == '__main__':
    main()