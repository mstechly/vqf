from preprocessing import create_clauses
import pdb
import pyquil.api as api
from grove.pyqaoa.qaoa import QAOA
from pyquil.paulis import PauliTerm, PauliSum
import pyquil.quil as pq
from pyquil.gates import X, I
import scipy.optimize
import numpy as np
from grove.pyvqe.vqe import VQE
from functools import reduce
from visualization import plot_energy_landscape, plot_variance_landscape


def perform_qaoa(clauses, steps=1, grid_size=None, visualize=True):
    cost_operators, mapping = create_operators_from_clauses(clauses)
    driver_operators = create_driver_operators(mapping)

    minimizer_kwargs = {'method': 'BFGS',
                            'options': {'ftol': 1e-1, 'xtol': 1e-1,
                                        'disp': False}}

    vqe_option = {'disp': print, 'return_all': True,
                  'samples': None}

    qubits = list(range(len(mapping)));

    qvm = api.QVMConnection()
    qaoa_inst = QAOA(qvm, 
                      qubits, 
                      steps=steps, 
                      init_betas=None, 
                      init_gammas=None,
                      cost_ham=cost_operators,
                      ref_ham=driver_operators, 
                      minimizer=scipy.optimize.minimize,
                      minimizer_kwargs=minimizer_kwargs,
                      rand_seed=None,
                      vqe_options=vqe_option, 
                      store_basis=True)

    if grid_size is None:
        grid_size = len(clauses) + len(qubits)
    betas, gammas = grid_search_angles(qaoa_inst, grid_size, visualize)
    qaoa_inst.betas = betas
    qaoa_inst.gammas = gammas
    betas, gammas = qaoa_inst.get_angles()
    most_frequent_string, sampling_results = qaoa_inst.get_string(betas, gammas, samples=10000)
    return most_frequent_string, mapping


def create_operators_from_clauses(clauses, verbose=False):
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
                if verbose:
                    print("Constant term", single_term)
                pauli_terms.append(PauliTerm("I", 0, int(single_term) / 2))
            elif len(single_term.free_symbols) == 1:
                if verbose:
                    print("Single term", single_term)
                symbol = list(single_term.free_symbols)[0]
                symbol_id = mapping[str(symbol)]
                pauli_terms.append(PauliTerm("I", symbol_id, 1/2))
                pauli_terms.append(PauliTerm("Z", symbol_id, -1/2))
            elif len(single_term.free_symbols) == 2:
                if verbose:
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
        if verbose:
            print("C:", clause_operator)
            print("C**2:", squared_clause_operator)
        operators.append(squared_clause_operator)


    return operators, mapping


def create_driver_operators(mapping):
    driver_operators = []
    
    for key, value in mapping.items():
        driver_operators.append(PauliSum([PauliTerm("X", value, -1.0)]))

    return driver_operators


def grid_search_angles(qaoa_inst, grid_size=5, visualize=False):
    best_betas = None
    best_gammas = None
    best_energy = np.inf

    # For some reasons np.meshgrid returns columns in order, where values in second
    # grow slower than in the first one. This is to fix it.
    if qaoa_inst.steps == 1:
        column_order = [0]
    else:
        column_order = [1, 0] + list(range(2, qaoa_inst.steps))

    new_indices = np.argsort(column_order)
    beta_ranges = [np.linspace(0, np.pi, grid_size)] * qaoa_inst.steps
    all_betas = np.vstack(np.meshgrid(*beta_ranges)).reshape(qaoa_inst.steps, -1).T
    all_betas = all_betas[:, column_order]

    gamma_ranges = [np.linspace(0, 2*np.pi, grid_size)] * qaoa_inst.steps
    all_gammas = np.vstack(np.meshgrid(*gamma_ranges)).reshape(qaoa_inst.steps, -1).T        
    all_gammas = all_gammas[:, column_order]


    vqe = VQE(qaoa_inst.minimizer, minimizer_args=qaoa_inst.minimizer_args,
                  minimizer_kwargs=qaoa_inst.minimizer_kwargs)
    cost_hamiltonian = reduce(lambda x, y: x + y, qaoa_inst.cost_ham)
    all_energies = []
    for betas in all_betas:
        for gammas in all_gammas:
            stacked_params = np.hstack((betas, gammas))
            program = qaoa_inst.get_parameterized_program()
            energy = vqe.expectation(program(stacked_params), cost_hamiltonian, None, qaoa_inst.qvm)
            all_energies.append(energy)
            print(betas, gammas, end="\r")
            if energy < best_energy:
                best_energy = energy
                best_betas = betas
                best_gammas = gammas
                print("Lowest energy:", best_energy)
                print("Angles:", best_betas, best_gammas)
            

    if visualize:
        if qaoa_inst.steps == 1:
            plot_energy_landscape(all_betas, all_gammas, np.array(all_energies))
        else:
            plot_variance_landscape(all_betas, all_gammas, np.array(all_energies))

    return best_betas, best_gammas