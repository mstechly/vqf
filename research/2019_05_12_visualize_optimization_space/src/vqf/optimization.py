import pyquil.api as api
from pyquil.api._qvm import ForestConnection, QVM
from pyquil.device import NxDevice
from pyquil.api._quantum_computer import QuantumComputer   
from pyquil.api._compiler import QVMCompiler
from grove.pyqaoa.qaoa import QAOA
from pyquil.paulis import PauliTerm, PauliSum
from pyquil.gates import X, I
from grove.pyvqe.vqe import VQE

import scipy.optimize
import numpy as np
import networkx as nx
from functools import reduce
from sympy import Add, Mul, Number
from itertools import product
import time

from visualization import plot_energy_landscape, plot_variance_landscape, plot_optimization_trajectory

import pdb


class OptimizationEngine(object):
    """
    The optimization engine for the VQF algorithm.

    This class takes a problem encoded as clauses, further encodes it into hamiltonian
    and solves it using QAOA.

    Args:
        clauses (list): List of clauses (sympy expressions) representing the problem.
        m (int): Number to be factored. Needed only for the purpose of tagging result files.
        steps (int, optional): Number of steps in the QAOA algorithm. Default: 1
        grid_size (int, optional): The resolution of the grid for grid search. Default: None
        tol (float, optional): Parameter of BFGS optimization method. Gradient norm must be less than tol before successful termination. Default:1e-5
        gate_noise (float, optional): Specifies gate noise for qvm. Default: None.
        verbose (bool): Boolean flag, if True, information about the execution will be printed to the console. Default: False
        visualize (bool): Flag indicating if visualizations should be created. Default: False

    Attributes:
        clauses (list): See Args.
        grid_size (int): See Args.
        mapping (dict): Maps variables into qubit indices.
        qaoa_inst (object): Instance of QAOA class from Grove.
        samples (int): If noise model is active, specifies how many samples we should take for any given quantum program.
        ax (object): Matplotlib `axis` object, used for plotting optimization trajectory.

    """
    def __init__(self, clauses, m=None, steps=1, grid_size=None, tol=1e-5, gate_noise=None, verbose=False, visualize=False):
        self.clauses = clauses
        self.m = m
        self.verbose = verbose
        self.visualize = visualize
        self.gate_noise = gate_noise
        if grid_size is None:
            self.grid_size = len(clauses) + len(qubits)
        else:
            self.grid_size = grid_size

        cost_operators, mapping = self.create_operators_from_clauses()
        self.mapping = mapping
        mixing_operators = self.create_mixing_operators()
        minimizer_kwargs = {'method': 'BFGS',
                                'options': {'gtol': tol, 'disp': False}}
        if self.verbose:
            print_fun = print
        else:
            print_fun = pass_fun

        qubits = list(range(len(mapping)));

        if gate_noise:
            self.samples = int(1e3)
            pauli_channel = [gate_noise] * 3
        else:
            self.samples = None
            pauli_channel = None
        connection = ForestConnection()
        qvm = QVM(connection=connection, gate_noise=pauli_channel)
        topology = nx.complete_graph(len(qubits))
        device = NxDevice(topology=topology)
        qc = QuantumComputer(name="my_qvm",
                       qam=qvm,
                       device=device,
                       compiler=QVMCompiler(
                           device=device,
                           endpoint=connection.compiler_endpoint))

        vqe_option = {'disp': print_fun, 'return_all': True,
                      'samples': self.samples}


        self.qaoa_inst = QAOA(qc, 
                          qubits, 
                          steps=steps, 
                          init_betas=None, 
                          init_gammas=None,
                          cost_ham=cost_operators,
                          ref_ham=mixing_operators, 
                          minimizer=scipy.optimize.minimize,
                          minimizer_kwargs=minimizer_kwargs,
                          rand_seed=None,
                          vqe_options=vqe_option, 
                          store_basis=True)

        self.ax = None

    def create_operators_from_clauses(self):
        """
        Creates cost hamiltonian from clauses.
        For details see section IIC from the article.
        """
        operators = []
        mapping = {}
        variable_counter = 0
        for clause in self.clauses:
            if clause == 0:
                continue
            variables = list(clause.free_symbols)
            for variable in variables:
                if str(variable) not in mapping.keys():
                    mapping[str(variable)] = variable_counter
                    variable_counter += 1
            pauli_terms = []
            quadratic_pauli_terms = []
            if type(clause) == Add:
                clause_terms = clause.args
            elif type(clause) == Mul:
                clause_terms = [clause]
            for single_term in clause_terms:
                if len(single_term.free_symbols) == 0:
                    if self.verbose:
                        print("Constant term", single_term)
                    pauli_terms.append(PauliTerm("I", 0, int(single_term)))
                elif len(single_term.free_symbols) == 1:
                    if self.verbose:
                        print("Single term", single_term)
                    multiplier = 1
                    if type(single_term) == Mul:
                        multiplier = int(single_term.args[0])
                    symbol = list(single_term.free_symbols)[0]
                    symbol_id = mapping[str(symbol)]
                    pauli_terms.append(PauliTerm("I", symbol_id, 1/2*multiplier))
                    pauli_terms.append(PauliTerm("Z", symbol_id, -1/2*multiplier))
                elif len(single_term.free_symbols) == 2 and type(single_term) == Mul:
                    if self.verbose:
                        print("Double term", single_term)
                    multiplier = 1
                    if isinstance(single_term.args[0], Number):
                        multiplier = int(single_term.args[0])
                    symbol_1 = list(single_term.free_symbols)[0]
                    symbol_2 = list(single_term.free_symbols)[1]
                    symbol_id_1 = mapping[str(symbol_1)]
                    symbol_id_2 = mapping[str(symbol_2)]
                    pauli_term_1 = PauliTerm("I", symbol_id_1, 1/2*multiplier) - PauliTerm("Z", symbol_id_1, 1/2*multiplier)
                    pauli_term_2 = PauliTerm("I", symbol_id_2, 1/2) - PauliTerm("Z", symbol_id_2, 1/2)
                    quadratic_pauli_terms.append(pauli_term_1 * pauli_term_2)
                else:
                    Exception("Terms of orders higher than quadratic are not handled.")

            clause_operator = PauliSum(pauli_terms)
            for quadratic_term in quadratic_pauli_terms:
                clause_operator += quadratic_term
            
            squared_clause_operator = clause_operator**2
            if self.verbose:
                print("C:", clause_operator)
                print("C**2:", squared_clause_operator)
            operators.append(squared_clause_operator)


        return operators, mapping

    def create_mixing_operators(self):
        """
        Creates mixing hamiltonian. (eq. 10)
        """

        mixing_operators = []
        
        for key, value in self.mapping.items():
            mixing_operators.append(PauliSum([PauliTerm("X", value, -1.0)]))

        return mixing_operators

    def perform_qaoa(self):
        """
        Finds optimal angles for QAOA.

        Returns:
            sampling_results (Counter): Counter, where each element represents a bitstring that has been obtained.
            mapping (dict): See class description.

        """
        betas, gammas = self.simple_grid_search_angles(save_data=True)
        # betas, gammas = self.step_by_step_grid_search_angles()
        self.qaoa_inst.betas = betas
        self.qaoa_inst.gammas = gammas
        betas, gammas = self.get_angles()
        _, sampling_results = self.qaoa_inst.get_string(betas, gammas, samples=10000)    
        return sampling_results, self.mapping

    def get_angles(self):
        """
        Finds optimal angles with the quantum variational eigensolver method.
        
        It's direct copy of the function `get_angles` from Grove. I decided to copy it here
        to access to the optimization trajectory (`angles_history`).
        Returns:
            best_betas, best_gammas (np.arrays): best values of the betas and gammas found. 

        """
        stacked_params = np.hstack((self.qaoa_inst.betas, self.qaoa_inst.gammas))
        vqe = VQE(self.qaoa_inst.minimizer, minimizer_args=self.qaoa_inst.minimizer_args,
                  minimizer_kwargs=self.qaoa_inst.minimizer_kwargs)
        cost_ham = reduce(lambda x, y: x + y, self.qaoa_inst.cost_ham)
        # maximizing the cost function!
        param_prog = self.qaoa_inst.get_parameterized_program()
        result = vqe.vqe_run(param_prog, cost_ham, stacked_params, qc=self.qaoa_inst.qc,
                             **self.qaoa_inst.vqe_options)
        best_betas = result.x[:self.qaoa_inst.steps]
        best_gammas = result.x[self.qaoa_inst.steps:]
        optimization_trajectory = result.iteration_params
        energy_history = result.expectation_vals

        if self.ax is not None and self.visualize and self.qaoa_inst.steps==1:
            plot_optimization_trajectory(self.ax, optimization_trajectory)
        return best_betas, best_gammas

    def simple_grid_search_angles(self, save_data=False):
        """
        Finds optimal angles for QAOA by performing grid search on all the angles.
        This is not recommended for higher values of steps parameter, 
        since it results in grid_size**(2*steps) evaluations.

        Returns:
            best_betas, best_gammas (np.arrays): best values of the betas and gammas found. 

        """
        best_betas = None
        best_gammas = None
        best_energy = np.inf

        # For some reasons np.meshgrid returns columns in order, where values in second
        # grow slower than in the first one. This a fix to it.
        if self.qaoa_inst.steps == 1:
            column_order = [0]
        else:
            column_order = [1, 0] + list(range(2, self.qaoa_inst.steps))

        new_indices = np.argsort(column_order)
        beta_ranges = [np.linspace(0, np.pi, self.grid_size)] * self.qaoa_inst.steps
        all_betas = np.vstack(np.meshgrid(*beta_ranges)).reshape(self.qaoa_inst.steps, -1).T
        all_betas = all_betas[:, column_order]

        gamma_ranges = [np.linspace(0, 2*np.pi, self.grid_size)] * self.qaoa_inst.steps
        all_gammas = np.vstack(np.meshgrid(*gamma_ranges)).reshape(self.qaoa_inst.steps, -1).T        
        all_gammas = all_gammas[:, column_order]


        vqe = VQE(self.qaoa_inst.minimizer, minimizer_args=self.qaoa_inst.minimizer_args,
                      minimizer_kwargs=self.qaoa_inst.minimizer_kwargs)
        cost_hamiltonian = reduce(lambda x, y: x + y, self.qaoa_inst.cost_ham)
        all_energies = []
        data_to_save = []
        if save_data:
            file_name = "_".join([str(self.m), "grid", str(self.grid_size), str(time.time())]) + ".csv"
        for betas in all_betas:
            for gammas in all_gammas:
                stacked_params = np.hstack((betas, gammas))
                program = self.qaoa_inst.get_parameterized_program()
                energy = vqe.expectation(program(stacked_params), cost_hamiltonian, self.samples, self.qaoa_inst.qc)
                all_energies.append(energy)
                if self.verbose:
                    print(betas, gammas, energy, end="\r")
                if save_data:
                    data_to_save.append(np.hstack([betas, gammas, energy]))
                if energy < best_energy:
                    best_energy = energy
                    best_betas = betas
                    best_gammas = gammas
                    if self.verbose:
                        print("Lowest energy:", best_energy)
                        print("Angles:", best_betas, best_gammas)
            if save_data:
                np.savetxt(file_name, np.array(data_to_save), delimiter=",")

        if self.visualize:
            if self.qaoa_inst.steps == 1:
                self.ax = plot_energy_landscape(all_betas, all_gammas, np.array(all_energies), log_legend=True)
            else:
                plot_variance_landscape(all_betas, all_gammas, np.array(all_energies))

        return best_betas, best_gammas

    def step_by_step_grid_search_angles(self):
        """
        Finds optimal angles for QAOA by performing "step-by-step" grid search.
        It finds optimal angles by performing grid search on the QAOA instance with steps=1.
        Then it fixes these angles and performs grid search on the second pair of angles.
        This method requires steps*grid_size**2 evaluations and hence is more suitable
        for higger values of steps.

        Returns:
            best_betas, best_gammas (np.arrays): best values of the betas and gammas found. 

        """

        max_step = self.qaoa_inst.steps
        self.qaoa_inst.betas = np.array([])
        self.qaoa_inst.gammas = np.array([])
        best_betas = np.array([])
        best_gammas = np.array([])
        for current_step in range(1, max_step+1):
            if self.verbose:
                print("step:", current_step, "\n")
            beta, gamma = self.one_step_grid_search(current_step)
            best_betas = np.append(best_betas, beta)
            best_gammas = np.append(best_gammas, gamma)
            self.qaoa_inst.betas = best_betas
            self.qaoa_inst.gammas = best_gammas

        return best_betas, best_gammas

    def one_step_grid_search(self, current_step):
        """
        Grid search on n-th pair of QAOA angles, where n=current_step.

        Args:
            current_step (int): specify on which layer do we perform search.

        Returns:
            best_beta, best_gamma (floats): best values of the beta and gamma found. 
        """
        self.qaoa_inst.steps = current_step
        best_beta = None
        best_gamma = None
        best_energy = np.inf

        fixed_betas = self.qaoa_inst.betas
        fixed_gammas = self.qaoa_inst.gammas
        beta_range = np.linspace(0, np.pi, self.grid_size)
        gamma_range = np.linspace(0, 2*np.pi, self.grid_size)

        vqe = VQE(self.qaoa_inst.minimizer, minimizer_args=self.qaoa_inst.minimizer_args,
                      minimizer_kwargs=self.qaoa_inst.minimizer_kwargs)
        cost_hamiltonian = reduce(lambda x, y: x + y, self.qaoa_inst.cost_ham)
        for beta in beta_range:
            for gamma in gamma_range:
                betas = np.append(fixed_betas, beta)
                gammas = np.append(fixed_gammas, gamma)
                stacked_params = np.hstack((betas, gammas))
                program = self.qaoa_inst.get_parameterized_program()
                energy = vqe.expectation(program(stacked_params), cost_hamiltonian, self.samples, self.qaoa_inst.qc)
                print(beta, gamma, end="\r")
                if energy < best_energy:
                    best_energy = energy
                    best_beta = beta
                    best_gamma = gamma

        return best_beta, best_gamma

