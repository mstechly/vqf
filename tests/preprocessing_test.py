import pytest
from src import preprocessing
from sympy import symbols
import pdb

def test_apply_z_rule_1():
    ## Given
    known_expressions = {}
    q, p, z = symbols('q p z')
    clause = q + p - 1 - 2*z
    ## When
    known_expressions = preprocessing.apply_z_rule_1(clause, known_expressions)
    ## Then
    assert known_expressions[z] == 0

    ## Given
    known_expressions = {}
    q_0, q_1, p_0, p_1, z_0, z_1 = symbols('q_0 q_1 p_0 p_1 z_0 z_1')
    clause = q_0 + q_1 + p_0 + p_1 - 2 * z_0 - 4*z_1 - 1
    ## When
    known_expressions = preprocessing.apply_z_rule_1(clause, known_expressions)
    ## Then
    assert len(known_expressions) == 1
    assert known_expressions[z_1] == 0

    ## Given
    known_expressions = {}
    q, p, z = symbols('q p z')
    clause = q + p - 2*z
    ## When
    known_expressions = preprocessing.apply_z_rule_1(clause, known_expressions)
    ## Then
    assert len(known_expressions) == 0
