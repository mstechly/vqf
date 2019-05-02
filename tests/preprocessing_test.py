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
    clause = q_0 + q_1 + p_0 + p_1 - 2*z_0 - 4*z_1 - 1
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

    ## Given
    known_expressions = {}
    q_0, q_1, q_2, z_0 = symbols('q_0 q_1 q_2 z_0')
    clause = q_0 + 2*q_1 - p_0 - 2*z_0
    ## When
    known_expressions = preprocessing.apply_z_rule_1(clause, known_expressions)
    ## Then
    assert len(known_expressions) == 0

    ## Given
    known_expressions = {}
    q_0, q_1, q_2, z_0 = symbols('q_0 q_1 q_2 z_0')
    clause = q_0 + 2*q_1 - p_0 - 4*z_0
    ## When
    known_expressions = preprocessing.apply_z_rule_1(clause, known_expressions)
    ## Then
    assert len(known_expressions) == 1
    assert known_expressions[z_0] == 0


def test_apply_z_rule_2():
    ## Given
    known_expressions = {}
    q, p, z = symbols('q p z')
    clause = q + p - 2*z
    ## When
    known_expressions = preprocessing.apply_z_rule_2(clause, known_expressions)
    ## Then
    assert known_expressions[p] == q
    assert known_expressions[z] == q

    ## Given
    known_expressions = {}
    q, p, z = symbols('q p z')
    clause = q + 2*p - 2*z
    ## When
    known_expressions = preprocessing.apply_z_rule_2(clause, known_expressions)
    ## Then
    assert known_expressions[q] == 0
    assert known_expressions[z] == p

    ## Given
    known_expressions = {}
    q, z = symbols('q z')
    clause = q + 1 - 2*z
    ## When
    known_expressions = preprocessing.apply_z_rule_2(clause, known_expressions)
    ## Then
    assert known_expressions[z] == 1
    assert known_expressions[q] == 1


def test_apply_rule_of_equality():
    ## Given
    known_expressions = {}
    x = symbols('x')
    clause = x - 1
    ## When
    known_expressions = preprocessing.apply_rule_of_equality(clause, known_expressions)
    ## Then
    assert known_expressions[x] == 1

    ## Given
    known_expressions = {}
    x = symbols('x')
    clause = x
    ## When
    known_expressions = preprocessing.apply_rule_of_equality(clause, known_expressions)
    ## Then
    assert known_expressions[x] == 0

    ## Given
    known_expressions = {}
    x, y = symbols('x y')
    clause = x*y - 1
    ## When
    known_expressions = preprocessing.apply_rule_of_equality(clause, known_expressions)
    ## Then
    assert known_expressions[x*y] == 1

def test_apply_rule_1():
    ## Given
    known_expressions = {}
    x, y = symbols('x y')
    clause = x * y - 1
    ## When
    known_expressions = preprocessing.apply_rule_1(clause, known_expressions)
    ## Then
    assert known_expressions[x] == 1
    assert known_expressions[y] == 1
