from sympy import simplify, diff, Expr, Derivative
from sympy.abc import x as X

from src.utils import bolzano


def bisection(function, a, b, tolerance):
    if not bolzano(function, a, b):
        raise Exception(
            f"invalid initial interval -> bolzano does not apply in this interval [{min(a, b)};{max(a, b)}]")
    iteration = 1
    simplified = simplify(function)
    results = {}
    while True:
        middle = (a + b) / 2
        error = abs(a - b) / 2
        results[iteration] = {
            "a": a,
            "b": b,
            "x": middle,
            "error": error
        }
        iteration += 1
        function_middle = simplified.subs(X, middle)
        if function_middle == 0 or error <= abs(tolerance):
            break
        elif bolzano(simplified, a, middle):
            b = middle
        else:
            a = middle
    return results


def newton_raphson(function, a, b, tolerance):
    if not bolzano(function, a, b):
        raise Exception(
            f"invalid initial interval -> bolzano does not apply in this interval [{min(a, b)};{max(a, b)}]")

    def select(expression: Expr, derivative: Derivative, a, b):
        function_a = expression.subs(X, a)
        derivative_a = derivative.subs(X, a)
        return a if function_a * derivative_a > 0 else b

    simplified = simplify(function)
    derivative_x = diff(simplified, X)
    results = {}
    n = select(simplified, derivative_x, a, b)
    iteration = 1
    while True:
        function_n = simplified.subs(X, n)
        derivative_n = derivative_x.subs(X, n)
        m = float(n - function_n / derivative_n)
        error = abs(m - n) if iteration > 1 else (b - a) / 2
        results[iteration] = {
            "a": a,
            "b": b,
            "x": float(m),
            "error": error
        }
        iteration += 1
        if error <= tolerance:
            break
        n = m
    return results


def regula_falsi(function: Expr, a, b, tolerance):
    if not bolzano(function, a, b):
        raise Exception(
            f"invalid initial interval -> bolzano does not apply in this interval [{min(a, b)};{max(a, b)}]")
    results = {}
    iteration = 1
    simplified = simplify(function)
    m = 0
    while True:
        function_a = float(simplified.subs(X, a))
        function_b = float(simplified.subs(X, b))
        f_difference = function_b - function_a
        v_difference = b - a
        n = a - (v_difference / f_difference) * function_a
        error = abs(n - m) if iteration > 1 else abs((b - a) / 2)
        results[iteration] = {
            "a": float(a),
            "b": float(b),
            "x": float(n),
            "error": float(error)
        }
        iteration += 1
        if error < tolerance:
            break
        elif bolzano(function, a, n):
            b = n
        else:
            a = n
        m = n
    return results


def secant(function: Expr, a, b, tolerance):
    if not bolzano(function, a, b):
        raise Exception(
            f"invalid initial interval -> bolzano does not apply in this interval [{min(a, b)};{max(a, b)}]")
    results = {}
    iteration = 1
    simplified = simplify(function)
    function_a = simplified.subs(X, a)
    function_b = simplified.subs(X, b)

    while True:
        n = b - ((b - a) / (function_b - function_a) * function_b)
        function_c = simplified.subs(X, n)
        error = abs(function_c - b)
        results[iteration] = {
            "a": float(a),
            "b": float(b),
            "x": float(function_c),
            "error": float(error),
        }
        if error <= tolerance:
            break
        # a -> b; b -> c; forget a
        a = b
        function_a = function_b
        b = function_c
        function_b = function_c
        iteration += 1
    return results
