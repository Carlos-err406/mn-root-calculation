from sympy import simplify, diff, Expr, Derivative, symbols
from src.utils import bolzano

x, y = symbols('x y')
MAX_ITERATIONS = 1000


def bisection(function: Expr, a, b, tolerance):
    if not bolzano(function, a, b):
        print(f"invalid initial interval -> bolzano does not apply in this interval [{min(a, b)};{max(a, b)}]")
        return None
    iteration = 1
    simplified = simplify(function)
    results = {}
    while iteration <= MAX_ITERATIONS:
        # Compute the midpoint of the interval
        middle = (a + b) / 2
        error = abs(a - b) / 2
        results[iteration] = {
            "a": a,
            "b": b,
            "x": middle,
            "error": error
        }
        iteration += 1
        # Check if the tolerance has been reached
        if error <= abs(tolerance):
            break
        elif bolzano(simplified, a, middle):  # Update the interval
            b = middle
        else:
            a = middle
    return results


def newton_raphson(function: Expr, a, b, tolerance):
    if not bolzano(function, a, b):
        print(f"invalid initial interval -> bolzano does not apply in this interval [{min(a, b)};{max(a, b)}]")
        return None

    def select(expression: Expr, derivative: Derivative, a, b):
        function_a = expression.subs(x, a)
        derivative_a = derivative.subs(x, a)
        return a if function_a * derivative_a > 0 else b

    simplified = simplify(function)
    derivative_x = diff(simplified, x)
    results = {}
    # Determine initial guess x0
    n = select(simplified, derivative_x, a, b)
    iteration = 1
    while iteration <= MAX_ITERATIONS:
        function_n = simplified.subs(x, n)
        derivative_n = derivative_x.subs(x, n)
        m = float(n - function_n / derivative_n)
        error = abs(m - n) if iteration > 1 else None
        results[iteration] = {
            "a": a,
            "b": b,
            "x": float(m),
            "error": error
        }
        if iteration > 1 and (error <= tolerance):
            break
        n = m
        iteration += 1
    return results


def regula_falsi(function: Expr, a, b, tolerance):
    if not bolzano(function, a, b):
        print(f"invalid initial interval -> bolzano does not apply in this interval [{min(a, b)};{max(a, b)}]")
        return None
    results = {}
    iteration = 1
    simplified = simplify(function)
    m = 0
    while iteration <= MAX_ITERATIONS:
        function_a = float(simplified.subs(x, a))
        function_b = float(simplified.subs(x, b))
        n = float(a - ((b - a) / (function_b - function_a)) * function_a)
        error = float(abs(n - m)) if iteration > 1 else None
        results[iteration] = {
            "a": float(a),
            "b": float(b),
            "x": float(n),
            "error": error
        }
        function_n = simplified.subs(x, n)
        if function_n == 0 or (iteration > 1 and error < tolerance):
            break
        elif bolzano(function, a, n):
            b = n
        else:
            a = n
        m = n
        iteration += 1
    return results


def secant(function: Expr, a, b, tolerance):
    if not bolzano(function, a, b):
        print(f"invalid initial interval -> bolzano does not apply in this interval [{min(a, b)};{max(a, b)}]")
        return None
    results = {}
    iteration = 1
    simplified = simplify(function)
    function_a = simplified.subs(x, a)
    function_b = simplified.subs(x, b)

    while iteration <= MAX_ITERATIONS:
        n = float(b - ((b - a) / (function_b - function_a)) * function_b)
        error = abs(n - b)
        results[iteration] = {
            "a": float(a),
            "b": float(b),
            "x": float(n),
            "error": float(error),
        }
        if error <= tolerance:
            break
        function_n = simplified.subs(x, n)
        # a -> b; b -> n; forget a
        a = b
        function_a = function_b
        b = n
        function_b = function_n
        iteration += 1
    return results
