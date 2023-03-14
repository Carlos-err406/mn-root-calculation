import json
from math import e, pi

from sympy import init_printing
from sympy import ln, sin
from sympy.abc import x

from src.methods import bisection, newton_raphson, regula_falsi, secant

init_printing()
methods = {
    "bisection": bisection,
    'newton-raphson': newton_raphson,
    'regula-falsi': regula_falsi,
    'secant': secant
}
ex1 = e ** -x - ln(x)
ex2 = e ** -x - sin(x)
ex3 = e ** -x ** 2 - x ** 2
ex4 = x * e ** x - 2
expressions = [ex1, ex2, ex3, ex4]

tolerance = 0.01
a = 0.1
b = 2

for expr in expressions:
    for name, method in methods.items():
        if expr == ex1:
            a = 0.1
            b = 2
        elif expr == ex2:
            a = 0.1
            b = pi / 3
        print(
            f"method: {name}\nexpression: f(x) = {expr}\ntolerance: {tolerance}\ninterval: [{min(a, b)}; {max(a, b)}]")
        result = method(expr, a, b, tolerance)
        print(json.dumps(result, indent=2))
        print("---------" * 8)
