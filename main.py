import json

from sympy.abc import x, y, symbols
from src.methods import bisection, newton_raphson, regula_falsi, secant

methods = {
    "bisection": bisection,
    'newton-raphson': newton_raphson,
    'regula-falsi': regula_falsi,
    'secants': secant
}

expression = x
tolerance = 0.05
a = -10
b = 10

for name, method in methods.items():
    print(
        f"method: {name}\nexpression: f(x) = {expression}\ntolerance: {tolerance}\ninterval: [{min(a, b)};{max(a, b)}]")
    try:
        result = method(expression, a, b, tolerance)
        print(json.dumps(result, indent=2))
    except Exception as e:
        print(e)

    print("\n\n")
