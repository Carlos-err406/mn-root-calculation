import json
from sympy import simplify
from sympy.abc import x as X


def bolzano(function, a, b) -> bool:
    """Apply Bolzano's theorem to a given function and interval

    Keyword arguments:

    ``function`` -- a function represented as a string

    ``a,b`` -- the values that represent the intervals
    Return: True if the theorem applies, False otherwise
    """
    simplified = simplify(function)
    function_a = simplified.subs(X, a)  # f(a)
    function_b = simplified.subs(X, b)  # f(b)
    return (function_a * function_b) < 0
