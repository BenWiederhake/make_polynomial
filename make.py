#!/usr/bin/env python3
"""
Makes a polynomial through arbitrary points with arbitrary derivatives, inspired by the Lagrange polynomial.
Write your requirements in the REQUIREMENTS global, and run the program to compute the exact polynomial.
"""


import sympy

LN2 = sympy.log(2)
S2 = sympy.sqrt(2)
X = sympy.Symbol("x")

REQUIREMENTS = [
    # (place, value, value-of-derivative-or-None),
    # Example: (42, 2, 1337) specifies that we desire f(42)=2 and f'(42)=1337
    # Example: (7, 1, None) specifies that we desire f(7)=1
    # The resulting function will be of degree d where d is the degree of requirements,
    # i.e. len(REQUIREMENTS)+sum(d is not None for _,_,d in REQUIREMENTS).
    (1, 0, 1 / LN2),
    (S2, sympy.Rational(1, 2), None),
    (2, 1, 1 / (2 * LN2)),
]


def compute_basis(i, requirements):
    """
    Let (x_i, y_i, y'_i) be the i-th entry in 'requirements'.
    
    This function returns a polynomial p such that:
    - p(x_i) = y_i
    - If y'_i is not None, then p'(x_i) = y'_i
    - For all j != i, it holds that p(x_j) = 0
    - For all j != i, if y'_j is not None, then it p'(x_j) = 0
    """
    x_i, y_i, yp_i = requirements[i]
    # We multiply polynomials that are one in the place of interest, and zero in all the "unrelated" places.
    akku = sympy.sympify(1)
    for j, entry in enumerate(requirements):
        if j == i:
            continue
        x_j, _y_j, yp_j = entry
        akku *= (X - x_j) / (x_i - x_j)
        if yp_j is None:
            continue
        # TODO: What about requirements that *only* affect the derivative, but not the value?
        akku *= (X - x_j) / (x_i - x_j)
    # Now 'akku' has value one at x_i, and an indeterminate derivative.
    assert akku.subs(X, x_i) == 1, (akku, i, requirements)
    if yp_i is None:
        akku *= y_i
        assert akku.subs(X, x_i) == y_i, (akku, i, requirements)
    else:
        old_derivative = sympy.diff(akku, X).subs(X, x_i)
        normalization = X * (yp_i - old_derivative * y_i)
        normalization += y_i - normalization.subs(X, x_i)
        akku *= normalization
        # Now 'akku' has value y_i at x_i, and the derivative yp_i.
        assert akku.subs(X, x_i) == y_i, akku
        assert sympy.diff(akku, X).subs(X, x_i) == yp_i, (akku, i, requirements)
    return akku


def compute_polynomial(requirements):
    places = [x for x, _, _ in requirements]
    assert len(places) == len(set(places)), f"Duplicate places in requirements!"
    akku = sympy.sympify(0)
    # We add polynomials that are zero in all the "unrelated" places, and have the right values in the place of interest.
    for i, e in enumerate(requirements):
        akku += compute_basis(i, requirements)
    return akku


def run():
    formula = compute_polynomial(REQUIREMENTS)
    print(formula)


if __name__ == "__main__":
    run()
