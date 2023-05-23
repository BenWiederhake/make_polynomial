#!/usr/bin/env python3

from make import compute_basis, compute_polynomial, LN2, S2, X
import sympy
import unittest


class BasisTests(unittest.TestCase):
    def run_tests(self, *test_specs):
        for i, (expected, *requirements) in enumerate(test_specs):
            with self.subTest(i=i, expected=expected):
                self.assertEqual(compute_basis(0, requirements), expected)

    def test_constant(self):
        self.run_tests(
            (0, (0, 0, None)),
            (0, (3, 0, None)),
            (1337, (42, 1337, None)),
        )

    # FIXME: Too lazy to write more unit tests


class PolynomialTests(unittest.TestCase):
    def run_tests(self, *test_specs):
        for i, (expected, *requirements) in enumerate(test_specs):
            with self.subTest(i=i, expected=expected):
                actual = compute_polynomial(requirements)
                print(actual)
                self.assertEqual(actual, expected)

    def test_constant(self):
        self.run_tests(
            (0, (0, 0, None)),
            (0, (3, 0, None)),
            (1337, (42, 1337, None)),
        )

    def test_linear_twoarg(self):
        self.run_tests(
            (-2 * X, (0, 0, None), (1, -2, None)),
            (0, (0, 0, None), (1, 0, None)),
            (-X + 3, (0, 3, None), (1, 2, None)),
            (sympy.Rational(-1, 2) * X + 8, (4, 6, None), (6, 5, None)),
        )

    def test_linear_onearg(self):
        self.run_tests(
            (X, (0, 0, 1)),
            (2 * X, (0, 0, 2)),
            (2 * X + 1, (0, 1, 2)),
            (2 * X + 1, (1, 3, 2)),
            (-3 * X + 21, (5, 6, -3)),
        )

    def test_cubic_hermite(self):
        self.run_tests(
            (0, (0, 0, 0), (1, 0, 0)),
            ((1 + 2 * X) * (1 - X) ** 2, (0, 1, 0), (1, 0, 0)),
            (X * (1 - X) ** 2, (0, 0, 1), (1, 0, 0)),
            (X * X * (3 - 2 * X), (0, 0, 0), (1, 1, 0)),
            (X * X * (X - 1), (0, 0, 0), (1, 0, 1)),
            # Cubic approximation of ln(x) on [1, 2]:
            ((2 - X) * (2 - X) * (X / LN2 - 1 / LN2) + (X - 1) * (X - 1) * (X * (-2 * LN2 + 1 / (2 * LN2)) - 1 / LN2 + 5 * LN2), (1, 0, 1 / LN2), (2, LN2, 1 / (2 * LN2))),
        )



if __name__ == "__main__":
    unittest.main()
