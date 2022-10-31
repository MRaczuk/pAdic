import pytest
from math import log, ceil
from sympy import isprime, nextprime
from hypothesis import given, settings, Verbosity, note, assume
from hypothesis.strategies import integers, composite
from padic import Padic

zero = lambda p: Padic.from_int(100, 0, p)
# choosen_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
#                   37, 97, 101, 103, 107, 109, 997, 2137,
#                   10**9 + 7]


@composite
def primes(draw, min_value=2, max_value=10**9):
    p = draw(integers(min_value=min_value, max_value=max_value))
    return nextprime(p)


@composite
def padics(draw):
    N = draw(integers(min_value=1, max_value=10**4))
    a = draw(integers())
    p = draw(primes())
    return Padic.from_int(N if N < 100 else ceil(N * log(2, p)), a, p)


@composite
def _padics(draw):
    N = draw(integers(min_value=1, max_value=10**4))
    a = draw(integers())
    return lambda p: Padic.from_int(N if N < 100 else ceil(N * log(2, p)), a, p)


@given(primes())
def test_from_int(p):
    assert 0 == zero(p)
    assert zero(p) == 0
    for i in range(1, 15):
        assert i == Padic.from_int(100, i, p)


@given(padics())
def test_add1(n):
    note(str(n))
    assert n - n == 0
    assert n + n - n == n


@given(_padics(), _padics(), _padics(), primes())
def test_add2(_x, _y, _z, p):
    x, y, z = _x(p), _y(p), _z(p)
    note(str(x) + " " + str(y) + " " + str(z))
    assert (x + y) + z == x + (y + z)
