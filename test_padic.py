from math import log, ceil
from sympy import nextprime
from hypothesis import given, note, assume
from hypothesis.strategies import integers, composite
from padic import Padic

zero = lambda p: Padic.from_int(0, p)
one = lambda p: Padic.from_int(1, p)
two = lambda p: Padic.from_int(2, p)
three = lambda p: Padic.from_int(3, p)
four = lambda p: Padic.from_int(4, p)
five = lambda p: Padic.from_int(5, p)
six = lambda p: Padic.from_int(6, p)
seven = lambda p: Padic.from_int(7, p)
eight = lambda p: Padic.from_int(8, p)
nine = lambda p: Padic.from_int(9, p)
ten = lambda p: Padic.from_int(10, p)
# pytest -s --hypothesis-show-statistics


@composite
def primes(draw, min_value=1, max_value=10**9):
    p = draw(integers(min_value=min_value, max_value=max_value))
    return nextprime(p)


@composite
def padics(draw):
    #N = draw(integers(min_value=1, max_value=10**4))
    a = draw(integers())
    p = draw(primes())
    return Padic.from_int(a, p)#, N if N < 100 else ceil(N * log(2, p)))


@composite
def _padics(draw):
    #N = draw(integers(min_value=1, max_value=10**4))
    a = draw(integers())
    return lambda p: Padic.from_int(a, p)#, N if N < 100 else ceil(N * log(2, p)))


@given(primes())
def test_from_int_1(p):
    assert 0 == zero(p)
    assert zero(p) == 0
    for i in range(1, 15):
        assert i == Padic.from_int(i, p, 100)


@given(primes())
def test_add_0(p):
    assert zero(p) + zero(p) == zero(p)
    assert one(p) + one(p) == two(p)
    assert two(p) + two(p) == four(p)
    assert one(p) + two(p) == three(p)
    assert two(p) + three(p) == five(p)
    assert three(p) + seven(p) == ten(p)


@given(padics())
def test_add_1(x):
    note(str(x))
    assert x + 0 == x
    assert 0 + x == x


@given(_padics(), _padics(), primes())
def test_add_2(_x, _y, p):
    x, y = _x(p), _y(p)
    note(str(x) + " " + str(y))
    assert x + y == y + x


@given(_padics(), _padics(), _padics(), primes())
def test_add_3(_x, _y, _z, p):
    x, y, z = _x(p), _y(p), _z(p)
    note(str(x) + " " + str(y) + " " + str(z))
    assert (x + y) + z == x + (y + z)


@given(primes())
def test_neg_sub_0(p):
    assert two(p) - one(p) == one(p)
    assert three(p) - two(p) == one(p)
    assert five(p) - three(p) == two(p)
    assert ten(p) - five(p) == five(p)
    assert seven(p) - two(p) == five(p)
    assert six(p) - zero(p) == six(p)
    assert four(p) - four(p) == zero(p)
    assert -zero(p) == zero(p)


@given(padics())
def test_neg_sub_1(x):
    note(str(x))
    assert -(-x) == x
    assert x - x == 0
    assert 0 == x - x
    assert 0 - x == -x
    assert x - 0 == x


@given(primes())
def test_mul_0(p):
    assert zero(p) * one(p) == zero(p)
    assert one(p) * two(p) == two(p)
    assert two(p) * three(p) == six(p)
    assert four(p) * two(p) == eight(p)
    assert two(p) * five(p) == ten(p)
    assert (-two(p)) * five(p) == -ten(p)
    assert three(p) * three(p) == nine(p)
    assert one(p) * one(p) == one(p)


@given(padics())
def test_mul_1(x):
    note(str(x))
    assert 1 * x == x
    assert -1 * x == -x
    assert 0 * x == 0


@given(_padics(), _padics(), primes())
def test_mul_2(_x, _y, p):
    x, y = _x(p), _y(p)
    note(str(x) + " " + str(y))
    assert x * y == y * x
    assert (-x) * y == -(x * y)


@given(_padics(), _padics(), _padics(), primes())
def test_mul_3(_x, _y, _z, p):
    x, y, z = _x(p), _y(p), _z(p)
    note(str(x) + " " + str(y) + " " + str(z))
    assert x * (y * z) == (x * y) * z
    assert x * (y + z) == x * y + x * z
    assert (x + y) * z == x * z + y * z


@given(primes())
def test_div_0(p):
    assert three(p) / one(p) == three(p)
    assert six(p) / two(p) == three(p)
    assert ten(p) / five(p) == two(p)
    assert ten(p) / two(p) == five(p)
    assert nine(p) / three(p) == three(p)
    assert four(p) / two(p) == two(p)
    assert eight(p) / (-four(p)) == -two(p)
    assert -six(p) / two(p) == -three(p)
    assert -eight(p) / (-four(p)) == two(p)
    assert zero(p) / ten(p) == zero(p)


@given(padics())
def test_div_1(x):
    note(str(x))
    assume(x != 0)
    assert x * (1/x) == 1
    assert 1/(x*x) == (1/x) * (1/x)
    assert x/1 == x
    assert 0/x == 0
    assert x/x == 1
    assert 1/(1/x) == x


@given(_padics(), _padics(), primes())
def test_div_2(_x, _y, p):
    x, y = _x(p), _y(p)
    note(str(x) + " " + str(y))
    assume(y != 0)
    assert (x/y) * y == x
    assert x * (1/y) == x/y
    assert (-x)/y == -(x/y)


@given(primes())
def test_mod_0(p):
    assert ten(p) % two(p) == zero(p)
    assert five(p) % three(p) == two(p) % three(p)
    assert seven(p) % one(p) == zero(p)
    assert four(p) % three(p) == one(p) % three(p)
    assert eight(p) % five(p) == three(p) % five(p)
    assert six(p) % zero(p) == six(p)
    assert six(p) % zero(p) != one(p) or p == 5


@given(padics())
def test_mod_1(x):
    assert x % x == 0


# FAILED test_padic.py::test_div_1 - assert (1 / (110000000 + O(2^64) * 110000000 + O(2^64))) == ((1 / 110000000 + O(2^64)) * (1 / 110000000 + O(2^64)))
# Zdaje się, że zachodzą problemy:
# a) wypisywanie na stringa dla 1/x * 1/x daje więcej cyfr niż powinno przez co nie wszystkie wskazane cyfry są poprawne
# b) dokładność obliczeń jest niższa niż oczekiwana przy ujemnej waluacji wyniku - N patrzy na liczbę cyfr
# zamiast na faktyczną dokładność - czyli przy 12 cyfrach po przecinku program błędnie traktuje 43 cyfrowe przybliżenie
# jako O(2^43) podczas gdy faktycznie daje to jedynie O(2^31)
