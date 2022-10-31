from numpy import base_repr
from math import log10, floor
import sys

sys.setrecursionlimit(2000)


class Padic:
    PRECISION = 30

    def __init__(self, N: int, v: int, s: int, p: int):
        # Assumes p is prime
        # Assumes s _|_ p or s == 0
        if v >= N:
            s = 0
        if s == 0:
            v = N
        self.N: int = N
        self.v: int = v
        self.s: int = s % p**(N - v)
        self.p: int = p

    def __abs__(self):
        return self.p**(-self.v)

    def __eq__(self, other):
        diff = self - other
        return diff.v >= Padic.PRECISION or diff.v == diff.N

    def __add__(self, other):
        if isinstance(other, Padic) and self.p == other.p:
            return Padic.from_int(min(self.N, other.N), self.center() + other.center(), self.p)
        if isinstance(other, int):
            return self + Padic.from_int(self.N, other, self.p)
        else:
            raise RuntimeError("Can't add " + str(self) + " to " + str(other))

    def __neg__(self):
        return Padic(self.N, self.v, -self.s, self.p)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, Padic) and self.p == other.p:
            return Padic.from_int(min(self.v + other.N, other.s + self.N), self.center() * other.center(), self.p)
        if isinstance(other, int):
            return self * Padic.from_int(self.N, other, self.p)
        else:
            raise RuntimeError("Can't multiply " + str(self) + " with " + str(other))

    def __truediv__(self, other):
        if isinstance(other, Padic) and self.p == other.p:
            N = min(self.v + other.N - 2 * other.v, self.N - other.v)
            v = self.v - other.v
            s = self.s * pow(other.s, -1, self.p**N)
            return Padic(N, v, s, self.p)
        if isinstance(other, int):
            return self / Padic.from_int(self.N, other, self.p)
        else:
            raise RuntimeError("Can't divide " + str(self) + " by " + str(other))

    def __str__(self):
        if self.p > 31:
            return f"{self.p}-adic + O({self.p}^{self.N})"
        out = base_repr(self.s, self.p)
        if self.v >= 0:
            return out + ''.join(['0']*self.v) + f' + O({self.p}^{self.N})'
        else:
            return out[:-self.v] + '.' + out[self.v:] + f' + O({self.p}^{self.N})'

    def __pow__(self, power, modulo=None):
        return NotImplemented

    def center(self):
        return self.s * (self.p ** self.v)

    @staticmethod
    def val(n, p):
        if isinstance(n, Padic) and n.p == p:
            return n.v
        if isinstance(n, int):
            if n == 0:
                return 10**9 + 1
            out = 0
            while n % p == 0:
                out += 1
                n //= p
            return out
        raise RuntimeError("Valuation undefined for " + str(n), type(n))

    @staticmethod
    def _digit_value(c: str):
        o = ord(c)
        if ord('0') <= o <= ord('9'):
            return o - ord('0')
        if ord('A') <= o <= ord('Z'):
            return o - ord('A') + 10
        raise RuntimeError("Couldn't assign digit value for: " + c)

    @staticmethod
    def from_string(string: str, p: int):
        v = 0
        N = len(string) + 1
        non_zero_occurred = False
        dot_occurred = False
        stop = 0
        for i in range(len(string) - 1, -1, -1):
            char = string[i]
            if char == '.' and not dot_occurred: #!@**(!*(WARUNEK, PD, PRACA
                v -= len(string) - 1 - i
                N = i + 1
                dot_occurred = True
                continue
            if not char.isalnum():
                raise RuntimeError("Cannot parse string: " + string + " to p-adic integer.")
            if not non_zero_occurred and char != '0':
                non_zero_occurred = True
                stop = i
                v += len(string) - 2 - i if dot_occurred else len(string) - 1 - i
            if Padic._digit_value(char) >= p:
                raise RuntimeError("Cannot parse string: " + string + " to p-adic integer.")
        s = 0
        for char in string[:stop + 1]:
            if char == '.':
                continue
            s *= p
            s += Padic._digit_value(char)
        if s == 0:
            return Padic(N, N, 0, p)
        return Padic(N, v, s, p)

    @staticmethod
    def from_int(N: int, a: int, p: int):
        if a == 0:
            return Padic(N, N, a, p)
        v = Padic.val(a, p)
        return Padic(N, v, a // (p ** v), p)

    @staticmethod
    def from_frac(N: int, a: int, b: int, p: int):
        return Padic.from_int(N, a, p) / Padic.from_int(N, b, p)


def gcd(a, b):
    if b == 0:
        return a
    return gcd(b, a % b)


class Rational:
    def __init__(self, p, q):
        g = gcd(p, q)
        self.num = p // g
        self.den = q // g

    def __add__(self, other):
        return Rational(self.num * other.den + self.den * other.num, self.den * other.den)

    def __mul__(self, other):
        return Rational(self.num * other.num, self.den * other.den)

    def __neg__(self):
        return Rational(-self.num, self.den)

    def __str__(self):
        return str(self.num) + "/" + str(self.den)


def series(a, n):  # To be changed
    def u(x):
        x = Rational(x, 1)
        out = Rational(0, 1)
        for k in range(n, -1, -1):
            out *= x
            out += a(k)
        return out
    return lambda x: u(x)


D = 2022
B = 2022
Padic.PRECISION = D
log = lambda x: -series(lambda n: Rational(1, n) if n != 0 else Rational(0, 1), B)(1-x)

# print(Padic.from_frac(D, 7, 9, 2))

# for B in range(1, 5):
#     log = lambda x: -series(lambda n: Rational(1, n) if n != 0 else Rational(0, 1), B)(1 - x)
#     a = log(-7)
#     print(B, a)
#     print(B, Padic.from_frac(D, a.num, a.den, 2))
b = log(9)
# c = log(9**873)
#print('a =', Padic.from_frac(D, a.num, a.den, 2))
#print('b =', Padic.from_frac(D, b.num, b.den, 2))
# print('c =', Padic.from_frac(D, c.num, c.den, 2))
for t in range(1, 1000, 2):
    a = log(-7 ** t)
    m = Padic.from_frac(D, a.num, a.den, 2)/Padic.from_frac(D, b.num, b.den, 2)
    print(m.center())
    print(floor(log10(m.center())), ((pow(9, m.center(), 2**2022) + 7**t) % 2**2022) == 0)

# 27656791 / 1101001100000001001010111
#            1101001100000001001010111

#8151772563618518648943966462536951830846266470431048222059177550454766167008929826580330113751948635750371679249228964047228313112386712764427710725611539027570970262641009887922388444251216507230471242220815426736884267321126113530362698791211779423552501769158800381255325932062580327491735668461947341939491759440720926409258629685901443058846113576890498327326279682462930726278095952869872802403720270025270520091579212693073172726098655890330077943469614120175775421420065275506344487673726779554922899460392081513937093064183421756298325057875213681291817590677960139018641807564639926920723976028759
#23200551212717227548968557797013065131823589055247993953759736430577449708331006004362337332799659616825426626965365436111354390756210951604493678197159095659131816199895381052172668104769335668618403560662416695812900169372185529169665436028387785371226961478030262317941316423730839032019739785364042787148634666679131871655573713518644354587109436602354728571410450543321111376186569814343605706405873173368795119408324211422673907340074932325476045403350029112386754848030730769023134749970024983297836122602502217821268825197751220004890868085420277128286503221659411786675294175520157019730529554399831
#9406539041756847046807308053134742191561476826476198934476973771241614959704713301958983122208134926176060090031550420077558861693335899453217164705287060451152064850668658499516885672235530360303481408220845011134636899912318924951785359136459332322983045598604939207079987304519482277947202888483746580609332371083751833981460805224961417648275017705207264737894668186530611528925813997135885503209007906732286960957992639349618783564319691235844266370528427368316346837649530333002243200724882134921855475239066108234479547058982467020302432146080577596880767141052428769399273054738402687952366349715205