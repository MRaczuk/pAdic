from numpy import base_repr


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
            return Padic(min(self.v + other.N, other.s + self.N), self.center() * other.center(), self.p)
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
        raise RuntimeError("Valuation undefined for " + str(n))

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