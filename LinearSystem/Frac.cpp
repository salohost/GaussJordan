#include "Frac.hpp"

unsigned GCD(unsigned u, unsigned v) {
    while ( v != 0) {
        unsigned r = u % v;
        u = v;
        v = r;
    }
    return u;
}

bool operator==(const Frac& a, const Frac& b) {
    return a.m() == b.m() && a.n() == b.n();
}

bool operator!=(const Frac& a, const Frac& b) {
    return !(a == b);
}

Frac operator+(const Frac& a, const Frac& b) {
    return Frac(a.m()*b.n() + b.m()*a.n(), a.n()*b.n());
}

Frac operator-(const Frac&a) {
    return Frac(-a.m(), a.n());
}

Frac operator-(const Frac& a, const Frac& b) {
    return a + (-b);
}

Frac operator*(const Frac& a, const Frac& b) {
    if(a.m() == 0 || b.m() == 0) return Frac(0);
    else return Frac(a.m() * b.m(), a.n() * b.n());
}

Frac operator/(const Frac& a, const Frac& b) {
    if(a.m() == 0 && b.m() != 0) return Frac(0);
    if(b.m() < 0) return Frac(-a.m() * b.n(), -a.n() * b.m());
    return Frac(a.m()*b.n(), a.n() * b.m());
}

ostream& operator << (ostream& out, const Frac& a) {
    if (a.n() == 1) return out << a.m();
    return out << "\\frac{" << a.m() << "}{" << a.n() << "}";
}

istream& operator >> (istream& in, Frac& a) {
    int m; unsigned n;
    in >> m >> n;
    a = Frac(m, n);
    return in;
}

