#pragma once

#include <cstdlib>
#include <iostream>

unsigned GCD(unsigned u, unsigned v);


class Frac;
bool operator==(const Frac& a, const Frac& b);
bool operator!=(const Frac& a, const Frac& b);
Frac operator+(const Frac& a, const Frac& b);
Frac operator-(const Frac&a);
Frac operator-(const Frac& a, const Frac& b);
Frac operator*(const Frac& a, const Frac& b);
Frac operator/(const Frac& a, const Frac& b);


using namespace std;
ostream& operator << (ostream& out, const Frac& a);
istream& operator >> (istream& in, Frac& a);

class Frac {
private:
    int _m;
    unsigned int _n;
    
    void simplify(int& m, unsigned& n) {
        unsigned gcd = GCD((unsigned) std::abs(m), n);
        m /= (int)gcd;
        n /= gcd;
    }
    
public:
    Frac(int m, unsigned n): _m(m), _n(n) {
        if(n == 0) {
            throw "Divide by zero";
        }
        simplify(_m, _n);
    }
    
    Frac(int m) : _m(m), _n(1) {};
    
    Frac() : _m(0), _n(1) {};
    
    int m() const { return _m;};
    unsigned n() const { return _n;};
    
    Frac operator=(const Frac& a) {
        _m = a.m();
        _n = a.n();
        return *this;
    }
    
    Frac operator/=(const Frac& a) {
        if(a.m() == 0) {
            throw "Divide by zero";
        }
        if(a.m() < 0) {
            _m = -_m;
        }
        _m *= a.n();
        _n *= std::abs(a.m());
        simplify(_m, _n);
        return *this;
    }
    
    
    Frac operator +=(const Frac& a) {
        _m *= a.n();
        _m += a.m() * _n;
        _n *= a.n();
        simplify(_m, _n);
        return *this;
    }
    
    Frac operator -=(const Frac& a) {
        *this += (-a);
        return *this;
    }
};


