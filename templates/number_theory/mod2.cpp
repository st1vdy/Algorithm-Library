template<int MOD> struct Z {
    int x;
    Z(int v = 0) : x(v % MOD) { if (x < 0) x += MOD; }
    Z(long long v = 0) : x(v % MOD) { if (x < 0) x += MOD; }
    Z operator - () const { return x ? MOD - x : 0; }
    Z operator + (const Z& r) { return Z(*this) += r; }
    Z operator - (const Z& r) { return Z(*this) -= r; }
    Z operator * (const Z& r) { return Z(*this) *= r; }
    Z operator / (const Z& r) { return Z(*this) /= r; }
    Z& operator += (const Z& r) {
        x += r.x;
        if (x >= MOD) x -= MOD;
        return *this;
    }
    Z& operator -= (const Z& r) {
        x -= r.x;
        if (x < 0) x += MOD;
        return *this;
    }
    Z& operator *= (const Z& r) {
        x = 1LL * x * r.x % MOD;
        return *this;
    }
    Z& operator /= (const Z& r) {
        int a = r.x, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        x = 1LL * x * u % MOD;
        if (x < 0) x += MOD;
        return *this;
    }
    Z& power(long long k) {
        int a = x; x = 1;
        while (k > 0) {
            if (k & 1) x = 1LL * x * a % MOD;
            a = 1LL * a * a % MOD;
            k >>= 1;
        }
        return *this;
    }
    bool operator == (const Z& r) { return this->x == r.x; }
    bool operator != (const Z& r) { return this->x != r.x; }
    friend constexpr istream& operator >> (istream& is, Z<MOD>& x) noexcept {
        is >> x.x;
        x.x %= MOD;
        if (x.x < 0) x.x += MOD;
        return is;
    }
    friend ostream& operator << (ostream& os, const Z<MOD>& x) {
        return os << x.x;
    }
};
constexpr int MOD = 1e9 + 7;
using mint = Z<MOD>;