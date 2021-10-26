#define db double
#ifndef ONLINE_JUDGE // 这三个函数是给MSVC用的，G++不需要
inline int __builtin_clz(int v) { // 返回前导0的个数
    return __lzcnt(v);
}
inline int __builtin_ctz(int v) { // 返回末尾0的个数
    if (v == 0) {
        return 0;
    }
    __asm {
        bsf eax, dword ptr[v];
    }
}
inline int __builtin_popcount(int v) { // 返回二进制中1的个数
    return __popcnt(v);
}
#endif
struct Complex {
    db real, imag;
    Complex(db x = 0, db y = 0) :real(x), imag(y) {}
    Complex& operator+=(const Complex& rhs) {
        real += rhs.real; imag += rhs.imag;
        return *this;
    }
    Complex& operator-=(const Complex& rhs) {
        real -= rhs.real; imag -= rhs.imag;
        return *this;
    }
    Complex& operator*=(const Complex& rhs) {
        db t_real = real * rhs.real - imag * rhs.imag;
        imag = real * rhs.imag + imag * rhs.real;
        real = t_real;
        return *this;
    }
    Complex& operator/=(double x) {
        real /= x, imag /= x;
        return *this;
    }
    friend Complex operator + (const Complex& a, const Complex& b) { return Complex(a) += b; }
    friend Complex operator - (const Complex& a, const Complex& b) { return Complex(a) -= b; }
    friend Complex operator * (const Complex& a, const Complex& b) { return Complex(a) *= b; }
    friend Complex operator / (const Complex& a, const db& b) { return Complex(a) /= b; }
    Complex power(long long p) const {
        assert(p >= 0);
        Complex a = *this, res = { 1,0 };
        while (p > 0) {
            if (p & 1) res = res * a;
            a = a * a;
            p >>= 1;
        }
        return res;
    }
    static long long val(double x) { return x < 0 ? x - 0.5 : x + 0.5; }
    inline long long Real() const { return val(real); }
    inline long long Imag() const { return val(imag); }
    Complex conj()const { return Complex(real, -imag); }
    explicit operator int()const { return Real(); }
    friend ostream& operator<<(ostream& stream, const Complex& m) {
        return stream << complex<db>(m.real, m.imag);
    }
};
constexpr int MOD = 998244353;
constexpr int Phi_MOD = 998244352;
inline int exgcd(int a, int md = MOD) {
    a %= md;
    if (a < 0) a += md;
    int b = md, u = 0, v = 1;
    while (a) {
        int t = b / a;
        b -= t * a; swap(a, b);
        u -= t * v; swap(u, v);
    }
    assert(b == 1);
    if (u < 0) u += md;
    return u;
}
inline int add(int a, int b) { return a + b >= MOD ? a + b - MOD : a + b; }
inline int sub(int a, int b) { return a - b < 0 ? a - b + MOD : a - b; }
inline int mul(int a, int b) { return 1LL * a * b % MOD; }
inline int powmod(int a, long long b) {
    int res = 1;
    while (b > 0) {
        if (b & 1) res = mul(res, a);
        a = mul(a, a);
        b >>= 1;
    }
    return res;
}

vector<int> inv, fac, ifac;
void prepare_factorials(int maximum) {
    inv.assign(maximum + 1, 1);
    // Make sure MOD is prime, which is necessary for the inverse algorithm below.
    for (int p = 2; p * p <= MOD; p++)
        assert(MOD % p != 0);
    for (int i = 2; i <= maximum; i++)
        inv[i] = mul(inv[MOD % i], (MOD - MOD / i));

    fac.resize(maximum + 1);
    ifac.resize(maximum + 1);
    fac[0] = ifac[0] = 1;

    for (int i = 1; i <= maximum; i++) {
        fac[i] = mul(i, fac[i - 1]);
        ifac[i] = mul(inv[i], ifac[i - 1]);
    }
}
namespace FFT {
    vector<Complex> roots = { Complex(0, 0), Complex(1, 0) };
    vector<int> bit_reverse;
    int max_size = 1 << 20;
    const long double pi = acosl(-1.0l);
    constexpr int FFT_CUTOFF = 150;
    inline bool is_power_of_two(int n) { return (n & (n - 1)) == 0; }
    inline int round_up_power_two(int n) {
        assert(n > 0);
        while (n & (n - 1)) {
            n = (n | (n - 1)) + 1;
        }
        return n;
    }
    // Given n (a power of two), finds k such that n == 1 << k.
    inline int get_length(int n) {
        assert(is_power_of_two(n));
        return __builtin_ctz(n);
    }
    // Rearranges the indices to be sorted by lowest bit first, then second lowest, etc., rather than highest bit first.
    // This makes even-odd div-conquer much easier.
    void bit_reorder(int n, vector<Complex>& values) {
        if ((int)bit_reverse.size() != n) {
            bit_reverse.assign(n, 0);
            int length = get_length(n);
            for (int i = 0; i < n; i++) {
                bit_reverse[i] = (bit_reverse[i >> 1] >> 1) + ((i & 1) << (length - 1));
            }
        }
        for (int i = 0; i < n; i++) {
            if (i < bit_reverse[i]) {
                swap(values[i], values[bit_reverse[i]]);
            }
        }
    }
    void prepare_roots(int n) {
        assert(n <= max_size);
        if ((int)roots.size() >= n)
            return;
        int length = get_length(roots.size());
        roots.resize(n);
        // The roots array is set up such that for a given power of two n >= 2, roots[n / 2] through roots[n - 1] are
        // the first half of the n-th primitive roots of MOD.
        while (1 << length < n) {
            for (int i = 1 << (length - 1); i < 1 << length; i++) {
                roots[2 * i] = roots[i];
                long double angle = pi * (2 * i + 1) / (1 << length);
                roots[2 * i + 1] = Complex(-cos(angle), -sin(angle));
            }
            length++;
        }
    }
    void fft_iterative(int N, vector<Complex>& values) {
        assert(is_power_of_two(N));
        prepare_roots(N);
        bit_reorder(N, values);
        for (int n = 1; n < N; n *= 2) {
            for (int start = 0; start < N; start += 2 * n) {
                for (int i = 0; i < n; i++) {
                    Complex& even = values[start + i];
                    Complex odd = values[start + n + i] * roots[n + i];
                    values[start + n + i] = even - odd;
                    values[start + i] = even + odd;
                }
            }
        }
    }
    vector<long long> multiply(vector<int> a, vector<int> b) { // 普通FFT
        int n = a.size();
        int m = b.size();
        if (min(n, m) < FFT_CUTOFF) {
            vector<long long> res(n + m - 1);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    res[i + j] += 1LL * a[i] * b[j];
                }
            }
            return res;
        }
        int N = round_up_power_two(n + m - 1);
        vector<Complex> tmp(N);
        for (int i = 0; i < a.size(); i++) tmp[i].real = a[i];
        for (int i = 0; i < b.size(); i++) tmp[i].imag = b[i];
        fft_iterative(N, tmp);
        for (int i = 0; i < N; i++) tmp[i] = tmp[i] * tmp[i];
        reverse(tmp.begin() + 1, tmp.end());
        fft_iterative(N, tmp);
        vector<long long> res(n + m - 1);
        for (int i = 0; i < res.size(); i++) {
            res[i] = tmp[i].imag / 2 / N + 0.5;
        }
        return res;
    }
    vector<int> mod_multiply(vector<int> a, vector<int> b, int lim = max_size) { // 任意模数FFT
        int n = a.size();
        int m = b.size();
        if (min(n, m) < FFT_CUTOFF) {
            vector<int> res(n + m - 1);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    res[i + j] += 1LL * a[i] * b[j] % MOD;
                    res[i + j] %= MOD;
                }
            }
            return res;
        }
        int N = round_up_power_two(n + m - 1);
        N = min(N, lim);
        vector<Complex> P(N);
        vector<Complex> Q(N);
        for (int i = 0; i < n; i++) {
            P[i] = Complex(a[i] >> 15, a[i] & 0x7fff);
        }
        for (int i = 0; i < m; i++) {
            Q[i] = Complex(b[i] >> 15, b[i] & 0x7fff);
        }
        fft_iterative(N, P);
        fft_iterative(N, Q);
        vector<Complex>A(N), B(N), C(N), D(N);
        for (int i = 0; i < N; i++) {
            Complex P2 = P[(N - i) & (N - 1)].conj();
            A[i] = (P2 + P[i]) * Complex(0.5, 0),
                B[i] = (P2 - P[i]) * Complex(0, 0.5);
            Complex Q2 = Q[(N - i) & (N - 1)].conj();
            C[i] = (Q2 + Q[i]) * Complex(0.5, 0),
                D[i] = (Q2 - Q[i]) * Complex(0, 0.5);
        }
        for (int i = 0; i < N; i++) {
            P[i] = (A[i] * C[i]) + (B[i] * D[i]) * Complex(0, 1),
                Q[i] = (A[i] * D[i]) + (B[i] * C[i]) * Complex(0, 1);
        }
        reverse(P.begin() + 1, P.end());
        reverse(Q.begin() + 1, Q.end());
        fft_iterative(N, P);
        fft_iterative(N, Q);
        for (int i = 0; i < N; i++) {
            P[i] /= N, Q[i] /= N;
        }
        int size = min(n + m - 1, lim);
        vector<int> res(size);
        for (int i = 0; i < size; i++) {
            long long ac = P[i].Real() % MOD, bd = P[i].Imag() % MOD,
                ad = Q[i].Real() % MOD, bc = Q[i].Imag() % MOD;
            res[i] = ((ac << 30) + bd + ((ad + bc) << 15)) % MOD;
        }
        return res.resize(n + m - 1), res;
    }
    vector<int> mod_inv(vector<int> a) { // 多项式逆
        int n = a.size();
        int N = round_up_power_two(a.size());
        a.resize(N * 2);
        vector<int> res(1);
        res[0] = exgcd(a[0]);
        for (int i = 2; i <= N; i <<= 1) {
            vector<int> tmp(a.begin(), a.begin() + i);
            int n = (i << 1);
            tmp = mod_multiply(tmp, mod_multiply(res, res, n), n);
            res.resize(i);
            for (int j = 0; j < i; j++) {
                res[j] = add(res[j], sub(res[j], tmp[j]));
            }
        }
        res.resize(n);
        return res;
    }
    vector<int> integral(vector<int> a) { // 多项式积分
        assert(a.size() <= inv.size());
        a.push_back(0);
        for (int i = (int)a.size() - 1; i >= 1; i--) {
            a[i] = mul(a[i - 1], inv[i]);
        }
        return a;
    }
    vector<int> differential(vector<int> a) { // 多项式求导
        for (int i = 0; i < (int)a.size() - 1; i++) {
            a[i] = mul(i + 1, a[i + 1]);
        }
        a.pop_back();
        return a;
    }
    vector<int> ln(vector<int> a) { // 多项式对数函数
        assert((int)a[0] == 1);
        auto b = mod_multiply(differential(a), mod_inv(a));
        b = integral(b);
        b[0] = 0;
        return b;
    }
    vector<int> exp(vector<int> a) { // 多项式指数函数
        int N = round_up_power_two(a.size());
        int n = a.size();
        a.resize(N * 2);
        vector<int> res{ 1 };
        for (int i = 2; i <= N; i <<= 1) {
            auto tmp = res;
            tmp.resize(i);
            tmp = ln(tmp);
            for (int j = 0; j < i; j++) {
                tmp[j] = sub(a[j], tmp[j]);
            }
            tmp[0] = add(tmp[0], 1);
            res.resize(i);
            res = mod_multiply(res, tmp, i << 1);
            fill(res.begin() + i, res.end(), 0);
        }
        res.resize(n);
        return res;
    }
    // Multiplies many polynomials whose total degree is n in O(n log^2 n).
    vector<int> mod_multiply_all(const vector<vector<int>>& polynomials) {
        if (polynomials.empty())
            return { 1 };
        struct compare_size {
            bool operator()(const vector<int>& x, const vector<int>& y) {
                return x.size() > y.size();
            }
        };
        priority_queue<vector<int>, vector<vector<int>>, compare_size> pq(polynomials.begin(), polynomials.end());
        while (pq.size() > 1) {
            vector<int> a = pq.top(); pq.pop();
            vector<int> b = pq.top(); pq.pop();
            pq.push(mod_multiply(a, b));
        }
        return pq.top();
    }
    tuple<int, int, bool> power_reduction(string s, int n) { // 多项式快速幂预处理
        int p = 0, q = 0; bool zero = false;
        for (int i = 0; i < s.length(); i++) {
            p = mul(p, 10);
            p = add(p, s[i] - '0');
            q = 1LL * q * 10 % Phi_MOD; // Phi_MOD 是MOD的欧拉函数值
            q = (q + s[i] - '0');
            if (q >= Phi_MOD) q -= Phi_MOD;
            if (q >= (int)n) zero = true;
        }
        return { p,q,zero };
    }
    vector<int> power(vector<int> a, string s) { // 多项式快速幂 a^s O(nlogn)
        int n = a.size();
        auto [p, q, zero] = power_reduction(s, (int)a.size()); // 不需要降幂的话可以省去这部分
        if (a[0] == 1) {
            auto res = ln(a);
            while ((int)res.size() > n) res.pop_back();
            for (auto& i : res) {
                i = mul(p, i);
            }
            res = exp(res);
            return res;
        } else {
            int mn = -1;
            vector<int> copy_a;
            for (int i = 0; i < (int)a.size(); i++) {
                if (a[i]) {
                    mn = i;
                    break;
                }
            }
            if ((mn == -1) || (mn && (zero || (1LL * mn * p > n)))) { // a中所有元素都是0 或 偏移过大
                return vector<int>(n, 0);
            }
            int inverse_amin = exgcd(a[mn]);
            for (int i = mn; i < n; i++) {
                copy_a.emplace_back(mul(a[i], inverse_amin));
            }
            copy_a = ln(copy_a);
            while ((int)copy_a.size() > n) copy_a.pop_back();
            for (auto& i : copy_a) {
                i = mul(i, p);
            }
            copy_a = exp(copy_a);
            vector<int> res(n, 0);
            // shift是偏移量 power_k 是a_min^q(q是扩展欧拉定理降出来的幂次)
            int shift = mn * p, power_k = powmod(a[mn], q);
            for (int i = 0; i + shift < n; i++) {
                res[i + shift] = mul(copy_a[i], power_k);
            }
            return res;
        }
    }
    vector<long long> sub_convolution(vector<int> a, vector<int> b) { // 减法卷积 只保留非负次项
        int n = b.size();
        reverse(b.begin(), b.end());
        auto res = multiply(a, b);
        return vector<long long>(res.begin() + n - 1, res.end());
    }
    int bostan_mori(vector<int> p, vector<int> q, long long n) { // [x^n]p(x)/q(x)  O(2/3dlog(d)log(n+1)) d是多项式度数
        int i;
        for (; n; n >>= 1) {
            auto r = q;
            for (i = 1; i < r.size(); i += 2) {
                r[i] = MOD - r[i];
            }
            p = mod_multiply(p, r);
            q = mod_multiply(q, r);
            for (i = (n & 1); i < p.size(); i += 2) {
                p[i / 2] = p[i];
            }
            p.resize(i / 2);
            for (i = 0; i < q.size(); i += 2) {
                q[i / 2] = q[i];
            }
            q.resize(i / 2);
        }
        return p[0];
    }
};