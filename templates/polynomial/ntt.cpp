constexpr int MOD = 998244353;
struct Z {
    int val;
    Z(long long v = 0) {
        if (v < 0) v = v % MOD + MOD;
        if (v >= MOD) v %= MOD;
        val = v;
    }
    static int mod_inv(int a, int m = MOD) {
        // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Example
        int g = m, r = a, x = 0, y = 1;
        while (r != 0) {
            int q = g / r;
            g %= r; swap(g, r);
            x -= q * y; swap(x, y);
        }
        return x < 0 ? x + m : x;
    }
    explicit operator int() const { return val; }
    explicit operator uint64_t() const { return val; }
    Z& operator+=(const Z& other) {
        val += other.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
    Z& operator-=(const Z& other) {
        val -= other.val;
        if (val < 0) val += MOD;
        return *this;
    }
    static unsigned fast_mod(uint64_t x, unsigned m = MOD) {
        return x % m;
        /*
        // Optimized mod for Codeforces 32-bit machines.
        // x must be less than 2^32 * m for this to work, so that x / m fits in a 32-bit integer.
        unsigned x_high = x >> 32, x_low = (unsigned) x;
        unsigned quot, rem;
        asm("divl %4\n"
            : "=a" (quot), "=d" (rem)
            : "d" (x_high), "a" (x_low), "r" (m));
        return rem;*/
    }
    Z& operator*=(const Z& other) {
        val = fast_mod((uint64_t)val * other.val);
        return *this;
    }
    Z& operator/=(const Z& other) { return *this *= other.inv(); }
    friend Z operator+(const Z& a, const Z& b) { return Z(a) += b; }
    friend Z operator-(const Z& a, const Z& b) { return Z(a) -= b; }
    friend Z operator*(const Z& a, const Z& b) { return Z(a) *= b; }
    friend Z operator/(const Z& a, const Z& b) { return Z(a) /= b; }
    Z& operator++() {
        val = val == MOD - 1 ? 0 : val + 1;
        return *this;
    }
    Z& operator--() {
        val = val == 0 ? MOD - 1 : val - 1;
        return *this;
    }
    Z operator++(int) { Z before = *this; ++* this; return before; }
    Z operator--(int) { Z before = *this; --* this; return before; }
    Z operator-() const { return val == 0 ? 0 : MOD - val; }
    bool operator==(const Z& other) const { return val == other.val; }
    bool operator!=(const Z& other) const { return val != other.val; }
    Z inv() const { return mod_inv(val); }
    Z pow(long long p) const {
        assert(p >= 0);
        Z a = *this, res = 1;
        while (p > 0) {
            if (p & 1) res *= a;
            a *= a;
            p >>= 1;
        }
        return res;
    }
    friend ostream& operator<<(ostream& stream, const Z& m) {
        return stream << m.val;
    }
};
vector<Z> inv, fac, ifac;
void prepare_factorials(int maximum) {
    if (maximum + 1 < inv.size()) return;
    inv.assign(maximum + 1, 1);

    // Make sure MOD is prime, which is necessary for the inverse algorithm below.
    for (int p = 2; p * p <= MOD; p++)
        assert(MOD % p != 0);
    for (int i = 2; i <= maximum; i++)
        inv[i] = inv[MOD % i] * (MOD - MOD / i);

    fac.resize(maximum + 1);
    ifac.resize(maximum + 1);
    fac[0] = ifac[0] = 1;

    for (int i = 1; i <= maximum; i++) {
        fac[i] = i * fac[i - 1];
        ifac[i] = inv[i] * ifac[i - 1];
    }
}
inline Z binom(int n, int m) {
    assert(n < fac.size());
    if (n < m || n < 0 || m < 0) return 0;
    return fac[n] * ifac[m] * ifac[n - m];
}
namespace NTT {
    using ull = uint64_t;
    vector<Z> roots = { 0, 1 };
    vector<int> bit_reverse;
    int max_size = -1;
    Z root;

    bool is_power_of_two(int n) { return (n & (n - 1)) == 0; }
    int round_up_power_two(int n) {
        assert(n > 0);
        while (n & (n - 1))
            n = (n | (n - 1)) + 1;
        return n;
    }
    // Given n (a power of two), finds k such that n == 1 << k.
    int get_length(int n) {
        assert(is_power_of_two(n));
        return __builtin_ctz(n);
    }
    // Rearranges the indices to be sorted by lowest bit first, then second lowest, etc., rather than highest bit first.
    // This makes even-odd div-conquer much easier.
    void bit_reorder(int n, vector<Z>& values) {
        if ((int)bit_reverse.size() != n) {
            bit_reverse.assign(n, 0);
            int length = get_length(n);

            for (int i = 0; i < n; i++)
                bit_reverse[i] = (bit_reverse[i >> 1] >> 1) + ((i & 1) << (length - 1));
        }
        for (int i = 0; i < n; i++)
            if (i < bit_reverse[i])
                swap(values[i], values[bit_reverse[i]]);
    }
    /* 找原根 998244353的原根是31 */
    void find_root() {
        int order = MOD - 1;
        max_size = 1;
        while (order % 2 == 0) {
            order /= 2;
            max_size *= 2;
        }
        root = 2;
        // Find a max_size-th primitive root of MOD.
        while (!((int)root.pow(max_size) == 1 && (int)root.pow(max_size / 2) != 1))
            root = root + 1;
    }

    void prepare_roots(int n) {
        if (max_size < 0) find_root();
        assert(n <= max_size);
        if ((int)roots.size() >= n) return;
        int len = get_length(roots.size());
        roots.resize(n);

        // The roots array is set up such that for a given power of two n >= 2, roots[n / 2] through roots[n - 1] are
        // the first half of the n-th primitive roots of MOD.
        while (1 << len < n) {
            // z is a 2^(length + 1)-th primitive root of MOD.
            Z z = root.pow(max_size >> (len + 1));
            for (int i = 1 << (len - 1); i < 1 << len; i++) {
                roots[2 * i] = roots[i];
                roots[2 * i + 1] = roots[i] * z;
            }
            len++;
        }
    }
    void fft_iterative(int N, vector<Z>& values) {
        assert(is_power_of_two(N));
        prepare_roots(N);
        bit_reorder(N, values);
        for (int n = 1; n < N; n *= 2)
            for (int start = 0; start < N; start += 2 * n)
                for (int i = 0; i < n; i++) {
                    Z even = values[start + i];
                    Z odd = values[start + n + i] * roots[n + i];
                    values[start + n + i] = even - odd;
                    values[start + i] = even + odd;
                }
    }

    constexpr int FFT_CUTOFF = 150;
    vector<Z> mod_multiply(vector<Z> left, vector<Z> right) {
        int n = left.size(), m = right.size();
        // Brute force when either n or m is small enough.
        if (min(n, m) < FFT_CUTOFF) {
            const ull ULL_BOUND = numeric_limits<ull>::max() - (ull)MOD * MOD;
            vector<ull> res(n + m - 1);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    res[i + j] += (ull)((int)left[i]) * ((int)right[j]);
                    if (res[i + j] > ULL_BOUND)
                        res[i + j] %= MOD;
                }
            }
            for (ull& x : res)
                if (x >= MOD)
                    x %= MOD;
            return vector<Z>(res.begin(), res.end());
        }

        int N = round_up_power_two(n + m - 1);
        left.resize(N);
        right.resize(N);
        fft_iterative(N, left);
        fft_iterative(N, right);
        Z inv_N = Z(N).inv();
        for (int i = 0; i < N; i++)
            left[i] *= right[i] * inv_N;
        reverse(left.begin() + 1, left.end());
        fft_iterative(N, left);
        left.resize(n + m - 1);
        return left;
    }
    vector<Z> sub_conv(vector<Z> a, vector<Z> b) {
        int n = b.size();
        reverse(b.begin(), b.end());
        auto res = mod_multiply(a, b);
        return vector<Z>(res.begin() + n - 1, res.end());
    }
    /* 
     * 多项式平移c个单位
     * a(x+c)=\sum x^j 1/j!\sum i!a_i c^{j-i}/(j-i)!
     * O(NlogN)
     */
    vector<Z> shift(vector<Z> a, Z c) {
        int N = round_up_power_two(a.size());
        prepare_factorials(N);
        vector<Z>tmp(a.size());
        Z pc = 1;
        for (int i = 0; i < tmp.size(); i++, pc *= c) {
            tmp[i] = pc * ifac[i];
            a[i] *= fac[i];
        }
        tmp = sub_conv(a, tmp);
        for (int i = 0; i < a.size(); i++) {
            tmp[i] *= ifac[i];
        }
        return tmp;
    }
    vector<Z> mod_power(const vector<Z>& v, int exponent) {
        assert(exponent >= 0);
        int nn = v.size();
        vector<Z> res = { 1 };
        if (exponent == 0)
            return res;
        for (int k = 31 - __builtin_clz(exponent); k >= 0; k--) {
            res = mod_multiply(res, res);
            if (exponent >> k & 1)
                res = mod_multiply(res, v);
            res.resize(nn);
        }
        return res;
    }
    vector<Z> mod_inv(vector<Z> a) {
        int n = a.size();
        int N = round_up_power_two(a.size());
        a.resize(N * 2);
        vector<Z>res = { a[0].inv() };
        for (int i = 2; i <= N; i <<= 1) {
            vector<Z>tmp(a.begin(), a.begin() + i);
            int n = (i << 1);
            res.resize(n);
            tmp.resize(n);
            fft_iterative(n, tmp);
            fft_iterative(n, res);
            Z inv_n = Z(n).inv();
            for (int j = 0; j < n; j++)
                res[j] = res[j] * (2 - tmp[j] * res[j]) * inv_n;
            reverse(res.begin() + 1, res.end());
            fft_iterative(n, res);
            fill(res.begin() + i, res.end(), 0);
        }
        res.resize(n);
        return res;
    }
    vector<Z> integral(vector<Z> a) {
        assert(a.size() <= inv.size());
        a.push_back(0);
        for (int i = (int)a.size() - 1; i >= 1; i--) {
            a[i] = a[i - 1] * inv[i];
        }
        return a;
    }
    vector<Z> differential(vector<Z> a) {
        assert(a.size());
        for (int i = 0; i < (int)a.size() - 1; i++) {
            a[i] = a[i + 1] * (i + 1);
        }
        a.pop_back();
        return a;
    }
    vector<Z> ln(vector<Z>a) {
        assert((int)a[0] == 1);
        auto b = mod_multiply(differential(a), mod_inv(a));
        b = integral(b);
        b[0] = 0;
        return b;
    }
    vector<Z> exp(vector<Z>a) {
        int N = round_up_power_two(a.size());
        int n = a.size();
        a.resize(N * 2);
        vector<Z> res{ 1 };
        for (int i = 2; i <= N; i <<= 1) {
            auto tmp = res;
            tmp.resize(i);
            tmp = ln(tmp);
            for (int j = 0; j < i; j++)tmp[j] = a[j] - tmp[j];
            tmp[0] += 1;
            res.resize(i);
            res = mod_multiply(res, tmp);
            fill(res.begin() + i, res.end(), 0);
        }
        res.resize(n);
        return res;
    }
    vector<Z> poly_div(vector<Z> a, vector<Z> b) {
        if (a.size() < b.size()) {
            return { 0 };
        }
        reverse(a.begin(), a.end());
        reverse(b.begin(), b.end());
        b.resize(a.size() - b.size() + 1);
        a.resize(b.size());
        auto d = mod_multiply(mod_inv(b), a);
        d.resize(b.size());
        reverse(d.begin(), d.end());
        return d;
    }
    /* poly A mod B */
    vector<Z> poly_mod(vector<Z> a, vector<Z> b) {
        auto res = mod_multiply(b, poly_div(a, b));
        a.resize(b.size() - 1);
        res.resize(b.size() - 1);
        for (int i = 0; i < a.size(); i++) {
            a[i] -= res[i];
        }
        return a;
    }
    /* Brute force calc f(k) */
    Z eval(const vector<Z>& f, int k) {
        Z ans;
        for (int i = f.size() - 1; i >= 0; i--)
            ans = (ans * k + f[i]);
        return ans;
    }
    vector<Z> mod_multiply_all(const vector<vector<Z>>& polynomials) {
        if (polynomials.empty())
            return { 1 };
        struct compare_size {
            bool operator()(const vector<Z>& x, const vector<Z>& y) {
                return x.size() > y.size();
            }
        };
        priority_queue<vector<Z>, vector<vector<Z>>, compare_size> pq(polynomials.begin(), polynomials.end());
        while (pq.size() > 1) {
            vector<Z> a = pq.top(); pq.pop();
            vector<Z> b = pq.top(); pq.pop();
            pq.push(mod_multiply(a, b));
        }
        return pq.top();
    }
    /*
     * 线性递推求第K项
     * a[n] = \sum seq[i] * a[n-i], calc a[k]
     * O(Nlog(N)log(K))
     */
    Z linear_seq(const vector<Z>& _init, vector<Z> seq, long long k) {
        reverse(seq.begin(), seq.end());
        for (auto& i : seq)i = -i; seq.push_back(1);
        vector<Z> b = seq;
        reverse(b.begin(), b.end());
        b = mod_inv(b);
        auto poly_mod = [&](vector<Z>a) {
            if (a.size() < seq.size()) {
                a.resize(seq.size() - 1);
                return a;
            }
            vector<Z>tmp = a;
            reverse(a.begin(), a.end());
            b.resize(a.size() - seq.size() + 1);
            a.resize(b.size());
            auto d = mod_multiply(b, a);
            d.resize(b.size());
            reverse(d.begin(), d.end());
            auto res = mod_multiply(seq, d);
            tmp.resize(seq.size() - 1);
            res.resize(seq.size() - 1);
            for (int i = 0; i < tmp.size(); i++) {
                tmp[i] -= res[i];
            }
            return tmp;
        };
        vector<Z> a{ 0,1 };
        vector<Z> res{ 1 };
        for (; k; k >>= 1) {
            if (k & 1)
                res = poly_mod(mod_multiply(res, a));
            a = poly_mod(mod_multiply(a, a));
        }
        Z ans = 0;
        for (int i = 0; i < _init.size(); i++) {
            ans += _init[i] * res[i];
        }
        return ans;
    }
    /*
     * 多点求值
     * 注意如果是1..n的多点求值不要用这个
     * O(Nlog(N)log(M)) 常数偏大
     */
    vector<Z> multi_eval(vector<Z> F, vector<Z> x) {
        vector<vector<Z>> base;
        function<void(int, int, int)> build = [&](int l, int r, int o) {
            if (r - l == 1) {
                base[o] = { 1,-x[l] };
                return;
            }
            int mid = (l + r) >> 1;
            build(l, mid, o << 1);
            build(mid, r, o << 1 | 1);
            base[o] = (mod_multiply(base[o << 1], base[o << 1 | 1]));
        };
        vector<Z> res(x.size());
        int n = max(x.size(), F.size());
        x.resize(n); F.resize(n + 1);
        base.resize(4 * n);
        build(0, n, 1);
        function<void(vector<Z>, int, int, int)>solve = [&](vector<Z> f, int l, int r, int o) {
            if (r - l == 1) {
                if (l < res.size())res[l] = f[0];
                return;
            }
            int mid = (l + r) >> 1;
            auto L = sub_conv(f, base[o << 1 | 1]);
            auto R = sub_conv(f, base[o << 1]);
            L.resize(mid - l);
            R.resize(r - mid);
            solve(L, l, mid, o << 1);
            solve(R, mid, r, o << 1 | 1);
        };
        solve(sub_conv(F, mod_inv(base[1])), 0, n, 1);
        return res;
    }
    /*
     * 多点插值 给定N个点(xi, yi) 插出一个N-1阶的多项式
     * O(Nlog(N)log(M)) 常数偏大
     */
    vector<Z> multi_inter(const vector<Z>& y, const vector<Z>& x) {
        assert(y.size() == x.size());
        vector<vector<Z>>base;
        function<void(int, int, int)> build = [&](int l, int r, int o) {
            if (r - l == 1) {
                base[o] = { -x[l],1 };
                return;
            }
            int mid = (l + r) >> 1;
            build(l, mid, o << 1);
            build(mid, r, o << 1 | 1);
            base[o] = (mod_multiply(base[o << 1], base[o << 1 | 1]));
        };
        int n = x.size();
        base.resize(4 * n);
        int s = clock();
        build(0, n, 1);
        vector<Z> pi = differential(base[1]);
        vector<Z> res = multi_eval(pi, x);
        for (int i = 0; i < n; i++)res[i] = y[i] / res[i];
        function<vector<Z>(int, int, int)>solve2 = [&](int l, int r, int o) {
            if (r - l == 1) {
                return vector<Z>({ res[l] });
            }
            int mid = (l + r) >> 1;
            vector<Z> L = mod_multiply(solve2(l, mid, o << 1), base[o << 1 | 1]);
            vector<Z> R = mod_multiply(solve2(mid, r, o << 1 | 1), base[o << 1]);
            int n = max(L.size(), R.size());
            L.resize(n); R.resize(n);
            for (int i = 0; i < n; i++) {
                L[i] += R[i];
            }
            return L;
        };
        auto ans = solve2(0, n, 1);
        ans.resize(n);
        return ans;
    }
    vector<Z> mod_power2(const vector<Z>& a, int exponent) {
        int n = a.size();
        auto res = ln(a);
        while ((int)res.size() > n) res.pop_back();
        for (auto& i : res) {
            i = exponent * i;
        }
        res = exp(res);
        return res;
    }
}
