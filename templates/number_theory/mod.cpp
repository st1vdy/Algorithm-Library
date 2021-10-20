namespace SimpleMod {
    constexpr int MOD = (int)1e9 + 7;
    inline int norm(long long a) { return (a % MOD + MOD) % MOD; }
    inline int add(int a, int b) { return a + b >= MOD ? a + b - MOD : a + b; }
    inline int sub(int a, int b) { return a - b < 0 ? a - b + MOD : a - b; }
    inline int mul(int a, int b) { return (int)((long long)a * b % MOD); }
    inline int powmod(int a, long long b) {
        int res = 1;
        while (b > 0) {
            if (b & 1) res = mul(res, a);
            a = mul(a, a);
            b >>= 1;
        }
        return res;
    }
    inline int inv(int a) {
        a %= MOD;
        if (a < 0) a += MOD;
        int b = MOD, u = 0, v = 1;
        while (a) {
            int t = b / a;
            b -= t * a; swap(a, b);
            u -= t * v; swap(u, v);
        }
        assert(b == 1);
        if (u < 0) u += MOD;
        return u;
    }
}