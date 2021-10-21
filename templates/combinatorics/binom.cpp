namespace BinomialCoefficient {
    vector<int> fac, ifac, iv;
    // 组合数预处理 option=1则还会预处理线性逆元
    void prepareFactorials(int maximum = 1000000, int option = 0) {
        fac.assign(maximum + 1, 0);
        ifac.assign(maximum + 1, 0);
        fac[0] = ifac[0] = 1;
        if (option) { // O(3n)
            iv.assign(maximum + 1, 1);
            for (int p = 2; p * p <= MOD; p++)
                assert(MOD % p != 0);
            for (int i = 2; i <= maximum; i++)
                iv[i] = mul(iv[MOD % i], (MOD - MOD / i));
            for (int i = 1; i <= maximum; i++) {
                fac[i] = mul(i, fac[i - 1]);
                ifac[i] = mul(iv[i], ifac[i - 1]);
            }
        } else { // O(2n + log(MOD))
            for (int i = 1; i <= maximum; i++)
                fac[i] = mul(fac[i - 1], i);
            ifac[maximum] = inv(fac[maximum]);
            for (int i = maximum; i; i--)
                ifac[i - 1] = mul(ifac[i], i);
        }
    }
    inline int binom(int n, int m) {
        if (n < m || n < 0 || m < 0) return 0;
        return mul(fac[n], mul(ifac[m], ifac[n - m]));
    }
}