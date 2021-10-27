namespace Berlekamp_Massey {
    typedef long long ll;
    constexpr ll MOD = 1e9 + 7;
    constexpr int N = 10010;
    ll res[N], base[N], _c[N], _md[N];
    vector<int> Md;
    ll powmod(ll a, ll b) {
        ll res = 1;
        while (b > 0) {
            if (b & 1) res = res * a % MOD;
            a = a * a % MOD;
            b >>= 1;
        }
        return res;
    }
    void mul(ll* a, ll* b, int k) {
        for (int i = 0; i < k + k; i++)
            _c[i] = 0;
        for (int i = 0; i < k; i++) {
            if (!a[i]) continue;
            for (int j = 0; j < k; j++) {
                _c[i + j] = (_c[i + j] + a[i] * b[j]) % MOD;
            }
        }
        for (int i = k + k - 1; i >= k; i--) {
            if (!_c[i]) continue;
            for (int j = 0; j < Md.size(); j++) {
                _c[i - k + Md[j]] = (_c[i - k + Md[j]] - _c[i] * _md[Md[j]]) % MOD;
            }
        }
        for (int i = 0; i < k; i++)
            a[i] = _c[i];
    }
    int solve(ll n, vector<int> a, vector<int> b) { //a系数 b初值 b[n+1]=a[0]*b[n]+...
        //printf("%d\n", (int)b.size());
        //for (int i = 0; i < b.size(); i++)
        //    printf("b[%d] = %d\n", i, b[i]);
        //printf("%d\n", (int)a.size());
        //printf("b[n]");
        //for (int i = 0; i < a.size(); i++) {
        //    if (!i)putchar('='); else putchar('+');
        //    printf("%d*b[n-%d]", a[i], i + 1);
        //}
        //puts("");
        ll ans = 0, pnt = 0;
        int k = a.size();
        for (int i = 0; i < k; i++) {
            _md[k - 1 - i] = -a[i];
        }
        _md[k] = 1;
        Md.clear();
        for (int i = 0; i < k; i++) {
            if (_md[i]) {
                Md.push_back(i);
            }
            res[i] = base[i] = 0;
        }
        res[0] = 1;
        while ((1LL << pnt) <= n) pnt++;
        for (int p = pnt; p >= 0; p--) {
            mul(res, res, k);
            if ((n >> p) & 1) {
                for (int i = k - 1; i >= 0; i--)
                    res[i + 1] = res[i];
                res[0] = 0;
                for (int j = 0; j < Md.size(); j++) {
                    res[Md[j]] = (res[Md[j]] - res[k] * _md[Md[j]]) % MOD;
                }
            }
        }
        for (int i = 0; i < k; i++)
            ans = (ans + res[i] * b[i]) % MOD;
        return (ans < 0 ? ans + MOD : ans);
    }
    vector<int> BM(vector<int> s) { // O(n^2)
        vector<int> C(1, 1), B(1, 1);
        int L = 0, m = 1, b = 1;
        for (int n = 0; n < (int)s.size(); n++) {
            ll d = 0;
            for (int i = 0; i <= L; i++)
                d = (d + (ll)C[i] * s[n - i]) % MOD;
            if (!d) {
                ++m;
            } else if (2 * L <= n) {
                auto T = C;
                ll c = MOD - d * powmod(b, MOD - 2) % MOD;
                while (C.size() < B.size() + m)
                    C.push_back(0);
                for (int i = 0; i < B.size(); i++)
                    C[i + m] = (C[i + m] + c * B[i]) % MOD;
                L = n + 1 - L; B = T; b = d; m = 1;
            } else {
                ll c = MOD - d * powmod(b, MOD - 2) % MOD;
                while (C.size() < B.size() + m) C.push_back(0);
                for (int i = 0; i < B.size(); i++) {
                    C[i + m] = (C[i + m] + c * B[i]) % MOD;
                }
                ++m;
            }
        }
        return C;
    }
    int work(vector<int> a, ll n) { // 这里的n不是数组大小 是求数列第n项的值
        vector<int> c = BM(a);      // 求第n项的复杂度为 O(k^2logn) k是递推数列大小
        c.erase(c.begin());
        for (int i = 0; i < c.size(); i++)
            c[i] = (MOD - c[i]) % MOD;
        return solve(n, c, vector<int>(a.begin(), a.begin() + (int)c.size()));
    }
}