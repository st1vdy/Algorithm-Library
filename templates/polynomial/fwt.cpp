constexpr int MOD = 998244353;
constexpr int MAXN = 1 << 21;
using ll = long long;
namespace FWT {
    int a[MAXN], b[MAXN], ans[MAXN];
    int n;
    ll x, y;
    void print() {
        for (int i = 0; i < n; ++i) cout << a[i] << ' '; cout << "\n";
        for (int i = 0; i < n; ++i) cout << b[i] << ' '; cout << "\n";
        for (int i = 0; i < n; ++i) cout << ans[i] << ' '; cout << "\n";
    }
    void geta(int* t) {
        for (int i = 0; i < n; ++i) a[i] = t[i];
    }
    void getb(int* t) {
        for (int i = 0; i < n; ++i) b[i] = t[i];
    }
    void Or(int* a, int p) {
        for (int i = 1; i < n; i <<= 1) {
            for (int j = 0; j < n; j += (i << 1)) {
                for (int k = 0; k < i; ++k) {
                    x = a[j + k]; y = a[i + j + k];
                    a[j + k] = x;
                    a[i + j + k] = (x * p + y) % MOD;
                }
            }
        }
    }
    void workOr(int* A, int* B, int* Ans, int nn) {
        for (n = 1; n < nn; n <<= 1);
        geta(A); getb(B);
        Or(a, 1);
        Or(b, 1);
        for (int i = 0; i < n; ++i) ans[i] = 1LL * a[i] * b[i] % MOD;
        Or(ans, -1);
        for (int i = 0; i < n; ++i) Ans[i] = (ans[i] + MOD) % MOD;
    }
    void And(int* a, int p) {
        for (int i = 1; i < n; i <<= 1) {
            for (int j = 0; j < n; j += (i << 1)) {
                for (int k = 0; k < i; ++k) {
                    x = a[j + k]; y = a[i + j + k];
                    a[j + k] = (x + y * p) % MOD;
                    a[i + j + k] = y;
                }
            }
        }
    }
    void workAnd(int* A, int* B, int* Ans, int nn) {
        for (n = 1; n < nn; n <<= 1);
        geta(A); getb(B);
        And(a, 1);
        And(b, 1);
        for (int i = 0; i < n; ++i) ans[i] = 1LL * a[i] * b[i] % MOD;
        And(ans, -1);
        for (int i = 0; i < n; ++i) Ans[i] = (ans[i] + MOD) % MOD;
    }
    void Xor(int* a, int p) {
        for (int i = 1; i < n; i <<= 1) {
            for (int j = 0; j < n; j += (i << 1)) {
                for (int k = 0; k < i; ++k) {
                    x = a[j + k]; y = a[i + j + k];
                    a[j + k] = (x + y) % MOD;
                    a[i + j + k] = (x - y + MOD) % MOD;
                    if (p == -1) {
                        if (a[j + k] & 1) a[j + k] += MOD;
                        if (a[i + j + k] & 1) a[i + j + k] += MOD;
                        a[j + k] >>= 1;
                        a[i + j + k] >>= 1;
                    }
                }
            }
        }
    }
    void workXor(int* A, int* B, int* Ans, int nn) {
        for (n = 1; n < nn; n <<= 1);
        geta(A); getb(B);
        Xor(a, 1);
        Xor(b, 1);
        for (int i = 0; i < n; ++i) ans[i] = 1LL * a[i] * b[i] % MOD;
        Xor(ans, -1);
        for (int i = 0; i < n; ++i) Ans[i] = (ans[i] + MOD) % MOD;
    }
}