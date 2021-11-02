// n,a,b,c>=1e9 && k1+k2<=10
// O(k^4log(MAXN))
using ll = long long;
const int K = 10;
const int mod = 1e9 + 7;
const ll lmod = ll(mod) << 32;
const int signs[2] = { 1, mod - 1 };
int invs[K + 2], binom[K + 2][K + 2], B[K + 1], polys[K + 1][K + 2];

int add(int a, int b) { return (a += b - mod) < 0 ? a + mod : a; }
ll add64(ll a, ll b) { return (a += b - lmod) < 0 ? a + lmod : a;}
int mul(int a, int b) { return 1LL * a * b % mod; }
int power_sum(int e, int x) {
    int ret = 0;
    for (int i = 0; i < e + 2; ++i)
        ret = add(mul(ret, x), polys[e][i]);
    return mul(ret, invs[e + 1]);
}

void init() {
    invs[0] = invs[1] = 1; B[0] = 1;
    for (int i = 2; i <= K + 1; ++i)
        invs[i] = mul(invs[mod % i], mod - mod / i);
    for (int i = 0; i <= K + 1; ++i) {
        binom[i][0] = 1;
        for (int j = 1; j <= i; ++j)
            binom[i][j] = add(binom[i - 1][j - 1], binom[i - 1][j]);
    }    
    for (int i = 1; i <= K; ++i) {
        int s = 0;
        for (int j = 0; j < i; ++j)
            s = add(s, mul(binom[i + 1][j], B[j]));
        B[i] = mul(mul(s, invs[i + 1]), signs[1]);
    }
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= i; ++j)
            polys[i][j] = mul(mul(binom[i + 1][j], B[j]), signs[j & 1]);
        polys[i][i + 1] = 0;
    }
    polys[0][1] = 1;
}

int euclidLike(int N, int a, int b, int c, int k1, int k2) {
    assert(N >= 0); assert(a >= 0); assert(b >= 0); assert(c >= 1); assert(k1 + k2 <= K);
    using T = tuple<int, int, int, int>;
    stack<T> stac;
    while (1) {
        stac.emplace(N, a, b, c);
        if (N < 0 || a == 0)
            break;
        if (a >= c) {
            a %= c;
        } else if (b >= c) {
            b %= c;
        } else {
            N = (ll(a) * N + b) / c - 1;
            b = c - 1 - b;
            swap(a, c);
        }
    }

    const int S = k1 + k2;
    static int curr[K + 1][K + 1] = {}, next[K + 1][K + 1] = {};
    while (!stac.empty()) {
        tie(N, a, b, c) = stac.top();
        stac.pop();
        if (N < 0) {
            ;
        } else if (a == 0) {
            int q = b / c;
            for (int e1 = 0; e1 <= S; ++e1) {
                int s = power_sum(e1, N);
                for (int e2 = 0; e2 <= S - e1; ++e2)
                    next[e1][e2] = s, s = mul(s, q);
            }
        } else if (a >= c || b >= c) {
            int q = (a >= c) ? a / c : b / c;
            int d = (a >= c) ? 1 : 0;
            for (int e1 = 0; e1 <= S; ++e1) {
                for (int e2 = 0; e2 <= S - e1; ++e2) {
                    ll s = 0;
                    int p = 1;
                    for (int i2 = 0; i2 <= e2; ++i2) {
                        s = add64(s, ll(p) * mul(binom[e2][i2], curr[e1 + i2 * d][e2 - i2]));
                        p = mul(p, q);
                    }
                    next[e1][e2] = s % mod;
                }
            }
        } else {
            static int cumu[K + 1][K + 1];
            for (int e2 = 0; e2 <= S - 1; ++e2) {
                for (int e1 = 0; e1 <= S - e2 - 1; ++e1) {
                    ll s = 0;
                    for (int j = 0; j <= e1 + 1; ++j) {
                        s = add64(s, ll(polys[e1][e1 + 1 - j]) * curr[e2][j]);
                    }
                    cumu[e1][e2] = mul(s % mod, invs[e1 + 1]);
                }
            }
            const int M = (ll(a) * N + b) / c;
            for (int e1 = 0; e1 <= S; ++e1) {
                int p = power_sum(e1, N);
                for (int e2 = 0; e2 <= S - e1; ++e2) {
                    ll t = 0;
                    for (int i2 = 0; i2 < e2; ++i2) {
                        t = add64(t, ll(cumu[e1][i2]) * binom[e2][i2]);
                    }
                    next[e1][e2] = add(p, mod - t % mod);
                    p = mul(p, M);
                }
            }
        }
        swap(curr, next);
    }
    return curr[k1][k2];
}