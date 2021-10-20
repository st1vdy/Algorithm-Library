/*
 * 最大团 Bron-Kerbosch algorithm
 * 最劣复杂度 O(3^(n/3))
 * 采用位运算形式实现
 * */
namespace Max_clique {
#define ll long long
#define TWOL(x) (1ll<<(x))
    const int N = 60;
    int n, m;      // 点数 边数
    int r = 0;     // 最大团大小
    ll G[N];       // 以二进制形式存图
    ll clique = 0; // 最大团 以二进制形式存储
    void BronK(int S, ll P, ll X, ll R) { // 调用时参数这样设置: 0, TWOL(n)-1, 0, 0
        if (P == 0 && X == 0) {
            if (r < S) {
                r = S;
                clique = R;
            }
        }
        if (P == 0) return;
        int u = __builtin_ctzll(P | X);
        ll c = P & ~G[u];
        while (c) {
            int v = __builtin_ctzll(c);
            ll pv = TWOL(v);
            BronK(S + 1, P & G[v], X & G[v], R | pv);
            P ^= pv; X |= pv; c ^= pv;
        }
    }
    void init() {
        cin >> n >> m;
        for (int i = 0; i < m; i++) {
            int u, v;
            cin >> u >> v;
            --u, --v;
            G[u] |= TWOL(v);
            G[v] |= TWOL(u);
        }
        BronK(0, TWOL(n)-1, 0, 0);
        cout << r << ' ' << clique << '\n';
    }
}