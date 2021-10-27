// UOJ117 返回欧拉回路的边集（负数代表走了反向边）
#include <bits/stdc++.h>
using namespace std;
// 有向图欧拉回路 任意点的入度=出度
vector<int> directed_euler_circuit(int n, int m, const vector<vector<pair<int, int>>>& g) {
    vector<int> d(n);
    for (const auto& A : g) {
        for (auto p : A) {
            d[p.first]++;
        }
    }
    for (int i = 0; i < n; i++) {
        if (g[i].size() != d[i]) {
            return {};
        }
    }
    vector<vector<pair<int, int>>::const_iterator> it(n);
    for (int i = 0; i < n; i++) it[i] = g[i].begin();
    vector<int> vis(m + 1), p;
    function<void(int)> dfs = [&](int u) {
        for (auto& nxt = it[u]; nxt != g[u].end();) {
            if (!vis[nxt->second]) {
                vis[nxt->second] = 1;
                int v = nxt->second;
                dfs(nxt->first);
                p.push_back(v);
            } else {
                nxt = next(nxt);
            }
        }
    };
    for (int i = 0; i < n; i++) {
        if (!g[i].empty()) {
            dfs(i);
            break;
        }
    }
    if (p.size() < m) return {};
    reverse(p.begin(), p.end());
    return p;
}
// 无向图欧拉回路 任意点的度数为偶数
vector<int> undirected_euler_circuit(int n, int m, const vector<vector<pair<int, int>>>& g) {
    for (const auto& A : g) {
        if (A.size() & 1) {
            return {};
        }
    }
    vector<vector<pair<int, int>>::const_iterator> it(n);
    for (int i = 0; i < n; i++) it[i] = g[i].begin();
    vector<int> vis(m + 1), p;
    function<void(int)> dfs = [&](int u) {
        for (auto& nxt = it[u]; nxt != g[u].end();) {
            if (!vis[abs(nxt->second)]) {
                vis[abs(nxt->second)] = 1;
                int v = nxt->second;
                dfs(nxt->first);
                p.push_back(v);
            } else {
                nxt = next(nxt);
            }
        }
    };
    for (int i = 0; i < n; i++) {
        if (!g[i].empty()) {
            dfs(i);
            break;
        }
    }
    if (p.size() < m) return {};
    reverse(p.begin(), p.end());
    return p;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);
    int t, n, m;
    cin >> t >> n >> m; // t=1是无向图 t=-1是有向图
    vector<vector<pair<int, int>>> G(n + 1);
    for (int i = 1, u, v; i <= m; i++) {
        cin >> u >> v;
        if (t == 1) G[v].push_back({ u, -i });
        G[u].push_back({ v, i });
    }
    auto p = t == 1 ? undirected_euler_circuit(n + 1, m, G) : directed_euler_circuit(n + 1, m, G);
    if (p.size() == m) {
        cout << "YES\n";
        for (int x : p)
            cout << x << " \n"[x == p.back()];
    } else {
        cout << "NO\n";
    }
    return 0;
}