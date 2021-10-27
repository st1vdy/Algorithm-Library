namespace Dijkstra {
#define ll long long
    static constexpr ll INF = 1e18;
    int n, m; // 点数 边数
    struct edge {
        int to; // 点
        ll val; // 边权
        edge(int to_ = 0, ll val_ = 0) :to(to_), val(val_) {}
        bool operator < (const edge& k) const { return val > k.val; }
    };
    vector<vector<edge>> g;
    void init() { // 建图操作需要根据题意修改
        cin >> n >> m;
        g.resize(n);
        for (int i = 0; i < m; i++) {
            int u, v, w;
            cin >> u >> v >> w;
            --u, --v;
            g[u].push_back(edge(v, w));
        }
    }
    ll dijkstra(int s, int t) { // 最短路
        vector<ll> dis(n, INF);
        vector<bool> vis(n, false);
        dis[s] = 0;
        priority_queue<edge> pq;
        pq.push(edge(s, 0));
        while (!pq.empty()) {
            auto top = pq.top();
            pq.pop();
            if (!vis[top.to]) {
                vis[top.to] = true;
                for (auto nxt : g[top.to]) {
                    if (!vis[nxt.to] && dis[nxt.to] > nxt.val + dis[top.to]) {
                        dis[nxt.to] = nxt.val + dis[top.to];
                        pq.push(edge(nxt.to, dis[nxt.to]));
                    }
                }
            }
        }
        return dis[t];
    }
#undef ll
}