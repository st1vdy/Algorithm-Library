constexpr int N = 601;
constexpr int inf = 0x3f3f3f3f;
int edge[N][N]; // 边权存这里
int dis[N], vis[N], bin[N];
int n, m;
int contract(int& s, int& t) {  // Find s,t
    memset(dis, 0, sizeof(dis));
    memset(vis, false, sizeof(vis));
    int i, j, k, mincut, maxc;
    for (i = 1; i <= n; i++) {
        k = -1;
        maxc = -1;
        for (j = 1; j <= n; j++) {
            if (!bin[j] && !vis[j] && dis[j] > maxc) {
                k = j;
                maxc = dis[j];
            }
        }
        if (k == -1) return mincut;
        s = t; t = k;
        mincut = maxc;
        vis[k] = true;
        for (j = 1; j <= n; j++) {
            if (!bin[j] && !vis[j]) {
                dis[j] += edge[k][j];
            }
        }
    }
    return mincut;
}

int stoerWagner() { // O(NM + N^2logN) <=> O(N^3)
    int mincut, i, j, s, t, ans;
    for (mincut = inf, i = 1; i < n; i++) {
        ans = contract(s, t);
        bin[t] = true;
        if (mincut > ans) mincut = ans;
        if (mincut == 0) return 0;
        for (j = 1; j <= n; j++) {
            if (!bin[j]) {
                edge[s][j] = (edge[j][s] += edge[j][t]);
            }
        }
    }
    return mincut;
}