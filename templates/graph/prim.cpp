/*
 * Prim Algorithm O(N^2 + M)
 * g是一个带边权的邻接表 下标从0开始
 * 如果图不连通则返回1
 */
long long prim(vector<vector<pair<int, int>>>& g) {
    int n = g.size();
    long long res = 0;
    vector<int> dis(n, INT_MAX), vis(n);
    dis[0] = 0;
    for (int i = 0; i < n; i++) {
        int mn = INT_MAX, p = -1;
        for (int j = 0; j < n; j++) {
            if (!vis[j] && dis[j] < mn) {
                mn = dis[j];
                p = j;
            }
        }
        if (p == -1) {
            return -1;
        }
        vis[p] = true;
        res += dis[p];
        for (auto& [j, k] : g[p]) {
            if (!vis[j] && k < dis[j]) {
                dis[j] = k;
            }
        }
    }
    return res;
}