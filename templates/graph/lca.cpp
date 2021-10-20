constexpr int SIZE = 200010;
constexpr int DEPTH = 21; // 最大深度 2^DEPTH - 1
int pa[SIZE][DEPTH], dep[SIZE];
vector<int> g[SIZE]; //邻接表
void dfs(int rt, int fin) { //预处理深度和祖先
    pa[rt][0] = fin;
    dep[rt] = dep[pa[rt][0]] + 1; //深度
    for (int i = 1; i < DEPTH; i++) { // rt 的 2^i 祖先等价于 rt 的 2^(i-1) 祖先 的 2^(i-1) 祖先
        pa[rt][i] = pa[pa[rt][i - 1]][i - 1];
    }
    int sz = g[rt].size();
    for (int i = 0; i < sz; ++i) {
        if (g[rt][i] == fin) continue;
        dfs(g[rt][i], rt);
    }
}

int LCA(int x, int y) {
    if (dep[x] > dep[y]) swap(x, y);
    int dif = dep[y] - dep[x];
    for (int j = 0; dif; ++j, dif >>= 1) {
        if (dif & 1) {
            y = pa[y][j]; //先跳到同一高度
        }
    }
    if (y == x) return x;
    for (int j = DEPTH - 1; j >= 0 && y != x; j--) { //从底往上跳
        if (pa[x][j] != pa[y][j]) { //如果当前祖先不相等 我们就需要更新
            x = pa[x][j];
            y = pa[y][j];
        }
    }
    return pa[x][0];
}