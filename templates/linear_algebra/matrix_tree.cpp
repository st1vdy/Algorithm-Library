/*
 * 矩阵树定理
 * 有向图：若 u->v 有一条权值为 w 的边 基尔霍夫矩阵 a[v][v] += w, a[v][u] -= w
 * 生成树数量为除去 根所在行和列 后的n-1阶行列式的值
 * 无向图：若 u->v 有一条权值为 w 的边 基尔霍夫矩阵 a[v][v] += w, a[v][u] -= w, a[u][u] += w, a[u][v] -= w
 * 生成树数量为除去 任意一行和列 后的n-1阶行列式的值
 * 无权图则边权默认为1
 * */
typedef long long ll;
typedef unsigned long long u64;
int a[SIZE][SIZE];
int gauss(int a[][SIZE], int n) { // 任意模数求行列式 O(n^2(n + log(mod)))
    int ans = 1;
    for (int i = 1; i <= n; i++) {
        int* x = 0, * y = 0;
        for (int j = i; j <= n; j++) {
            if (a[j][i] && (x == NULL || a[j][i] < x[i])) {
                x = a[j];
            }
        }
        if (x == 0) {
            return 0;
        }
        for (int j = i; j <= n; j++) {
            if (a[j] != x && a[j][i]) {
                y = a[j];
                for (;;) {
                    int v = md - y[i] / x[i], k = i;
                    for (; k + 3 <= n; k += 4) {
                        y[k + 0] = (y[k + 0] + u64(x[k + 0]) * v) % md;
                        y[k + 1] = (y[k + 1] + u64(x[k + 1]) * v) % md;
                        y[k + 2] = (y[k + 2] + u64(x[k + 2]) * v) % md;
                        y[k + 3] = (y[k + 3] + u64(x[k + 3]) * v) % md;
                    }
                    for (; k <= n; ++k) {
                        y[k] = (y[k] + u64(x[k]) * v) % md;
                    }
                    if (!y[i]) break;
                    swap(x, y);
                }
            }
        }
        if (x != a[i]) {
            for (int k = i; k <= n; k++) {
                swap(x[k], a[i][k]);
            }
            ans = md - ans;
        }
        ans = 1LL * ans * a[i][i] % md;
    }
    return ans;
}