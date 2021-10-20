/*
 * 高斯-约旦消元法
 * 可以修改为解异或方程组 修改策略为
 * a+b -> a^b
 * a-b -> a^b
 * a*b -> a&b
 * a/b -> a*(b==1)
 * */
constexpr double eps = 1e-7;
double a[SIZE][SIZE], ans[SIZE];
void gauss(int n) {
    vector<bool> vis(n, false);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (vis[j]) continue;
            if (fabs(a[j][i]) > eps) {
                vis[i] = true;
                for (int k = 0; k <= n; k++) swap(a[i][k], a[j][k]);
                break;
            }
        }
        if (fabs(a[i][i]) < eps) continue;
        for (int j = 0; j <= n; j++) {
            if (i != j && fabs(a[j][i]) > eps) {
                double res = a[j][i] / a[i][i];
                for (int k = 0; k <= n; k++) a[j][k] -= a[i][k] * res;
            }
        }
    }
}

int check(int n) { // 解系检测
    int status = 1;
    for (int i = n - 1; i >= 0; i--) {
        if (fabs(a[i][i]) < eps && fabs(a[i][n]) > eps) return -1; // 无解
        if (fabs(a[i][i]) < eps && fabs(a[i][n]) < eps) status = 0; // 无穷解
        ans[i] = a[i][n] / a[i][i];
        if (fabs(ans[i]) < eps) ans[i] = 0;
    }
    return status; // 唯一解或无穷解
}