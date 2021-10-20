constexpr int SIZE = 1001;
bitset<SIZE> a[SIZE];
int ans[SIZE];
void gauss(int n) { // bitset版高斯消元 用于求解异或线性方程组
    bitset<SIZE> vis;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (vis[j]) continue;
            if (a[j][i]) {
                vis.set(i);
                swap(a[i], a[j]);
                break;
            }
        }
        if (!a[i][i]) continue;
        for (int j = 0; j <= n; j++) {
            if (i != j && (a[j][i] & a[i][i])) {
                a[j] ^= a[i];
            }
        }
    }
}