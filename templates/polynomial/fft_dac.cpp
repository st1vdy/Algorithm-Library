/*
 * given n and g[1], g[2], ..., g[n], f[0] = 1
 * f[i] = sum_{j=1..i}f[i - j]g[j]
 * get f[1], f[2], ..., f[n] in O(nlog^2n)
 */
int n;
vector<int> f, g;
void dfs(int l, int r, int sgn) {
    int mid = (l + r) >> 1;
    if (l != r) {
        dfs(l, mid, 1);
        dfs(mid + 1, r, 0);
    }
    if (sgn) {
        int seg = r - l;
        vector<int> slice(f.begin() + l - 1, f.begin() + r);
        vector<int> cg(g.begin(), g.begin() + min(2 + 2 * seg, n));
        auto h = FFT::mod_multiply(slice, cg);
        for (int i = r; i <= min(r + seg, n - 1); i++) {
            f[i] = add(f[i], h[seg + i - r + 1]);
        }
    }
}