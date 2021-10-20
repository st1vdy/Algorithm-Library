namespace KM {
    typedef long long ll;
    const int maxn = 510;
    const int inf = 1e9;
    int vx[maxn], vy[maxn], lx[maxn], ly[maxn], slack[maxn];
    int w[maxn][maxn]; // 以上为权值类型
    int pre[maxn], left[maxn], right[maxn], NL, NR, N;
    void match(int& u) {
        for (; u; std::swap(u, right[pre[u]]))
            left[u] = pre[u];
    }
    void bfs(int u) {
        static int q[maxn], front, rear;
        front = 0; vx[q[rear = 1] = u] = true;
        while (true) {
            while (front < rear) {
                int u = q[++front];
                for (int v = 1; v <= N; ++v) {
                    int tmp;
                    if (vy[v] || (tmp = lx[u] + ly[v] - w[u][v]) > slack[v])
                        continue;
                    pre[v] = u;
                    if (!tmp) {
                        if (!left[v]) return match(v);
                        vy[v] = vx[q[++rear] = left[v]] = true;
                    } else slack[v] = tmp;
                }
            }
            int a = inf;
            for (int i = 1; i <= N; ++i)
                if (!vy[i] && a > slack[i]) a = slack[u = i];
            for (int i = 1; i <= N; ++i) {
                if (vx[i]) lx[i] -= a;
                if (vy[i]) ly[i] += a;
                else slack[i] -= a;
            }
            if (!left[u]) return match(u);
            vy[u] = vx[q[++rear] = left[u]] = true;

        }

    }
    void exec() {
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                slack[j] = inf;
                vx[j] = vy[j] = false;
            }
            bfs(i);
        }
    }
    ll work(int nl, int nr) { // NL , NR 为左右点数, 返回最大权匹配的权值和
        NL = nl; NR = nr;
        N = std::max(NL, NR);
        for (int u = 1; u <= N; ++u)
            for (int v = 1; v <= N; ++v)
                lx[u] = std::max(lx[u], w[u][v]);
        exec();
        ll ans = 0;
        for (int i = 1; i <= N; ++i)
            ans += lx[i] + ly[i];
        return ans;
    }
    void output() { // 输出左边点与右边哪个点匹配, 没有匹配输出0
        for (int i = 1; i <= NL; ++i)
            printf("%d ", (w[i][right[i]] ? right[i] : 0));
        printf("\n");
    }
}