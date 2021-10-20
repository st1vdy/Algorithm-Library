bool is_prime[SIZE];
int prime[SIZE], f[SIZE], g[SIZE]; // f[i] 表示 i 的约数和
int getSigma(int n) {
    g[1] = f[1] = 1; is_prime[1] = true;
    int p = 0;
    for (int i = 2; i <= n; i++) {
        if (!is_prime[i]) prime[p++] = i, f[i] = g[i] = i + 1;
        for (int j = 0; j < p && prime[j] * i <= n; j++) {
            is_prime[prime[j] * i] = true;
            if (!(i % prime[j])) {
                g[i * prime[j]] = g[i] * prime[j] + 1;
                f[i * prime[j]] = f[i] / g[i] * g[i * prime[j]];
                break;
            }
            f[i * prime[j]] = f[i] * f[prime[j]];
            g[i * prime[j]] = 1 + prime[j];
        }
    }
    return p;
}