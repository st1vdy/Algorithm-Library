bool is_prime[SIZE];
int prime[SIZE], mu[SIZE]; // mu[i] 表示 i 的莫比乌斯函数值
int getMu(int n) { // 线性筛莫比乌斯函数
    mu[1] = 1; is_prime[1] = true;
    int p = 0;
    for (int i = 2; i <= n; i++) {
        if (!is_prime[i]) prime[p++] = i, mu[i] = -1;
        for (int j = 0; j < p && prime[j] * i <= n; j++) {
            is_prime[prime[j] * i] = true;
            if (!(i % prime[j])) {
                mu[i * prime[j]] = 0;
                break;
            }
            mu[i * prime[j]] = -mu[i];
        }
    }
    return p;
}