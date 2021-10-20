bool is_prime[SIZE];
int prime[SIZE], phi[SIZE]; // phi[i] 表示 i 的欧拉函数值
int Phi(int n) { // 线性筛素数的同时线性求欧拉函数
    phi[1] = 1; is_prime[1] = true;
    int p = 0;
    for (int i = 2; i <= n; i++) {
        if (!is_prime[i]) prime[p++] = i, phi[i] = i - 1;
        for (int j = 0; j < p && prime[j] * i <= n; j++) {
            is_prime[prime[j] * i] = true;
            if (!(i % prime[j])) {
                phi[i * prime[j]] = phi[i] * prime[j];
                break;
            }
            phi[i * prime[j]] = phi[i] * (prime[j] - 1);
        }
    }
    return p;
}