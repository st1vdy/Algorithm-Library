bool is_prime[SIZE];
int prime[SIZE], num[SIZE]; // num[i] 表示 i 的质因子数
int getPrimeFactors(int n) { // 线性筛质因子数
    is_prime[1] = true;
    int p = 0;
    for (int i = 2; i <= n; i++) {
        if (!is_prime[i]) prime[p++] = i, num[i] = 1;
        for (int j = 0; j < p && prime[j] * i <= n; j++) {
            is_prime[prime[j] * i] = true;
            if (!(i % prime[j])) {
                num[i * prime[j]] = num[i];
                break;
            }
            num[i * prime[j]] = num[i] + 1;
        }
    }
    return p;
}