bool is_prime[SIZE];
int prime[SIZE], d[SIZE], num[SIZE]; // d[i] 表示 i 的因子数  num[i] 表示 i 的最小质因子出现次数
int getFactors(int n) { // 线性筛因子数
    d[1] = 1; is_prime[1] = true;
    int p = 0;
    for (int i = 2; i <= n; i++) {
        if (!is_prime[i]) prime[p++] = i, d[i] = 2, num[i] = 1;
        for (int j = 0; j < p && prime[j] * i <= n; j++) {
            is_prime[prime[j] * i] = true;
            if (!(i % prime[j])) {
                num[i * prime[j]] = num[i] + 1;
                d[i * prime[j]] = d[i] / num[i * prime[j]] * (num[i * prime[j]] + 1);
                break;
            }
            num[i * prime[j]] = 1;
            d[i * prime[j]] = d[i] + d[i];
        }
    }
    return p;
}