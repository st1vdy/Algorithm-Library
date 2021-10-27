vector<int> prime;  // 求 n 的欧拉函数需要先把 <= ceil(sqrt(n)) 的素数筛出
long long phi(long long n) { // O(sqrt(N)/log(N))
    long long res = n;
    for (int i = 0; i < (int)prime.size(); i++) {
        long long p = prime[i];
        if (p * p > n) break;
        if (n % p == 0) {
            res = res / p * (p - 1);
            while (n % p == 0) n /= p;
        }
    }
    if (n > 1) res = res / n * (n - 1);
    return res;
}