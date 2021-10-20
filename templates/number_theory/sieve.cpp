vector<bool> isPrime; // true 表示非素数  false 表示是素数
vector<int> prime; // 保存素数
int sieve(int n) {
    isPrime.resize(n + 1, false);
    isPrime[0] = isPrime[1] = true;
    for (int i = 2; i <= n; i++) {
        if (!isPrime[i]) prime.emplace_back(i);
        for (int j = 0; j < (int)prime.size() && prime[j] * i <= n; j++) {
            isPrime[prime[j] * i] = true;
            if (!(i % prime[j])) break;
        }
    }
    return (int)prime.size();
}