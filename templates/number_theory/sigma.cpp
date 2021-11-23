long long getSigma(long long n, long long mod = 1e18) {
    if (n == 1) return 1;
    long long ans = 1;
    for (long long i = 2; i * i <= n; ++i) {
        long long cnt = 1;
        while (!(n % i)) {
            cnt = (cnt * i + 1) % mod;
            n /= i;
        }
        ans *= cnt;
    }
    return (n == 1 ? ans : ans * (n + 1) % mod);
}