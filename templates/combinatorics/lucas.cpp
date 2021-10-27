long long p; // C(n, m) % p
long long lucas(long long n, long long m) { // O(p + log(n))
    if (m == 0) return 1;
    return (binom(n % p, m % p) * lucas(n / p, m / p)) % p; // binom(x, y) 在线性预处理组合数之后是 O(1) 的
}