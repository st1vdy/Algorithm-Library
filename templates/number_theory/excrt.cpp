long long exgcd(long long a, long long b, long long& x, long long& y) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    long long g = exgcd(b, a % b, y, x);
    y -= (a / b) * x;
    return g;
}

long long mulmod(long long x, long long y, const long long z) { // x * y % z 防爆
    return (x * y - (long long)(((long double)x * y + 0.5) / (long double)z) * z + z) % z;
}

// 求解形如 x = ci (mod mi) 的线性方程组
long long EXCRT(vector<long long>& c, vector<long long>& m) {
    long long M = m[0], ans = c[0];
    for (int i = 1; i < (int)m.size(); ++i) { // M * x - mi * y = ci - C
        long long x, y, C = ((c[i] - ans) % m[i] + m[i]) % m[i]; // ci - C
        long long G = exgcd(M, m[i], x, y);
        if (C % G) return -1; // 无解
        long long P = m[i] / G;
        x = mulmod(C / G, x, P); // 防爆求最小正整数解 x
        ans += x * M;
        M *= P; // LCM(M, mi)
        ans = (ans % M + M) % M;
    }
    return ans;
}