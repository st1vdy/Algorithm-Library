long long exgcd(long long a, long long b, long long& x, long long& y) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    long long g = exgcd(b, a % b, y, x);
    y -= (a / b) * x;
    return g;
}

ll x, y; // 最小非负整数解
bool solve(ll a, ll b, ll c) { // ax+by=c
    ll g = gcd(a, b);
    if (c % g) return false;
    a /= g, b /= g, c /= g;
    bool flag = false;
    if (b < 0) b = -b, flag = true;
    exgcd(a, b, x, y);
    x = (x * c % b + b) % b;
    if (flag) b = -b;
    y = (c - a * x) / b;
    if (!c) x = y = 0; // ax+by=0
    return true;
}