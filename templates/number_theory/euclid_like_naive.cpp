ll f(ll a, ll b, ll c, ll n) { // O(log n)
    if (a == 0)
        return (n + 1) * (b / c);
    if (a >= c || b >= c)
        return (f(a % c, b % c, c, n) + (a / c) * n * (n + 1) / 2 + (b / c) * (n + 1));
    ll m = (a * n + b) / c;
    return (n * m - f(c, c - b - 1, a, m - 1));
}