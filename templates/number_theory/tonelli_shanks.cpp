using ll = long long;
inline ll mulmod(ll x, ll y, const ll z) {
    return (x * y - (ll)(((long double)x * y + 0.5) / (long double)z) * z + z) % z;
}
inline ll powmod(ll a, ll b, const ll mo) {
    ll s = 1;
    for (; b; b >>= 1, a = mulmod(a, a, mo)) if (b & 1) s = mulmod(s, a, mo);
    return s;
}
ll tonelliShanks(ll n, ll p) { // O(log p)
    if (n == 0) return 0;
    if (p == 2) return (n & 1) ? 1 : -1;
    if (powmod(n, p >> 1, p) != 1) return -1;
    if (p & 2) return powmod(n, p + 1 >> 2, p);
    int s = __builtin_ctzll(p ^ 1);
    ll q = p >> s, z = 2;
    for (; powmod(z, p >> 1, p) == 1; ++z);
    ll c = powmod(z, q, p);
    ll r = powmod(n, q + 1 >> 1, p);
    ll t = powmod(n, q, p), tmp;
    for (int m = s, i; t != 1;) {
        for (i = 0, tmp = t; tmp != 1; i++) tmp = tmp * tmp % p;
        for (; i < --m;) c = c * c % p;
        r = r * c % p;
        c = c * c % p;
        t = t * c % p;
    }
    return r;
}