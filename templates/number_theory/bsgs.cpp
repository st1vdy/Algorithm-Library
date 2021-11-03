ll bsgs(ll a, ll b, ll m) { // a ^ x = b mod m
    ll n = (ll)sqrt((double)m) + 1, base = 1, val = 1;
    map<int, int> mp; // 可以换 unordered_map
    b %= m;
    for (int i = 0; i < n; ++i) {
        mp[b * base % m] = i;
        base = (base * a) % m;
    }
    a = base;
    if (!a) return b ? -1 : 1;
    for (int i = 0; i <= n; ++i) {
        int j = (!mp.count(val) ? -1 : mp[val]);
        if (j >= 0 && i * n >= j) return i * n - j;
        val = (val * a) % m;
    }
    return -1; // 无解
}