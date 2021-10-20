namespace Pollard_Rho {
    typedef long long ll;
    vector<ll> ans; // 存储质因子的数组
    inline ll gcd(ll a, ll b) { ll c; while (b) c = a % b, a = b, b = c; return a; }
    inline ll mulmod(ll x, ll y, const ll z) {
        return (x * y - (ll)(((long double)x * y + 0.5) / (long double)z) * z + z) % z;
    }
    inline ll powmod(ll a, ll b, const ll mo) {
        ll s = 1;
        for (; b; b >>= 1, a = mulmod(a, a, mo)) if (b & 1) s = mulmod(s, a, mo);
        return s;
    }
    bool isPrime(ll p) { // Miller-Rabin O(klog^3(n)) k为素性测试轮数
        const int lena = 10, a[lena] = { 2,3,5,7,11,13,17,19,23,29 };
        if (p == 2) return true;
        if (p == 1 || !(p & 1) || (p == 46856248255981ll)) return false;
        ll D = p - 1;
        while (!(D & 1)) D >>= 1;
        for (int i = 0; i < lena && a[i] < p; i++) {
            ll d = D, t = powmod(a[i], d, p);
            if (t == 1) continue;
            for (; d != p - 1 && t != p - 1; d <<= 1) t = mulmod(t, t, p);
            if (d == p - 1) return false;
        }
        return true;
    }
    void reportFactor(ll n) { // 得到一个素因子
        ans.emplace_back(n); // 存储素因子
    }
    ll ran() { return rand(); } // 随机数
    void getFactor(ll n) { // Pollard-Rho O(n ^ 1/4)
        if (n == 1) return;
        if (isPrime(n)) { reportFactor(n); return; }
        while (true) {
            ll c = ran() % n, i = 1, x = ran() % n, y = x, k = 2;
            do {
                ll d = gcd(n + y - x, n);
                if (d != 1 && d != n) { getFactor(d); getFactor(n / d); return; }
                if (++i == k) y = x, k <<= 1;
                x = (mulmod(x, x, n) + c) % n;
            } while (y != x);
        }
    }
}
using namespace Pollard_Rho;