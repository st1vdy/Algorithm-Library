/*
 * 求解 ax+by+c=0 模板CF710D
 * 返回 [xl, xr] x [yl, yr] 的解数
 * 若至少有一组解 则x是一个合法解
 * 可以根据x推出y 但要注意b=0/a=0等特殊情况
 * 特别注意！！！设置边界时的取整问题！！！
 * x/y 向上取整时(x+y-1)/y 向下取整时floor(1.0*x/y)
 * */
ll a, b, c, xl, xr, yl, yr;
ll x, y, d;
ll exgcd(ll a, ll b, ll& x, ll& y) {
    if (!b) return x = 1, y = 0, a;
    ll d = exgcd(b, a % b, x, y), t = x;
    return x = y, y = t - a / b * y, d;
}
ll solve(ll a, ll b, ll c, ll xl, ll xr, ll yl, ll yr) {
    if (xl > xr) return 0;
    if (yl > yr) return 0;
    if (!a && !b) {
        if (c) return 0;
        return (xr - xl + 1) * (yr - yl + 1);
    }
    if (!b) {
        swap(a, b);
        swap(xl, yl);
        swap(xr, yr);
    }
    if (!a) {
        if (c % b) return 0;
        ll y = -c / b;
        if (y < yl || y > yr) return 0;
        return xr - xl + 1;
    }
    d = exgcd((a % abs(b) + abs(b)) % abs(b), abs(b), x, y);
    if (c % d) return 0;
    x = (x % abs(b) + abs(b)) % abs(b) * ((((-c) % abs(b)) + abs(b)) % abs(b) / d) % abs(b / d);
    d = abs(b / d);
    ll kl = (xl - x) / d - 3, kr = (xr - x) / d + 3;
    while (x + kl * d < xl) kl++;
    while (x + kr * d > xr) kr--;
    ll A = (-yl * b - a * x - c) / (a * d), B = (-yr * b - a * x - c) / (a * d);
    if (A > B) swap(A, B);
    kl = max(kl, A - 3);
    kr = min(kr, B + 3);
    while (kl <= kr) {
        ll y = (-c - a * x - a * d * kl) / b;
        if (yl <= y && y <= yr) break;
        kl++;
    }
    while (kl <= kr) {
        ll y = (-c - a * x - a * d * kr) / b;
        if (yl <= y && y <= yr) break;
        kr--;
    }
    if (kl > kr) return 0;
    return kr - kl + 1;
}