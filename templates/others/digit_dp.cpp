#include <bits/stdc++.h>
using namespace std;
using ll = long long;
ll a[20], dp[20][20];
ll dfs(int len, int las, int maxi, int lead) {
    if (len == 0) return 1;
    if (!maxi && !lead && dp[len][las] != -1) return dp[len][las];
    ll sum = 0;
    int ma = 9;
    if (maxi) ma = a[len];
    if (lead) {
        for (int i = 0; i <= ma; i++) {
            if (i == 0) sum += dfs(len - 1, i, 0, lead);
            else if (i == ma && maxi) sum += dfs(len - 1, i, maxi, 0);
            else sum += dfs(len - 1, i, 0, 0);
        }
    } else {
        for (int i = 0; i <= ma; i++) {
            if (abs(i - las) < 2) continue;
            if (i == ma && maxi) sum += dfs(len - 1, i, maxi, 0);
            else sum += dfs(len - 1, i, 0, 0);
        }
    }
    if (maxi == 0 && lead == 0) dp[len][las] = sum;
    return sum;
}
ll sol(ll x) {
    int cnt = 0;
    a[++cnt] = x % 10;
    x /= 10;
    while (x) {
        a[++cnt] = x % 10;
        x /= 10;
    }
    return dfs(cnt, 0, 1, 1);
}

int main() {
    ll l, r;
    cin >> l >> r;
    memset(dp, -1, sizeof(dp));
    cout << sol(r) - sol(l - 1) << "\n";
    return 0;
}
