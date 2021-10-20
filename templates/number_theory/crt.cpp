// 求解形如 x = ci (mod mi) 的线性方程组 (mi, mj)必须两两互质
long long CRT(vector<long long>& c, vector<long long>& m) {
    long long M = m[0], ans = 0;
    for (int i = 1; i < (int)m.size(); ++i) M *= m[i];
    for (int i = 0; i < (int)m.size(); ++i) { // Mi * ti * ci
        long long mi = M / m[i];
        long long ti = inv(mi, m[i]); // mi 模 m[i] 的逆元
        ans = (ans + mi * ti % M * c[i] % M) % M;
    }
    ans = (ans + M) % M; // 返回模 M 意义下的唯一解
    return ans;
}