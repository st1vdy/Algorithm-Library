template<typename T> struct fenwickTree {
    int n, hbit;
    vector<T> tree;
    fenwickTree(int n_ = 0) : n(n_), tree(n_ + 1), hbit(log2(n_) + 1) {}
    int lowbit(int x) { return x & (-x); }
    int size() { return n; }
    void add(int pos, int x) { // pos位置加上x
        for (; pos <= n; pos += lowbit(pos)) {
            tree[pos] += x;
        }
    }
    T query(int pos) { // 查询pos位置的前缀和 即a[1] + a[2] + ... + a[pos]
        T res = 0;
        for (; pos > 0; pos -= lowbit(pos)) {
            res += tree[pos];
        }
        return res;
    }
    T sum(int l, int r) { // [l, r]区间查询
        return query(r) - query(l - 1);
    }
    int kth(int k) { // 第k大元素
        int ans = 0, cnt = 0;
        for (int i = hbit; i >= 0; i--) {
            ans += (1 << i);
            if (ans > n || cnt + tree[ans] >= k) ans -= (1 << i);
            else cnt += tree[ans];
        }
        return ++ans;
    }
};