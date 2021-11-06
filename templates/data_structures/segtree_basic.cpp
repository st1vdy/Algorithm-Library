namespace SegmentTree {
    static constexpr int SIZE = 200001;
    int a[SIZE];
    int mx[SIZE << 2], lazy[SIZE << 2];
 
    void pushup(int rt) {
        mx[rt] = max(mx[rt << 1], mx[rt << 1 | 1]);
    }
 
    void pushdown(int rt) {
        if (lazy[rt]) {
            lazy[rt << 1] += lazy[rt];
            lazy[rt << 1 | 1] += lazy[rt];
            mx[rt << 1] += lazy[rt];
            mx[rt << 1 | 1] += lazy[rt];
            lazy[rt] = 0;
        }
    }
 
    void build(int rt, int l, int r) {
        if (l == r) {
            mx[rt] = a[l];
            lazy[rt] = 0;
            return;
        }
        int mid = (l + r) >> 1;
        pushdown(rt);
        build(rt << 1, l, mid), build(rt << 1 | 1, mid + 1, r);
        pushup(rt);
    }
 
    void rangeAdd(int rt, int l, int r, int ql, int qr, int val) {
        if (ql <= l && qr >= r) {
            mx[rt] += val;
            lazy[rt] += val;
            return;
        }
        int mid = (l + r) >> 1;
        pushdown(rt);
        if (ql <= mid) rangeAdd(rt << 1, l, mid, ql, qr, val);
        if (qr > mid) rangeAdd(rt << 1 | 1, mid + 1, r, ql, qr, val);
        pushup(rt);
    }
 
    int rangeMax(int rt, int l, int r, int ql, int qr) {
        if (ql <= l && qr >= r) return mx[rt];
        int mid = (l + r) >> 1;
        int res = 0;
        pushdown(rt);
        if (ql <= mid) res = max(res, rangeMax(rt << 1, l, mid, ql, qr));
        if (qr > mid) res = max(res, rangeMax(rt << 1 | 1, mid + 1, r, ql, qr));
        pushup(rt);
        return res;
    }
}