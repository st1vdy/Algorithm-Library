template<typename T> struct xorTrie {
    int HIGHBIT, cnt;
    vector<vector<int>> nxt;
    vector<bool> vis;
    xorTrie(int n_ = 0, int highbit_ = 30) : HIGHBIT(highbit_), cnt(0) {
        int size_ = upperBoundEstimate(n_);
        nxt.resize(size_, vector<int>(2, 0));
        vis.resize(size_, false);
    }
    int upperBoundEstimate(int n) { // 求内存上界
        int hbit = log2(n);
        return n * (HIGHBIT - hbit + 1) + (1 << (hbit + 1)) - 1;
    }
    void insert(T x) { // 插入
        int p = 0;
        for (int i = HIGHBIT; ~i; i--) {
            int s = ((x >> i) & 1);
            if (!nxt[p][s]) nxt[p][s] = ++cnt;
            p = nxt[p][s];
        }
        vis[p] = true;
    }
    bool find(T x) { // 查询
        int p = 0;
        for (int i = HIGHBIT; ~i; i--) {
            int s = ((x >> i) & 1);
            if (!nxt[p][s]) return false;
            p = nxt[p][s];
        }
        return vis[p];
    }
};