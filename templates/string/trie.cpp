struct trie {
    int cnt;
    vector<vector<int>> nxt;
    vector<bool> vis;
    /* 初始化的时候size需要设置为字符串总长之和 26是字符集大小 */
    trie(int size_ = 0) :cnt(0), vis(size_, false), nxt(size_, vector<int>(26, 0)) {}
    void insert(string s) {  // 插入字符串
        int p = 0;
        for (int i = 0; i < (int)s.length(); i++) {
            int c = s[i] - 'a';
            if (!nxt[p][c]) nxt[p][c] = ++cnt;
            p = nxt[p][c];
        }
        vis[p] = true;
    }
    bool find(string s) {  // 查找字符串
        int p = 0;
        for (int i = 0; i < (int)s.length(); i++) {
            int c = s[i] - 'a';
            if (!nxt[p][c]) return false;
            p = nxt[p][c];
        }
        return vis[p];
    }
};