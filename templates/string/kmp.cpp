namespace KMP {
    vector<int> getPrefixTable(string s) { // 求前缀表
        int n = s.length();
        vector<int> nxt(n, 0);
        for (int i = 1; i < n; i++) {
            int j = nxt[i - 1];
            while (j > 0 && s[i] != s[j]) {
                j = nxt[j - 1];
            }
            if (s[i] == s[j]) j++;
            nxt[i] = j;
        }
        return nxt;
    }

    vector<int> kmp(string s, string t) { // 返回所有匹配位置的集合
        int n = s.length(), m = t.length();
        vector<int> res;
        vector<int> nxt = getPrefixTable(t);
        for (int i = 0, j = 0; i < n; i++) {
            while (j > 0 && j < m && s[i] != t[j]) {
                j = nxt[j - 1];
            }
            if (s[i] == t[j]) j++;
            if (j == m) {
                res.push_back(i + 1 - m);
                j = nxt[m - 1];
            }
        }
        return res;
    }
}