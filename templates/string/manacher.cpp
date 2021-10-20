namespace Manacher {
    static constexpr int SIZE = 1e5 + 5; // 预设为原串长度
    int len = 1; // manacher 预处理后字符串的长度
    char stk[SIZE << 1]; // manacher预处理字符串 需要2倍空间+1
    void init(string s) { // 初始化stk
        stk[0] = '*'; len = 1;
        for (int i = 0; i < s.length(); ++i) {
            stk[len++] = s[i];
            stk[len++] = '*';
        }
    }
    int manacher() { // 返回最长回文子串长度
        vector<int> rad(len << 1); // 存储每个点作为对称中心可拓展的最大半径
        int md = 0; // 最远回文串对称中心下标
        for (int i = 1; i < len; ++i) {
            int& r = rad[i] = 0;
            if (i <= md + rad[md]) {
                r = min(rad[2 * md - i], md + rad[md] - i);
            }
            while (i - r - 1 >= 0 && i + r + 1 < len &&
                stk[i - r - 1] == stk[i + r + 1]) ++r;
            if (i + r >= md + rad[md]) md = i;
        }
        int res = 0;
        for (int i = 0; i < len; ++i) {
            if (rad[i] > res) {
                res = rad[i];
            }
        }
        return res;
    }
}