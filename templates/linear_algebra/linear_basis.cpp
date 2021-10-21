struct linearBasis {
    /* 线性基性质：
     * 1.若a[i]!=0（即主元i存在）
     *   则线性基中只有a[i]的第i位是1（只存在一个主元）
     *   且此时a[i]的最高位就是第i位
     * 2.将数组a插入线性基 假设有|B|个元素成功插入
     *   则数组a中每个不同的子集异或和都出现 2^(n-|B|) 次
     * */
    static const int MAXL = 60;
    long long a[MAXL + 1];
    int id[MAXL + 1];
    int zero;
    /* 0的标志位 =1则表示0可以被线性基表示出来
     * 求第k大元素时 需要注意题意中线性基为空时是否可以表示0
     * 默认不可以表示 有必要时进行修改即可
     * */
    linearBasis() {
        zero = 0;
        fill(a, a + MAXL + 1, 0);
    }
    long long& operator[] (int k) { return a[k]; }
    bool insert(long long x) {
        for (int j = MAXL; ~j; j--) {
            if (!(x & (1LL << j))) { // 如果 x 的第 j 位为 0，则跳过
                continue;
            }
            if (a[j]) { // 如果 a[j] != 0，则用 a[j] 消去 x 的第 j 位上的 1
                x ^= a[j];
            } else { // 找到插入位置
                for (int k = 0; k < j; k++) {
                    if (x & (1LL << k)) { // 如果x存在某个低位线性基的主元k则消去
                        x ^= a[k];
                    }
                }
                for (int k = j + 1; k <= MAXL; k++) {
                    if (a[k] & (1LL << j)) { // 如果某个高位线性基存在主元j则消去
                        a[k] ^= x;
                    }
                }
                a[j] = x;
                return true;
            }
        }
        zero = 1;
        return false;
    }
    long long query_max() { // 最大值
        long long res = 0;
        for (int i = MAXL; ~i; i--) {
            res ^= a[i];
        }
        return res;
    }
    long long query_max(long long x) { // 线性基异或x的最大值
        for (int i = MAXL; ~i; i--) {
            if ((x ^ a[i]) > x) {
                x ^= a[i];
            }
        }
        return x;
    }
    long long query_min() { // 最小值
        for (int i = 0; i < MAXL; i++) {
            if (a[i]) {
                return a[i];
            }
        }
        return -1; // 线性基为空
    }
    long long query_min(long long x) { // 线性基异或x的最小值
        for (int i = MAXL; ~i; i--) {
            if ((x ^ a[i]) < x) {
                x ^= a[i];
            }
        }
        return x;
    }
    int count(long long x) { // 元素 x 能否被线性基内元素表示
        int res = 0;
        vector<long long> b(MAXL + 1);
        for (int i = 0; i <= MAXL; i++) {
            b[i] = a[i];
        }
        res = this->insert(x);
        for (int i = 0; i <= MAXL; i++) {
            a[i] = b[i];
        }
        return !res; // 成功插入则无法表示 失败则可以表示
    }
    int size() { // 线性基有效元素数量
        int res = 0;
        for (int i = 0; i <= MAXL; i++) {
            if (a[i]) {
                res++;
            }
        }
        return res;
    }
    long long kth_element(long long k) { // 第k大元素
        vector<long long> b;
        for (int i = 0; i <= MAXL; i++) {
            if (a[i]) {
                b.push_back(a[i]);
            }
        }
        if (zero) {
            if (--k == 0) {
                return 0;
            }
        }
        if (k >= (1LL << this->size())) { // k超过了线性基可以表示的最大数量
            return -1;
        }
        long long res = 0;
        for (int i = 0; i <= MAXL; i++) {
            if (k & (1LL << i)) {
                res ^= b[i];
            }
        }
        return res;
    }
    long long rank(long long x) { // 元素x在线性基内的排名（默认不考虑0）
        vector<long long> b;
        for (int i = 0; i <= MAXL; i++) {
            if (a[i]) {
                b.push_back(1LL << i);
            }
        }
        long long res = 0;
        for (int i = 0; i < (int)b.size(); i++) {
            if (x & b[i]) {
                res |= (1LL << i);
            }
        }
        return res;
    }
    void clear() {
        zero = 0;
        fill(a, a + MAXL + 1, 0);
    }
};