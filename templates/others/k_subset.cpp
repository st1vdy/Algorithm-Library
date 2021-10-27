int n, k;
int s = (1 << k) - 1;
while (s < (1 << n)) { // O(binom(n, k))
    // 每次取出s就是一个大小为k的子集
    int x = s & -s, y = s + x;
    s = (((s & ~y) / x) >> 1) | y;
}