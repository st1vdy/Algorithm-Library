/* O(3^n)遍历子集 */
for (int i = 0; i < (1 << n); i++) { // i是当前需要遍历的全集
    for (int l = i;; l = i & (l - 1)) { // l是i的子集
        int r = i - l; // l+r=i
        if (!l) break;
    }
}