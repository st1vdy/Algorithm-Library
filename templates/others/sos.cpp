/* 高维前缀和/子集前缀和 */
for (int i = 0; i < n; i++) { // O(n2^n)
    for (int j = 0; j < (1 << n); j++) {
        if (j & (1 << i)) {
            pre[j] += pre[j ^ (1 << i)];
        }
    }
}