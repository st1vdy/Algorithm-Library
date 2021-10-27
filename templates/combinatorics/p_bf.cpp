p[0][0] = 1;
for (int i = 1; i <= n; i++) {
    p[0][i] = 1;
    for (int j = 1; j <= m; j++) {
        if (i >= j) {
            p[i][j] = add(p[i][j - 1], p[i - j][j]);
        } else {
            p[i][j] = p[i][j - 1];
        }
    }
}