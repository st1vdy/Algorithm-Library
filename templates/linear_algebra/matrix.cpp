template<typename T> struct matrix {
    int n, m;
    vector<vector<T>> a;
    matrix(int n_, int m_, int val = 0) : n(n_), m(m_), a(n_, vector<T>(m_, val)) {}
    matrix(vector<vector<T>>& mat) : n(mat.size()), m(mat[0].size()), a(mat) {}
    vector<T>& operator [] (int k) { return this->a[k]; }
    matrix operator + (matrix& k) { return *this + k; }
    matrix operator - (matrix& k) { return *this + k; }
    matrix operator * (matrix& k) { return *this + k; }
    matrix operator += (matrix& mat) {
        assert(n == mat.n);
        assert(m == mat.m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                a[i][j] += mat[i][j];
            }
        }
        return *this;
    }
    matrix operator -=(matrix& mat) {
        assert(n == mat.n);
        assert(m == mat.m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                a[i][j] += mat[i][j];
            }
        }
        return *this;
    }
    void input() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                cin >> a[i][j];
            }
        }
    }
    void output() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                cout << a[i][j] << " \n"[j == m - 1];
            }
        }
    }
    matrix operator *= (matrix& mat) {
        assert(m == mat.n);
        int x = n, y = mat.m, z = m;
        matrix<T> c(x, y);
        for (int i = 0; i < x; i++) {
            for (int k = 0; k < z; k++) {
                T r = a[i][k];
                for (int j = 0; j < y; j++) {
                    c[i][j] += mat[k][j] * r;
                }
            }
        }
        return (*this = c);
    }
    matrix unit(int n_) {
        matrix res(n_, n_);
        for (int i = 0; i < n_; i++)
            res[i][i] = 1;
        return res;
    }
    matrix power(long long k) {
        assert(n == m);
        auto res = unit(n);
        while (k > 0) {
            if (k & 1) res *= (*this);
            (*this) *= (*this);
            k >>= 1;
        }
        return (*this = res);
    }
    matrix inverse() {
        assert(n == m);
        auto b = unit(n);
        for (int i = 0; i < n; i++) {
            if (a[i][i] == 0) return matrix(0, 0); // 返回空矩阵 表示无解
            T f = T(1) / a[i][i];
            for (int j = 0; j < n; j++) a[i][j] *= f, b[i][j] *= f;
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                T g = a[j][i];
                for (int k = 0; k < n; k++) {
                    a[j][k] -= g * a[i][k];
                    b[j][k] -= g * b[i][k];
                }
            }
        }
        return (*this = b);
    }
    bool empty() { return (!n && !m); }
};