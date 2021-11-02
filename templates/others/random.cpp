int rd(int l, int r) {
    mt19937_64 gen(chrono::steady_clock::now().time_since_epoch().count());
    int p = uniform_int_distribution<int>(l, r)(gen);
    return p;
}