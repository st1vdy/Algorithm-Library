struct dsu {
private:
    // number of nodes
    int n;
    // root node: -1 * component size
    // otherwise: parent
    std::vector<int> pa;
public:
    dsu(int n_ = 0) : n(n_), pa(n_, -1) {}
    // find node x's parent
    int find(int x) {
        return pa[x] < 0 ? x : pa[x] = find(pa[x]);
    }
    // merge node x and node y 
    // if x and y had already in the same component, return false, otherwise return true
    // Implement (union by size) + (path compression)
    bool unite(int x, int y) {
        int xr = find(x), yr = find(y);
        if (xr != yr) {
            if (-pa[xr] < -pa[yr]) std::swap(xr, yr);
            pa[xr] += pa[yr];
            pa[yr] = xr; // y -> x
            return true;
        }
        return false;
    }
    // size of the connected component that contains the vertex x
    int size(int x) {
        return -pa[find(x)];
    }
};