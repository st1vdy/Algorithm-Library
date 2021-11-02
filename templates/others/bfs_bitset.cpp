#include <bits/stdc++.h>
using namespace std;
const int SIZE = 1001;
bitset<SIZE> g[SIZE], vis, now;
int dis[SIZE][SIZE];
 
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);
    memset(dis, -1, sizeof(dis));
    int n;
    cin >> n;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            int x;
            cin >> x;
            if (x) g[i].set(j);
        }
    }
    for (int i = 1; i <= n; i++) {
        vis.reset(); // 清空已遍历的数组
        vis.set(i);
        queue<int> q;
        q.push(i);
        dis[i][i] = 0;
        while (!q.empty()) {
            auto top = q.front();
            q.pop();
            now = g[top] ^ (g[top] & vis); // 去掉已经遍历到的节点
             // 本方法的关键: O(n/w) 遍历bitset
            for (int to = now._Find_first(); to != now.size(); to = now._Find_next(to)) {
                dis[i][to] = dis[i][top] + 1;
                q.push(to);
            }
            vis |= now; // 更新已遍历的节点
        }
    }
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            cout << dis[i][j] << " \n"[j == n];
        }
    }
    return 0;
}