namespace Geometry {
#define db long double
#define pi acos(-1.0)
    constexpr db eps = 1e-7;
    int sign(db k) {
        if (k > eps) return 1;
        else if (k < -eps) return -1;
        return 0;
    }
    int cmp(db k1, db k2) { // k1 < k2 : -1, k1 == k2 : 0, k1 > k2 : 1
        return sign(k1 - k2);
    }
    int inmid(db k1, db k2, db k3) { // k3 在 [k1, k2] 内
        return sign(k1 - k3) * sign(k2 - k3) <= 0;
    }

    struct point { // 点类
        db x, y;
        point() {}
        point(db x_, db y_) :x(x_), y(y_) {}
        point operator + (const point& k) const { return point(k.x + x, k.y + y); }
        point operator - (const point& k) const { return point(x - k.x, y - k.y); }
        point operator * (db k) const { return point(x * k, y * k); }
        point operator / (db k1) const { return point(x / k1, y / k1); }
        point turn(db k1) { return point(x * cos(k1) - y * sin(k1), x * sin(k1) + y * cos(k1)); } // 逆时针旋转
        point turn90() { return point(-y, x); } // 逆时针方向旋转 90 度
        db len() { return sqrt(x * x + y * y); } // 向量长度
        db len2() { return x * x + y * y; } // 向量长度的平方
        db getPolarAngle() { return atan2(y, x); } // 向量极角
        db dis(point k) { return ((*this) - k).len(); } // 到点k的距离
        point unit() { db d = len(); return point(x / d, y / d); } // 单位向量
        point getdel() { // 将向量的方向调整为指向第一/四象限 包括y轴正方向
            if (sign(x) == -1 || (sign(x) == 0 && sign(y) == -1))
                return (*this) * (-1);
            else return (*this);
        }
        bool operator < (const point& k) const { // 水平序排序 x坐标为第一关键字,y坐标第二关键字
            return x == k.x ? y < k.y : x < k.x;
        }
        bool operator == (const point& k) const { return cmp(x, k.x) == 0 && cmp(y, k.y) == 0; }
        bool getP() const { // 判断点是否在上半平面 含x负半轴 不含x正半轴及零点
            return sign(y) == 1 || (sign(y) == 0 && sign(x) == -1);
        }
        void input() { cin >> x >> y; }
    };
    db cross(point k1, point k2) { return k1.x * k2.y - k1.y * k2.x; } // 向量 k1,k2 的叉积
    db dot(point k1, point k2) { return k1.x * k2.x + k1.y * k2.y; }   // 向量 k1,k2 的点积
    db rad(point k1, point k2) { // 向量 k1,k2 之间的有向夹角
        return atan2(cross(k1, k2), dot(k1, k2));
    }
    int inmid(point k1, point k2, point k3) { // k1 k2 k3共线时 判断点 k3 是否在线段 k1k2 上
        return inmid(k1.x, k2.x, k3.x) && inmid(k1.y, k2.y, k3.y);
    }
    int compareAngle(point k1, point k2) { // 比较向量 k1,k2 的角度大小 角度按照atan2()函数定义
        // k1 < k2 返回 1, k1 >= k2 返回 0
        return k1.getP() < k2.getP() || (k1.getP() == k2.getP() && sign(cross(k1, k2)) > 0);
    }
    point proj(point k1, point k2, point q) { // q 到直线 k1,k2 的投影
        point k = k2 - k1; return k1 + k * (dot(q - k1, k) / k.len2());
    }
    point reflect(point k1, point k2, point q) { return proj(k1, k2, q) * 2 - q; } // q 关于直线 k1,k2 的对称点
    int counterclockwise(point k1, point k2, point k3) { // k1 k2 k3 逆时针1 顺时针-1 否则0
        return sign(cross(k2 - k1, k3 - k1));
    }
    int checkLL(point k1, point k2, point k3, point k4) { // 判断直线 k1k2 和直线k3k4 是否相交
        // 即判断直线 k1k2 和 k3k4 是否平行 平行返回0 不平行返回1
        return sign(cross(k2 - k1, k4 - k3)) != 0;
    }
    point getLL(point k1, point k2, point k3, point k4) { // 求 k1k2 k3k4 两直线交点
        db w1 = cross(k1 - k3, k4 - k3), w2 = cross(k4 - k3, k2 - k3);
        return (k1 * w2 + k2 * w1) / (w1 + w2);
    }
    int intersect(db l1, db r1, db l2, db r2) { // 判断 [l1,r1] 和 [l2, r2] 是否相交
        if (l1 > r1) swap(l1, r1);
        if (l2 > r2) swap(l2, r2);
        return cmp(r1, l2) != -1 && cmp(r2, l1) != -1;
    }
    int checkSS(point k1, point k2, point k3, point k4) { // 判断线段 k1k2 和线段 k3k4 是否相交
        return intersect(k1.x, k2.x, k3.x, k4.x) && intersect(k1.y, k2.y, k3.y, k4.y) &&
            sign(cross(k3 - k1, k4 - k1)) * sign(cross(k3 - k2, k4 - k2)) <= 0 &&
            sign(cross(k1 - k3, k2 - k3)) * sign(cross(k1 - k4, k2 - k4)) <= 0;
    }
    db disSP(point k1, point k2, point q) { // 点 q 到线段 k1k2 的最短距离
        point k3 = proj(k1, k2, q);
        if (inmid(k1, k2, k3)) return q.dis(k3);
        else return min(q.dis(k1), q.dis(k2));
    }
    db disLP(point k1, point k2, point q) { // 点 q 到直线 k1k2 的最短距离
        point k3 = proj(k1, k2, q);
        return q.dis(k3);
    }
    db disSS(point k1, point k2, point k3, point k4) { // 线段 k1k2 和线段 k3k4 的最短距离
        if (checkSS(k1, k2, k3, k4)) return 0;
        else return min(min(disSP(k1, k2, k3), disSP(k1, k2, k4)),
            min(disSP(k3, k4, k1), disSP(k3, k4, k2)));
    }
    bool onLine(point k1, point k2, point q) { // 判断点 q 是否在直线 k1k2 上
        return sign(cross(k1 - q, k2 - q)) == 0;
    }
    bool onSegment(point k1, point k2, point q) { // 判断点 q 是否在线段 k1k2 上
        if (!onLine(k1, k2, q)) return false; // 如果确定共线 要删除这个特判
        return inmid(k1, k2, q);
    }
    void polarAngleSort(vector<point>& p, point t) { // p为待排序点集 t为极角排序中心
        sort(p.begin(), p.end(), [&](const point& k1, const point& k2) {
            return compareAngle(k1 - t, k2 - t);
        });
    }

    struct line { // 直线 / 线段类
        point p[2];
        line() {}
        line(point k1, point k2) { p[0] = k1, p[1] = k2; }
        point& operator [] (int k) { return p[k]; }
        point dir() { return p[1] - p[0]; } // 向量 p[0] -> p[1]
        bool include(point k) { // 判断点是否在直线上
            return sign(cross(p[1] - p[0], k - p[0])) > 0;
        }
        bool includeS(point k) { // 判断点是否在线段上
            return onSegment(p[0], p[1], k);
        }
        line push(db len) { // 向外（左手边）平移 len 个单位
            point delta = (p[1] - p[0]).turn90().unit() * len;
            return line(p[0] - delta, p[1] - delta);
        }
    };
    bool parallel(line k1, line k2) { // 判断是否平行
        return sign(cross(k1.dir(), k2.dir())) == 0;
    }
    bool sameLine(line k1, line k2) { // 判断是否共线
        return parallel(k1, k2) && parallel(k1, line(k2.p[0], k1.p[0]));
    }
    bool sameDir(line k1, line k2) { // 判断向量 k1 k2 是否同向
        return parallel(k1, k2) && sign(dot(k1.dir(), k2.dir())) == 1;
    }
    bool operator < (line k1, line k2) {
        if (sameDir(k1, k2)) return k2.include(k1[0]);
        return compareAngle(k1.dir(), k2.dir());
    }
    bool checkLL(line k1, line k2) {
        return checkLL(k1[0], k1[1], k2[0], k2[1]);
    }
    point getLL(line k1, line k2) {  // 求 k1 k2 两直线交点 不要忘了判平行!
        return getLL(k1[0], k1[1], k2[0], k2[1]);
    }
    bool checkpos(line k1, line k2, line k3) {  // 判断是否三线共点
        return k3.include(getLL(k1, k2));
    }

    struct circle { // 圆类
        point o;
        double r;
        circle() {}
        circle(point o_, double r_) : o(o_), r(r_) {}
        int inside(point k) {  // 判断点 k 和圆的位置关系
            return cmp(r, o.dis(k)); // 圆外:-1, 圆上:0, 圆内:1
        }
    };
    int checkposCC(circle k1, circle k2) { // 返回两个圆的公切线数量
        if (cmp(k1.r, k2.r) == -1) swap(k1, k2);
        db dis = k1.o.dis(k2.o);
        int w1 = cmp(dis, k1.r + k2.r), w2 = cmp(dis, k1.r - k2.r);
        if (w1 > 0) return 4; // 外离
        else if (w1 == 0) return 3; // 外切
        else if (w2 > 0) return 2;  // 相交
        else if (w2 == 0) return 1; // 内切
        else return 0; // 内离(包含)
    }
    vector<point> getCL(circle k1, point k2, point k3) { // 求直线 k2k3 和圆 k1 的交点
        // 沿着 k2->k3 方向给出 相切给出两个
        point k = proj(k2, k3, k1.o);
        db d = k1.r * k1.r - (k - k1.o).len2();
        if (sign(d) == -1) return {};
        point del = (k3 - k2).unit() * sqrt(max((db)0.0, d));
        return { k - del,k + del };
    }
    vector<point> getCC(circle k1, circle k2) { // 求圆 k1 和圆 k2 的交点
        // 沿圆 k1 逆时针给出, 相切给出两个
        int pd = checkposCC(k1, k2); if (pd == 0 || pd == 4) return {};
        db a = (k2.o - k1.o).len2(), cosA = (k1.r * k1.r + a -
            k2.r * k2.r) / (2 * k1.r * sqrt(max(a, (db)0.0)));
        db b = k1.r * cosA, c = sqrt(max((db)0.0, k1.r * k1.r - b * b));
        point k = (k2.o - k1.o).unit(), m = k1.o + k * b, del = k.turn90() * c;
        return { m - del,m + del };
    }
    vector<point> tangentCP(circle k1, point k2) { // 点 k2 到圆 k1 的切点 沿圆 k1 逆时针给出
        db a = (k2 - k1.o).len(), b = k1.r * k1.r / a, c = sqrt(max((db)0.0, k1.r * k1.r - b * b));
        point k = (k2 - k1.o).unit(), m = k1.o + k * b, del = k.turn90() * c;
        return { m - del,m + del };
    }
    vector<line> tangentOutCC(circle k1, circle k2) {
        int pd = checkposCC(k1, k2);
        if (pd == 0) return {};
        if (pd == 1) {
            point k = getCC(k1, k2)[0];
            return { line(k,k) };
        }
        if (cmp(k1.r, k2.r) == 0) {
            point del = (k2.o - k1.o).unit().turn90().getdel();
            return { line(k1.o - del * k1.r,k2.o - del * k2.r),
                line(k1.o + del * k1.r,k2.o + del * k2.r) };
        } else {
            point p = (k2.o * k1.r - k1.o * k2.r) / (k1.r - k2.r);
            vector<point> A = tangentCP(k1, p), B = tangentCP(k2, p);
            vector<line> ans; for (int i = 0; i < A.size(); i++)
                ans.push_back(line(A[i], B[i]));
            return ans;
        }
    }
    vector<line> tangentInCC(circle k1, circle k2) {
        int pd = checkposCC(k1, k2);
        if (pd <= 2) return {};
        if (pd == 3) {
            point k = getCC(k1, k2)[0];
            return { line(k, k) };
        }
        point p = (k2.o * k1.r + k1.o * k2.r) / (k1.r + k2.r);
        vector<point> A = tangentCP(k1, p), B = tangentCP(k2, p);
        vector<line> ans;
        for (int i = 0; i < (int)A.size(); i++) ans.push_back(line(A[i], B[i]));
        return ans;
    }
    vector<line> tangentCC(circle k1, circle k2) { // 求两圆公切线
        int flag = 0;
        if (k1.r < k2.r) swap(k1, k2), flag = 1;
        vector<line> A = tangentOutCC(k1, k2), B = tangentInCC(k1, k2);
        for (line k : B) A.push_back(k);
        if (flag) for (line& k : A) swap(k[0], k[1]);
        return A;
    }
    db getAreaCT(circle k1, point k2, point k3) { // 圆 k1 与三角形 k2k3k1.o 的有向面积交
        point k = k1.o; k1.o = k1.o - k; k2 = k2 - k; k3 = k3 - k;
        int pd1 = k1.inside(k2), pd2 = k1.inside(k3);
        vector<point> A = getCL(k1, k2, k3);
        if (pd1 >= 0) {
            if (pd2 >= 0) return cross(k2, k3) / 2;
            return k1.r * k1.r * rad(A[1], k3) / 2 + cross(k2, A[1]) / 2;
        } else if (pd2 >= 0) {
            return k1.r * k1.r * rad(k2, A[0]) / 2 + cross(A[0], k3) / 2;
        } else {
            int pd = cmp(k1.r, disSP(k2, k3, k1.o));
            if (pd <= 0) return k1.r * k1.r * rad(k2, k3) / 2;
            return cross(A[0], A[1]) / 2 + k1.r * k1.r * (rad(k2, A[0]) + rad(A[1], k3)) / 2;
        }
    }
    db getAreaCC(circle k1, circle k2) { // 两圆面积交
        db d = k1.o.dis(k2.o);
        if (cmp(d, k1.r + k2.r) >= 0) return 0; // 两圆相离
        if (cmp(k1.r, k2.r) == -1) swap(k1, k2);
        if (cmp(k1.r - k2.r, d) >= 0) return pi * k2.r * k2.r; // 圆k1包含k2
        db g1 = acos((k1.r * k1.r + d * d - k2.r * k2.r) / (2 * k1.r * d));
        db g2 = acos((k2.r * k2.r + d * d - k1.r * k1.r) / (2 * k2.r * d));
        return g1 * k1.r * k1.r + g2 * k2.r * k2.r - k1.r * d * sin(g1);
    }
    circle getCircleOut(point k1, point k2, point k3) { // 三角形外接圆
        db a1 = k2.x - k1.x, b1 = k2.y - k1.y, c1 = (a1 * a1 + b1 * b1) / 2;
        db a2 = k3.x - k1.x, b2 = k3.y - k1.y, c2 = (a2 * a2 + b2 * b2) / 2;
        db d = a1 * b2 - a2 * b1;
        point o(k1.x + (c1 * b2 - c2 * b1) / d, k1.y + (a1 * c2 - a2 * c1) / d);
        return circle(o, k1.dis(o));
    }
    circle getCircleIn(point k1, point k2, point k3) {  // 三角形内切圆
        db a = k1.dis(k2), b = k2.dis(k3), c = k3.dis(k1);
        db len = a + b + c;
        db r = abs(cross(k1 - k2, k1 - k3)) / len;
        point o((k1.x * b + k2.x * c + k3.x * a) / len, (k1.y * b + k2.y * c + k3.y * a) / len);
        return circle(o, r);
    }
    circle minCircleCovering(vector<point> A) { // 最小圆覆盖 O(n)随机增量法
        // random_shuffle(A.begin(), A.end()); // <= C++14
        auto seed = chrono::steady_clock::now().time_since_epoch().count();
        default_random_engine e(seed);
        shuffle(A.begin(), A.end(), e); // >= C++11
        circle ans = circle(A[0], 0);
        for (int i = 1; i < A.size(); i++) {
            if (ans.inside(A[i]) == -1) {
                ans = circle(A[i], 0);
                for (int j = 0; j < i; j++) {
                    if (ans.inside(A[j]) == -1) {
                        ans.o = (A[i] + A[j]) / 2;
                        ans.r = ans.o.dis(A[i]);
                        for (int k = 0; k < j; k++) {
                            if (ans.inside(A[k]) == -1)
                                ans = getCircleOut(A[i], A[j], A[k]);
                        }
                    }
                }
            }
        }
        return ans;
    }

    typedef vector<point> polygon;
    db area(polygon p) { // 多边形有向面积
        if (p.size() < 3) return 0;
        db ans = 0;
        for (int i = 1; i < p.size() - 1; i++)
            ans += cross(p[i] - p[0], p[i + 1] - p[0]);
        return 0.5L * ans;
    }
    
    int checkConvexP(polygon p, point a) { // O(logn)判断点是否在凸包内 2内部 1边界 0外部
        // 必须保证凸多边形是一个水平序凸包且不能退化
        // 退化情况 比如凸包退化成线段 可使用 onSegment() 函数特判
        auto check = [&](int x) {
            int ccw1 = counterclockwise(p[0], a, p[x]),
                ccw2 = counterclockwise(p[0], a, p[x + 1]);
            if (ccw1 == -1 && ccw2 == -1) return 2;
            else if (ccw1 == 1 && ccw2 == 1) return 0;
            else if (ccw1 == -1 && ccw2 == 1) return 1;
            else return 1;
        };
        if (counterclockwise(p[0], a, p[1]) <= 0 && counterclockwise(p[0], a, p.back()) >= 0) {
            int l = 1, r = p.size() - 2, mid;
            while (l <= r) {
                mid = (l + r) >> 1;
                int chk = check(mid);
                if (chk == 1) l = mid + 1;
                else if (chk == -1) r = mid;
                else break;
            }
            int res = counterclockwise(p[mid], a, p[mid + 1]);
            if (res < 0) return 2;
            else if (res == 0) return 1;
            else return 0;
        } else {
            return 0;
        }
    }
    int checkPolyP(vector<point> p, point q) { // O(n)判断点是否在一般多边形内
        // 必须保证简单多边形的点按逆时针给出 返回 2 内部 1 边界 0 外部
        int pd = 0, n = p.size();
        for (int i = 0; i < n; i++) {
            point u = p[i], v = p[(i + 1) % n];
            if (onSegment(u, v, q)) return 1;
            if (cmp(u.y, v.y) > 0) swap(u, v);
            if (cmp(u.y, q.y) >= 0 || cmp(v.y, q.y) < 0) continue;
            if (sign(cross(u - v, q - v)) < 0) pd ^= 1;
        }
        return pd << 1;
    }
    db convexDiameter(polygon p) { // 0(n)旋转卡壳求凸包直径 / 平面最远点对的平方
        int n = p.size(); // 请保证多边形是凸包
        db ans = 0;
        for (int i = 0, j = n < 2 ? 0 : 1; i < j; i++) {
            for (;; j = (j + 1) % n) {
                ans = max(ans, (p[i] - p[j]).len2());
                if (sign(cross(p[i + 1] - p[i], p[(j + 1) % n] - p[j])) <= 0) break;
            }
        }
        return ans;
    }
    polygon convexHull(polygon A, int flag = 1) { // 凸包 flag=0 不严格 flag=1 严格
        int n = A.size(); polygon ans(n + n);
        sort(A.begin(), A.end()); int now = -1;
        for (int i = 0; i < A.size(); i++) {
            while (now > 0 && sign(cross(ans[now] - ans[now - 1], A[i] - ans[now - 1])) < flag)
                now--;
            ans[++now] = A[i];
        }
        int pre = now;
        for (int i = n - 2; i >= 0; i--) {
            while (now > pre && sign(cross(ans[now] - ans[now - 1], A[i] - ans[now - 1])) < flag)
                now--;
            ans[++now] = A[i];
        }
        ans.resize(now);
        return ans;
    }
    bool checkConvexHull(polygon p) { // 检测多边形是否是凸包（可以有三点共线）
        int sgn, n = p.size(), i = 0; // 如果三点共线不算凸包 去掉ccw=0的情况
        for (;; i++) { // 这一步是为了防止第一步遇到共线的三个点
            sgn = counterclockwise(p[i], p[(i + 1) % n], p[(i + 2) % n]);
            if (sgn) break;
        }
        for (; i < n; i++) {
            int ccw = counterclockwise(p[i], p[(i + 1) % n], p[(i + 2) % n]);
            if (ccw && ccw != sgn) {
                return false;
            }
        }
        return true;
    }
    polygon convexCut(polygon A, point k1, point k2) { // 半平面 k1k2 切凸包 A
        int n = A.size();  // 保留所有满足 k1 -> p -> k2 为逆时针方向的点
        A.push_back(A[0]); // 保留的点可能有重点
        polygon ans;
        line cut(k1, k2);
        for (int i = 0; i < n; i++) {
            int ccw1 = counterclockwise(k1, k2, A[i]);
            int ccw2 = counterclockwise(k1, k2, A[i + 1]);
            if (ccw1 >= 0) ans.push_back(A[i]);
            if (ccw1 * ccw2 <= 0) {
                if (sameLine(cut, line(A[i], A[i + 1]))) { // 半平面恰好切到凸包上某条边
                    ans.push_back(A[i]);
                    ans.push_back(A[i + 1]);
                } else {
                    ans.push_back(getLL(k1, k2, A[i], A[i + 1]));
                }
            }
        }
        return ans;
    }

    vector<line> getHL(vector<line>& L) { // 求半平面交 逆时针方向存储
        sort(L.begin(), L.end());
        deque<line> q;
        for (int i = 0; i < (int)L.size(); ++i) {
            if (i && sameDir(L[i], L[i - 1])) continue;
            while (q.size() > 1 && !checkpos(q[q.size() - 2], q[q.size() - 1], L[i])) q.pop_back();
            while (q.size() > 1 && !checkpos(q[1], q[0], L[i])) q.pop_front();
            q.push_back(L[i]);
        }
        while (q.size() > 2 && !checkpos(q[q.size() - 2], q[q.size() - 1], q[0])) q.pop_back();
        while (q.size() > 2 && !checkpos(q[1], q[0], q[q.size() - 1])) q.pop_front();
        vector<line> ans;
        for (int i = 0; i < q.size(); ++i) ans.push_back(q[i]);
        return ans;
    }

    db closestPoint(vector<point>& A, int l, int r) { // 最近点对, 先要按照 x 坐标排序
        if (r - l <= 5) {
            db ans = 1e20;
            for (int i = l; i <= r; ++i)
                for (int j = i + 1; j <= r; j++)
                    ans = min(ans, A[i].dis(A[j]));
            return ans;
        }
        int mid = l + r >> 1;
        db ans = min(closestPoint(A, l, mid), closestPoint(A, mid + 1, r));
        vector<point> B;
        for (int i = l; i <= r; i++)
            if (abs(A[i].x - A[mid].x) <= ans)
                B.push_back(A[i]);
        sort(B.begin(), B.end(), [&](const point& k1, const point& k2) {
            return k1.y < k2.y;
        });
        for (int i = 0; i < B.size(); i++)
            for (int j = i + 1; j < B.size() && B[j].y - B[i].y < ans; j++)
                ans = min(ans, B[i].dis(B[j]));
        return ans;
    }
}
using namespace Geometry;