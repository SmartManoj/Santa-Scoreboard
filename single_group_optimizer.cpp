// Single Group Optimizer - Optimizes only one group (n value) specified by environment variable
// Compile: g++ -O3 -march=native -std=c++17 -fopenmp -o single_group_optimizer single_group_optimizer.cpp

#include <bits/stdc++.h>
#include <omp.h>
#include <cstdlib>
using namespace std;

constexpr int MAX_N = 200;
constexpr int NV = 15;
constexpr double PI = 3.14159265358979323846;

alignas(64) const long double TX[NV] = {0,0.125,0.0625,0.2,0.1,0.35,0.075,0.075,-0.075,-0.075,-0.35,-0.1,-0.2,-0.0625,-0.125};
alignas(64) const long double TY[NV] = {0.8,0.5,0.5,0.25,0.25,0,0,-0.2,-0.2,0,0,0.25,0.25,0.5,0.5};

struct FastRNG {
    uint64_t s[2];
    FastRNG(uint64_t seed = 42) {
        s[0] = seed ^ 0x853c49e6748fea9bULL;
        s[1] = (seed * 0x9e3779b97f4a7c15ULL) ^ 0xc4ceb9fe1a85ec53ULL;
    }
    inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    inline uint64_t next() {
        uint64_t s0 = s[0], s1 = s[1], r = s0 + s1;
        s1 ^= s0; s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); s[1] = rotl(s1, 37);
        return r;
    }
    inline long double rf() { return (next() >> 11) * 0x1.0p-53L; }
    inline long double rf2() { return rf() * 2.0L - 1.0L; }
    inline int ri(int n) { return next() % n; }
    inline long double gaussian() {
        long double u1 = rf() + 1e-10L, u2 = rf();
        return sqrtl(-2.0L * logl(u1)) * cosl(2.0L * PI * u2);
    }
};

struct Poly {
    long double px[NV], py[NV];
    long double x0, y0, x1, y1;
};

inline void getPoly(long double cx, long double cy, long double deg, Poly& q) {
    long double rad = deg * (PI / 180.0L);
    long double s = sinl(rad), c = cosl(rad);
    long double minx = 1e9L, miny = 1e9L, maxx = -1e9L, maxy = -1e9L;
    for (int i = 0; i < NV; i++) {
        long double x = TX[i] * c - TY[i] * s + cx;
        long double y = TX[i] * s + TY[i] * c + cy;
        q.px[i] = x; q.py[i] = y;
        if (x < minx) minx = x; if (x > maxx) maxx = x;
        if (y < miny) miny = y; if (y > maxy) maxy = y;
    }
    q.x0 = minx; q.y0 = miny; q.x1 = maxx; q.y1 = maxy;
}

inline bool pip(long double px, long double py, const Poly& q) {
    bool in = false;
    int j = NV - 1;
    for (int i = 0; i < NV; i++) {
        if ((q.py[i] > py) != (q.py[j] > py) &&
            px < (q.px[j] - q.px[i]) * (py - q.py[i]) / (q.py[j] - q.py[i]) + q.px[i])
            in = !in;
        j = i;
    }
    return in;
}

inline bool segInt(long double ax, long double ay, long double bx, long double by,
                   long double cx, long double cy, long double dx, long double dy) {
    long double d1 = (dx-cx)*(ay-cy) - (dy-cy)*(ax-cx);
    long double d2 = (dx-cx)*(by-cy) - (dy-cy)*(bx-cx);
    long double d3 = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
    long double d4 = (bx-ax)*(dy-ay) - (by-ay)*(dx-ax);
    return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0));
}

inline bool overlap(const Poly& a, const Poly& b) {
    if (a.x1 < b.x0 || b.x1 < a.x0 || a.y1 < b.y0 || b.y1 < a.y0) return false;
    for (int i = 0; i < NV; i++) {
        if (pip(a.px[i], a.py[i], b)) return true;
        if (pip(b.px[i], b.py[i], a)) return true;
    }
    for (int i = 0; i < NV; i++) {
        int ni = (i + 1) % NV;
        for (int j = 0; j < NV; j++) {
            int nj = (j + 1) % NV;
            if (segInt(a.px[i], a.py[i], a.px[ni], a.py[ni],
                      b.px[j], b.py[j], b.px[nj], b.py[nj])) return true;
        }
    }
    return false;
}

struct Cfg {
    int n;
    long double x[MAX_N], y[MAX_N], a[MAX_N];
    Poly pl[MAX_N];
    long double gx0, gy0, gx1, gy1;

    inline void upd(int i) { getPoly(x[i], y[i], a[i], pl[i]); }
    inline void updAll() { for (int i = 0; i < n; i++) upd(i); updGlobal(); }

    inline void updGlobal() {
        gx0 = gy0 = 1e9L; gx1 = gy1 = -1e9L;
        for (int i = 0; i < n; i++) {
            if (pl[i].x0 < gx0) gx0 = pl[i].x0;
            if (pl[i].x1 > gx1) gx1 = pl[i].x1;
            if (pl[i].y0 < gy0) gy0 = pl[i].y0;
            if (pl[i].y1 > gy1) gy1 = pl[i].y1;
        }
    }

    inline bool hasOvl(int i) const {
        for (int j = 0; j < n; j++)
            if (i != j && overlap(pl[i], pl[j])) return true;
        return false;
    }

    inline bool anyOvl() const {
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (overlap(pl[i], pl[j])) return true;
        return false;
    }

    inline long double side() const { return max(gx1 - gx0, gy1 - gy0); }
    inline long double score() const { long double s = side(); return s * s / n; }

    void getBoundary(vector<int>& b) const {
        b.clear();
        long double eps = 0.01L;
        for (int i = 0; i < n; i++) {
            if (pl[i].x0 - gx0 < eps || gx1 - pl[i].x1 < eps ||
                pl[i].y0 - gy0 < eps || gy1 - pl[i].y1 < eps)
                b.push_back(i);
        }
    }
};

// Squeeze
Cfg squeeze(Cfg c) {
    long double cx = (c.gx0 + c.gx1) / 2.0L, cy = (c.gy0 + c.gy1) / 2.0L;
    for (long double scale = 0.9995L; scale >= 0.98L; scale -= 0.0005L) {
        Cfg trial = c;
        for (int i = 0; i < c.n; i++) {
            trial.x[i] = cx + (c.x[i] - cx) * scale;
            trial.y[i] = cy + (c.y[i] - cy) * scale;
        }
        trial.updAll();
        if (!trial.anyOvl()) c = trial;
        else break;
    }
    return c;
}

// Compaction
Cfg compaction(Cfg c, int iters) {
    long double bs = c.side();
    for (int it = 0; it < iters; it++) {
        long double cx = (c.gx0 + c.gx1) / 2.0L, cy = (c.gy0 + c.gy1) / 2.0L;
        bool improved = false;
        for (int i = 0; i < c.n; i++) {
            long double ox = c.x[i], oy = c.y[i];
            long double dx = cx - c.x[i], dy = cy - c.y[i];
            long double d = sqrtl(dx*dx + dy*dy);
            if (d < 1e-6L) continue;
            for (long double step : {0.02L, 0.008L, 0.003L, 0.001L, 0.0004L}) {
                c.x[i] = ox + dx/d * step; c.y[i] = oy + dy/d * step; c.upd(i);
                if (!c.hasOvl(i)) {
                    c.updGlobal();
                    if (c.side() < bs - 1e-12L) { bs = c.side(); improved = true; ox = c.x[i]; oy = c.y[i]; }
                    else { c.x[i] = ox; c.y[i] = oy; c.upd(i); }
                } else { c.x[i] = ox; c.y[i] = oy; c.upd(i); }
            }
        }
        c.updGlobal();
        if (!improved) break;
    }
    return c;
}

// Local search
Cfg localSearch(Cfg c, int maxIter) {
    long double bs = c.side();
    const long double steps[] = {0.01L, 0.004L, 0.0015L, 0.0006L, 0.00025L, 0.0001L};
    const long double rots[] = {5.0L, 2.0L, 0.8L, 0.3L, 0.1L};
    const int dx[] = {1,-1,0,0,1,1,-1,-1};
    const int dy[] = {0,0,1,-1,1,-1,1,-1};

    for (int iter = 0; iter < maxIter; iter++) {
        bool improved = false;
        for (int i = 0; i < c.n; i++) {
            long double cx = (c.gx0 + c.gx1) / 2.0L, cy = (c.gy0 + c.gy1) / 2.0L;
            long double ddx = cx - c.x[i], ddy = cy - c.y[i];
            long double dist = sqrtl(ddx*ddx + ddy*ddy);
            if (dist > 1e-6L) {
                for (long double st : steps) {
                    long double ox = c.x[i], oy = c.y[i];
                    c.x[i] += ddx/dist * st; c.y[i] += ddy/dist * st; c.upd(i);
                    if (!c.hasOvl(i)) { c.updGlobal(); if (c.side() < bs - 1e-12L) { bs = c.side(); improved = true; }
                        else { c.x[i]=ox; c.y[i]=oy; c.upd(i); c.updGlobal(); } }
                    else { c.x[i]=ox; c.y[i]=oy; c.upd(i); }
                }
            }
            for (long double st : steps) {
                for (int d = 0; d < 8; d++) {
                    long double ox=c.x[i], oy=c.y[i];
                    c.x[i] += dx[d]*st; c.y[i] += dy[d]*st; c.upd(i);
                    if (!c.hasOvl(i)) { c.updGlobal(); if (c.side() < bs - 1e-12L) { bs = c.side(); improved = true; }
                        else { c.x[i]=ox; c.y[i]=oy; c.upd(i); c.updGlobal(); } }
                    else { c.x[i]=ox; c.y[i]=oy; c.upd(i); }
                }
            }
            for (long double rt : rots) {
                for (long double da : {rt, -rt}) {
                    long double oa = c.a[i]; c.a[i] += da;
                    while (c.a[i] < 0) c.a[i] += 360.0L;
                    while (c.a[i] >= 360.0L) c.a[i] -= 360.0L;
                    c.upd(i);
                    if (!c.hasOvl(i)) { c.updGlobal(); if (c.side() < bs - 1e-12L) { bs = c.side(); improved = true; }
                        else { c.a[i]=oa; c.upd(i); c.updGlobal(); } }
                    else { c.a[i]=oa; c.upd(i); }
                }
            }
        }
        if (!improved) break;
    }
    return c;
}

// Swap move operator
bool swapTrees(Cfg& c, int i, int j) {
    if (i == j || i >= c.n || j >= c.n) return false;
    swap(c.x[i], c.x[j]);
    swap(c.y[i], c.y[j]);
    swap(c.a[i], c.a[j]);
    c.upd(i);
    c.upd(j);
    return !c.hasOvl(i) && !c.hasOvl(j);
}

// SA optimization (Enhanced with swap moves)
Cfg sa_opt(Cfg c, int iter, long double T0, long double Tm, uint64_t seed) {
    FastRNG rng(seed);
    Cfg best = c, cur = c;
    long double bs = best.side(), cs = bs, T = T0;
    long double alpha = powl(Tm / T0, 1.0L / iter);
    int noImp = 0;

    for (int it = 0; it < iter; it++) {
        int mt = rng.ri(11);
        long double sc = T / T0;
        bool valid = true;

        if (mt == 0) {
            int i = rng.ri(c.n);
            long double ox = cur.x[i], oy = cur.y[i];
            cur.x[i] += rng.gaussian() * 0.5L * sc;
            cur.y[i] += rng.gaussian() * 0.5L * sc;
            cur.upd(i);
            if (cur.hasOvl(i)) { cur.x[i]=ox; cur.y[i]=oy; cur.upd(i); valid=false; }
        }
        else if (mt == 1) {
            int i = rng.ri(c.n);
            long double ox = cur.x[i], oy = cur.y[i];
            long double bcx = (cur.gx0+cur.gx1)/2.0L, bcy = (cur.gy0+cur.gy1)/2.0L;
            long double dx = bcx - cur.x[i], dy = bcy - cur.y[i];
            long double d = sqrtl(dx*dx + dy*dy);
            if (d > 1e-6L) { cur.x[i] += dx/d * rng.rf() * 0.6L * sc; cur.y[i] += dy/d * rng.rf() * 0.6L * sc; }
            cur.upd(i);
            if (cur.hasOvl(i)) { cur.x[i]=ox; cur.y[i]=oy; cur.upd(i); valid=false; }
        }
        else if (mt == 2) {
            int i = rng.ri(c.n);
            long double oa = cur.a[i];
            cur.a[i] += rng.gaussian() * 80.0L * sc;
            while (cur.a[i] < 0) cur.a[i] += 360.0L;
            while (cur.a[i] >= 360.0L) cur.a[i] -= 360.0L;
            cur.upd(i);
            if (cur.hasOvl(i)) { cur.a[i]=oa; cur.upd(i); valid=false; }
        }
        else if (mt == 3) {
            int i = rng.ri(c.n);
            long double ox=cur.x[i], oy=cur.y[i], oa=cur.a[i];
            cur.x[i] += rng.rf2() * 0.5L * sc;
            cur.y[i] += rng.rf2() * 0.5L * sc;
            cur.a[i] += rng.rf2() * 60.0L * sc;
            while (cur.a[i] < 0) cur.a[i] += 360.0L;
            while (cur.a[i] >= 360.0L) cur.a[i] -= 360.0L;
            cur.upd(i);
            if (cur.hasOvl(i)) { cur.x[i]=ox; cur.y[i]=oy; cur.a[i]=oa; cur.upd(i); valid=false; }
        }
        else if (mt == 4) {
            vector<int> boundary; cur.getBoundary(boundary);
            if (!boundary.empty()) {
                int i = boundary[rng.ri(boundary.size())];
                long double ox=cur.x[i], oy=cur.y[i], oa=cur.a[i];
                long double bcx = (cur.gx0+cur.gx1)/2.0L, bcy = (cur.gy0+cur.gy1)/2.0L;
                long double dx = bcx - cur.x[i], dy = bcy - cur.y[i];
                long double d = sqrtl(dx*dx + dy*dy);
                if (d > 1e-6L) { cur.x[i] += dx/d * rng.rf() * 0.7L * sc; cur.y[i] += dy/d * rng.rf() * 0.7L * sc; }
                cur.a[i] += rng.rf2() * 50.0L * sc;
                while (cur.a[i] < 0) cur.a[i] += 360.0L;
                while (cur.a[i] >= 360.0L) cur.a[i] -= 360.0L;
                cur.upd(i);
                if (cur.hasOvl(i)) { cur.x[i]=ox; cur.y[i]=oy; cur.a[i]=oa; cur.upd(i); valid=false; }
            } else valid = false;
        }
        else if (mt == 5) {
            long double factor = 1.0L - rng.rf() * 0.004L * sc;
            long double cx = (cur.gx0 + cur.gx1) / 2.0L, cy = (cur.gy0 + cur.gy1) / 2.0L;
            Cfg trial = cur;
            for (int i = 0; i < c.n; i++) { trial.x[i] = cx + (cur.x[i] - cx) * factor; trial.y[i] = cy + (cur.y[i] - cy) * factor; }
            trial.updAll();
            if (!trial.anyOvl()) cur = trial; else valid = false;
        }
        else if (mt == 6) {
            int i = rng.ri(c.n);
            long double ox=cur.x[i], oy=cur.y[i];
            long double levy = powl(rng.rf() + 0.001L, -1.3L) * 0.008L;
            cur.x[i] += rng.rf2() * levy; cur.y[i] += rng.rf2() * levy; cur.upd(i);
            if (cur.hasOvl(i)) { cur.x[i]=ox; cur.y[i]=oy; cur.upd(i); valid=false; }
        }
        else if (mt == 7 && c.n > 1) {
            int i = rng.ri(c.n), j = (i + 1) % c.n;
            long double oxi=cur.x[i], oyi=cur.y[i], oxj=cur.x[j], oyj=cur.y[j];
            long double dx = rng.rf2() * 0.3L * sc, dy = rng.rf2() * 0.3L * sc;
            cur.x[i]+=dx; cur.y[i]+=dy; cur.x[j]+=dx; cur.y[j]+=dy;
            cur.upd(i); cur.upd(j);
            if (cur.hasOvl(i) || cur.hasOvl(j)) { cur.x[i]=oxi; cur.y[i]=oyi; cur.x[j]=oxj; cur.y[j]=oyj; cur.upd(i); cur.upd(j); valid=false; }
        }
        else if (mt == 10 && c.n > 1) {
            int i = rng.ri(c.n), j = rng.ri(c.n);
            Cfg old = cur;
            if (!swapTrees(cur, i, j)) {
                cur = old;
                valid = false;
            }
        }
        else {
            int i = rng.ri(c.n);
            long double ox=cur.x[i], oy=cur.y[i];
            cur.x[i] += rng.rf2() * 0.002L; cur.y[i] += rng.rf2() * 0.002L; cur.upd(i);
            if (cur.hasOvl(i)) { cur.x[i]=ox; cur.y[i]=oy; cur.upd(i); valid=false; }
        }

        if (!valid) { noImp++; T *= alpha; if (T < Tm) T = Tm; continue; }

        cur.updGlobal();
        long double ns = cur.side();
        long double delta = ns - cs;

        if (delta < 0 || rng.rf() < expl(-delta / T)) {
            cs = ns;
            if (ns < bs) { bs = ns; best = cur; noImp = 0; }
            else noImp++;
        } else { cur = best; cs = bs; noImp++; }

        if (noImp > 200) { T = min(T * 5.0L, T0); noImp = 0; }
        T *= alpha;
        if (T < Tm) T = Tm;
    }
    return best;
}

// Perturb
Cfg perturb(Cfg c, long double str, FastRNG& rng) {
    Cfg original = c;
    int np = max(1, (int)(c.n * 0.08L + str * 3.0L));
    for (int k = 0; k < np; k++) {
        int i = rng.ri(c.n);
        c.x[i] += rng.gaussian() * str * 0.5L;
        c.y[i] += rng.gaussian() * str * 0.5L;
        c.a[i] += rng.gaussian() * 30.0L;
        while (c.a[i] < 0) c.a[i] += 360.0L;
        while (c.a[i] >= 360.0L) c.a[i] -= 360.0L;
    }
    c.updAll();
    for (int iter = 0; iter < 150; iter++) {
        bool fixed = true;
        for (int i = 0; i < c.n; i++) {
            if (c.hasOvl(i)) {
                fixed = false;
                long double cx = (c.gx0+c.gx1)/2.0L, cy = (c.gy0+c.gy1)/2.0L;
                long double dx = c.x[i] - cx, dy = c.y[i] - cy;
                long double d = sqrtl(dx*dx + dy*dy);
                if (d > 1e-6L) { c.x[i] += dx/d*0.02L; c.y[i] += dy/d*0.02L; }
                c.a[i] += rng.rf2() * 15.0L;
                while (c.a[i] < 0) c.a[i] += 360.0L;
                while (c.a[i] >= 360.0L) c.a[i] -= 360.0L;
                c.upd(i);
            }
        }
        if (fixed) break;
    }
    c.updGlobal();
    if (c.anyOvl()) return original;
    return c;
}

// PARALLEL optimization
Cfg optimizeParallel(Cfg c, int iters, int restarts) {
    Cfg globalBest = c;
    long double globalBestSide = c.side();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        FastRNG rng(42 + tid * 1000 + c.n);
        Cfg localBest = c;
        long double localBestSide = c.side();

        #pragma omp for schedule(dynamic)
        for (int r = 0; r < restarts; r++) {
            Cfg start;
            if (r == 0) {
                start = c;
            }
            else if (r % 4 == 0 && r < restarts / 2) {
                start = c;
                long double angleOffset = (r / 4) * 45.0L;
                for (int i = 0; i < start.n; i++) {
                    start.a[i] += angleOffset;
                    while (start.a[i] >= 360.0L) start.a[i] -= 360.0L;
                }
                start.updAll();
                if (start.anyOvl()) {
                    start = perturb(c, 0.02L + 0.02L * (r % 8), rng);
                    if (start.anyOvl()) continue;
                }
            }
            else {
                start = perturb(c, 0.02L + 0.02L * (r % 8), rng);
                if (start.anyOvl()) continue;
            }

            uint64_t seed = 42 + r * 1000 + tid * 100000 + c.n;
            Cfg o = sa_opt(start, iters, 3.0L, 0.0000005L, seed);
            o = squeeze(o);
            o = compaction(o, 50);
            o = localSearch(o, 80);

            if (!o.anyOvl() && o.side() < localBestSide) {
                localBestSide = o.side();
                localBest = o;
            }
        }

        #pragma omp critical
        {
            if (!localBest.anyOvl() && localBestSide < globalBestSide) {
                globalBestSide = localBestSide;
                globalBest = localBest;
            }
        }
    }

    globalBest = squeeze(globalBest);
    globalBest = compaction(globalBest, 80);
    globalBest = localSearch(globalBest, 150);

    if (globalBest.anyOvl()) return c;
    return globalBest;
}

map<int, Cfg> loadCSV(const string& fn) {
    map<int, Cfg> cfg;
    ifstream f(fn);
    if (!f) return cfg;
    string ln; getline(f, ln);
    map<int, vector<tuple<int,long double,long double,long double>>> data;
    while (getline(f, ln)) {
        size_t p1=ln.find(','), p2=ln.find(',',p1+1), p3=ln.find(',',p2+1);
        string id=ln.substr(0,p1), xs=ln.substr(p1+1,p2-p1-1), ys=ln.substr(p2+1,p3-p2-1), ds=ln.substr(p3+1);
        if(!xs.empty() && xs[0]=='s') xs=xs.substr(1);
        if(!ys.empty() && ys[0]=='s') ys=ys.substr(1);
        if(!ds.empty() && ds[0]=='s') ds=ds.substr(1);
        int n=stoi(id.substr(0,3)), idx=stoi(id.substr(4));
        data[n].push_back({idx, stold(xs), stold(ys), stold(ds)});
    }
    for (auto& [n,v] : data) {
        Cfg c; c.n = n;
        for (auto& [i,x,y,d] : v) if (i < n) { c.x[i]=x; c.y[i]=y; c.a[i]=d; }
        c.updAll();
        cfg[n] = c;
    }
    return cfg;
}

void saveCSV(const string& fn, const map<int, Cfg>& cfg) {
    ofstream f(fn);
    f << fixed << setprecision(17) << "id,x,y,deg\n";
    for (int n = 1; n <= 200; n++) {
        if (cfg.count(n)) {
            const Cfg& c = cfg.at(n);
            for (int i = 0; i < n; i++)
                f << setfill('0') << setw(3) << n << "_" << i << ",s" << c.x[i] << ",s" << c.y[i] << ",s" << c.a[i] << "\n";
        }
    }
}

int main(int argc, char** argv) {
    string in="submission.csv", out="submission.csv";

    // Get group number from environment variable
    const char* groupEnv = getenv("GROUP_NUMBER");
    if (!groupEnv) {
        printf("Error: GROUP_NUMBER environment variable not set\n");
        printf("Usage: export GROUP_NUMBER=<n> && ./single_group_optimizer\n");
        return 1;
    }
    int targetN = stoi(groupEnv);
    if (targetN < 1 || targetN > 200) {
        printf("Error: GROUP_NUMBER must be between 1 and 200, got %d\n", targetN);
        return 1;
    }

    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a=="-i" && i+1<argc) in=argv[++i];
        else if (a=="-o" && i+1<argc) out=argv[++i];
    }

    int numThreads = omp_get_max_threads();
    printf("Single Group Optimizer (%d threads)\n", numThreads);
    printf("Target group: n=%d\n", targetN);
    printf("Loading %s...\n", in.c_str());

    auto cfg = loadCSV(in);
    if (cfg.empty()) { printf("No data!\n"); return 1; }
    printf("Loaded %d configs\n", (int)cfg.size());

    if (!cfg.count(targetN)) {
        printf("Error: Group n=%d not found in input file\n", targetN);
        return 1;
    }

    auto globalStart = chrono::high_resolution_clock::now();

    Cfg bestOverall = cfg[targetN];
    long double bestOverallScore = bestOverall.score();
    long double initialScore = bestOverallScore;

    printf("Initial score for n=%d: %.12Lf\n", targetN, initialScore);
    if (bestOverall.anyOvl()) printf("WARNING: Initial config has overlaps!\n");

    int noImprovementLevels = 0;
    const int maxNoImprovementLevels = 5;  // Break after 5 consecutive iters levels with no improvement
    int totalIterations = 0;
    bool converged = false;

    for (int iters = 100; iters < 10000 && !converged; iters += 100) {
        bool improvedThisLevel = false;
        int noImprovementRestarts = 0;
        const int maxNoImprovementRestarts = 10;

        for (int restarts = 10; restarts < 1000; restarts += 10) {
            totalIterations++;

            auto t0 = chrono::high_resolution_clock::now();

            printf("\n[%d] iters=%d, restarts=%d\n", totalIterations, iters, restarts);

            Cfg c = bestOverall;  // Start from best known
            Cfg o = optimizeParallel(c, iters, max(4, restarts));

            bool o_ovl = o.anyOvl();
            bool c_ovl = c.anyOvl();

            // Keep better solution
            if (!c_ovl && o_ovl) {
                o = c;
            } else if (!c_ovl && !o_ovl && o.side() > c.side() + 1e-14L) {
                o = c;
            } else if (c_ovl && o_ovl && o.side() > c.side() + 1e-14L) {
                o = c;
            }

            long double newScore = o.score();

            auto t1 = chrono::high_resolution_clock::now();
            long double el = chrono::duration_cast<chrono::milliseconds>(t1-t0).count() / 1000.0L;

            // Check for improvement
            if (newScore < bestOverallScore - 1e-12L) {
                long double improvement = (bestOverallScore - newScore) / bestOverallScore * 100.0L;
                printf("NEW BEST: %.12Lf -> %.12Lf (%.4Lf%% improvement) [%.1Lfs]\n",
                       bestOverallScore, newScore, improvement, el);
                bestOverallScore = newScore;
                bestOverall = o;
                improvedThisLevel = true;
                noImprovementRestarts = 0;

                // Save intermediate result
                cfg[targetN] = bestOverall;
                saveCSV(out, cfg);
            } else {
                noImprovementRestarts++;
                printf("No improvement (%d/%d) score=%.12Lf [%.1Lfs]\n",
                       noImprovementRestarts, maxNoImprovementRestarts, newScore, el);
                if (noImprovementRestarts >= maxNoImprovementRestarts) {
                    printf("Moving to next iters level...\n");
                    break;
                }
            }
        }

        // Track per-iters-level improvement
        if (improvedThisLevel) {
            noImprovementLevels = 0;
            printf("\n=== iters=%d: IMPROVED ===\n", iters);
        } else {
            noImprovementLevels++;
            printf("\n=== iters=%d: No improvement (%d/%d levels) ===\n",
                   iters, noImprovementLevels, maxNoImprovementLevels);
            if (noImprovementLevels >= maxNoImprovementLevels) {
                printf("\nConverged: No improvement for %d consecutive iters levels\n", maxNoImprovementLevels);
                converged = true;
            }
        }
    }

    auto globalEnd = chrono::high_resolution_clock::now();
    long double totalTime = chrono::duration_cast<chrono::milliseconds>(globalEnd-globalStart).count() / 1000.0L;

    printf("\n========================================\n");
    printf("FINAL RESULT for n=%d\n", targetN);
    printf("Initial: %.12Lf\n", initialScore);
    printf("Final:   %.12Lf\n", bestOverallScore);
    if (bestOverallScore < initialScore - 1e-10L) {
        printf("Total improvement: %.4Lf%%\n", (initialScore - bestOverallScore) / initialScore * 100.0L);
    } else {
        printf("No overall improvement\n");
    }
    printf("Total iterations: %d\n", totalIterations);
    printf("Total time: %.1Lfs (with %d threads)\n", totalTime, numThreads);
    printf("========================================\n");

    cfg[targetN] = bestOverall;
    saveCSV(out, cfg);
    printf("Saved to %s\n", out.c_str());
    return 0;
}

