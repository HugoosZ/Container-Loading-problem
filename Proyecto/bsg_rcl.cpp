// ============================================================================
// BSG-RCL: Beam Search Greedy con Restricted Candidate List para el SCLP
// Proyecto CIT3352 - Algoritmos Exactos y Metaheuristicas
//
// FUENTE DEL ALGORITMO BASE:
//   El diseño del algoritmo base (BSG-CLP) proviene del paper:
//   I. Araya, M.-C. Riff, "A beam search approach to the container loading
//   problem", Computers & Operations Research 43 (2014) 100-107.
//   Esta implementacion fue escrita desde cero a partir del pseudocodigo
//   del paper (Algoritmos 1, 2 y 3, y elementos K1-K6).
//
// CONTRIBUCION PROPIA:
//   Operador Expand-RCL: la seleccion determinista de los w mejores bloques
//   (por f(b,r)) se generaliza a una seleccion estocastica sesgada tipo
//   GRASP: se construye una lista restringida de candidatos (RCL) con los
//   K = ceil(gamma*w) mejores bloques y se muestrean w sin reemplazo con
//   pesos lineales por rango. Con gamma = 1 el algoritmo se reduce
//   EXACTAMENTE al BSG-CLP original (baseline determinista).
//
// Compilar:  g++ -O2 -std=c++17 -o bsg_rcl bsg_rcl.cpp
// Uso:       ./bsg_rcl <archivo_instancias> [opciones]   (ver --help)
// ============================================================================
#include <bits/stdc++.h>
#include <cassert>
using namespace std;
typedef long long ll;
typedef chrono::steady_clock Clock;

static const ll NEG = LLONG_MIN / 4;

// ----------------------------------------------------------------------------
// 1. INSTANCIA (parser formato BR / OR-Library thpack)
// ----------------------------------------------------------------------------
struct BoxType {
    int d[3];     // dimensiones originales (d1,d2,d3)
    bool f[3];    // flag[i]=true si la dimension i puede ir en el eje vertical
    int count;    // cantidad disponible
    ll vol;       // volumen de una caja
};
struct Orient { int l, w, h; };  // una orientacion concreta (largo,ancho,alto)

struct Instance {
    int id = 0; long seed = 0;
    int L = 0, W = 0, H = 0;            // dimensiones del contenedor
    vector<BoxType> types;
    vector<vector<Orient>> orients;     // orientaciones factibles por tipo
    ll cvol = 0;
};

// Orientaciones factibles: la dimension colocada en el eje vertical (h)
// debe tener flag=1. Horizontal siempre permitido (semantica BR).
static vector<Orient> genOrients(const BoxType& t) {
    static const int P[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    vector<Orient> res;
    for (auto& p : P) {
        if (!t.f[p[2]]) continue;                     // restriccion de orientacion
        Orient o{t.d[p[0]], t.d[p[1]], t.d[p[2]]};
        bool dup = false;
        for (auto& q : res) if (q.l==o.l && q.w==o.w && q.h==o.h) { dup = true; break; }
        if (!dup) res.push_back(o);
    }
    return res;
}

// Formato BR (thpack): 1ra linea = numero de problemas P. Por problema:
//   id seed | L W H | n | n lineas: idx d1 f1 d2 f2 d3 f3 cantidad
static vector<Instance> parseFile(const string& path) {
    ifstream in(path);
    if (!in) { cerr << "ERROR: no se pudo abrir " << path << "\n"; exit(1); }
    int P; in >> P;
    vector<Instance> v;
    for (int p = 0; p < P; p++) {
        Instance I;
        in >> I.id >> I.seed >> I.L >> I.W >> I.H;
        int n; in >> n;
        I.types.resize(n);
        for (int i = 0; i < n; i++) {
            int idx, f0, f1, f2; auto& t = I.types[i];
            in >> idx >> t.d[0] >> f0 >> t.d[1] >> f1 >> t.d[2] >> f2 >> t.count;
            t.f[0] = f0; t.f[1] = f1; t.f[2] = f2;
            t.vol = (ll)t.d[0] * t.d[1] * t.d[2];
        }
        if (!in) { cerr << "ERROR: formato de instancia invalido en problema " << p+1 << "\n"; exit(1); }
        I.cvol = (ll)I.L * I.W * I.H;
        I.orients.resize(n);
        for (int i = 0; i < n; i++) I.orients[i] = genOrients(I.types[i]);
        v.push_back(move(I));
    }
    return v;
}

// ----------------------------------------------------------------------------
// 2. GEOMETRIA: cuboides y cover representation (K1)
// ----------------------------------------------------------------------------
struct Cuboid {
    int x1, y1, z1, x2, y2, z2;
    int lx() const { return x2 - x1; }
    int ly() const { return y2 - y1; }
    int lz() const { return z2 - z1; }
    ll  vol() const { return (ll)lx() * ly() * lz(); }
    bool containsC(const Cuboid& o) const {
        return x1<=o.x1 && y1<=o.y1 && z1<=o.z1 && x2>=o.x2 && y2>=o.y2 && z2>=o.z2;
    }
    bool operator==(const Cuboid& o) const {
        return x1==o.x1 && y1==o.y1 && z1==o.z1 && x2==o.x2 && y2==o.y2 && z2==o.z2;
    }
};
static bool overlapC(const Cuboid& a, const Cuboid& b) {
    return a.x1<b.x2 && b.x1<a.x2 && a.y1<b.y2 && b.y1<a.y2 && a.z1<b.z2 && b.z1<a.z2;
}

// ----------------------------------------------------------------------------
// 3. BLOQUES (K2): simples + generales (combinacion por pares)
// ----------------------------------------------------------------------------
struct Block {
    int l, w, h;
    ll vol;                          // volumen UTILIZADO (suma de cajas)
    vector<pair<int,int>> need;      // (tipo, cantidad), ordenado por tipo
    ll cub() const { return (ll)l * w * h; }
};

struct BlockGenerator {
    const Instance& I;
    int max_bl, min_fr;
    unordered_set<string> seen;
    vector<Block> out;
    BlockGenerator(const Instance& i, int mb, int mf) : I(i), max_bl(mb), min_fr(mf) {}

    string key(const Block& b) {
        string s = to_string(b.l) + "x" + to_string(b.w) + "x" + to_string(b.h) + "|";
        for (auto& p : b.need) { s += to_string(p.first) + ":" + to_string(p.second) + ","; }
        return s;
    }
    bool add(Block b) {
        if ((int)out.size() >= max_bl) return false;
        if (b.l > I.L || b.w > I.W || b.h > I.H) return false;
        sort(b.need.begin(), b.need.end());
        string k = key(b);
        if (seen.count(k)) return false;
        seen.insert(k);
        out.push_back(move(b));
        return true;
    }
    bool needFeasible(const vector<pair<int,int>>& need) {
        for (auto& p : need) if (p.second > I.types[p.first].count) return false;
        return true;
    }
    // Bloques simples: arreglos nx x ny x nz de un tipo en una orientacion
    void genSimple() {
        for (int t = 0; t < (int)I.types.size(); t++) {
            const auto& bt = I.types[t];
            for (const auto& o : I.orients[t]) {
                for (int nx = 1; nx * o.l <= I.L; nx++) {
                    if (nx > bt.count) break;
                    for (int ny = 1; ny * o.w <= I.W; ny++) {
                        if (nx * ny > bt.count) break;
                        for (int nz = 1; nz * o.h <= I.H; nz++) {
                            int q = nx * ny * nz;
                            if (q > bt.count) break;
                            Block b{nx*o.l, ny*o.w, nz*o.h, (ll)q * bt.vol, {{t, q}}};
                            add(move(b));
                            if ((int)out.size() >= max_bl) return;
                        }
                    }
                }
            }
        }
    }
    // Bloques generales: combinacion por pares a lo largo de cada eje
    // (Fanslau & Bortfeldt 2010; usado por Zhu et al. 2012 y Araya & Riff 2014).
    // Simplificacion declarada: el pool de combinacion se limita a los
    // POOL_CAP bloques de mayor volumen y a un maximo de rondas.
    void genGeneral() {
        const int POOL_CAP = 1200, MAX_ROUNDS = 3;
        for (int round = 0; round < MAX_ROUNDS && (int)out.size() < max_bl; round++) {
            vector<Block> pool = out;
            sort(pool.begin(), pool.end(), [](const Block&a, const Block&b){ return a.vol > b.vol; });
            if ((int)pool.size() > POOL_CAP) pool.resize(POOL_CAP);
            bool added = false;
            for (size_t i = 0; i < pool.size() && (int)out.size() < max_bl; i++)
                for (size_t j = 0; j < pool.size() && (int)out.size() < max_bl; j++)
                    for (int ax = 0; ax < 3; ax++) {
                        const Block &a = pool[i], &b = pool[j];
                        int nl, nw, nh;
                        if      (ax == 0) { nl = a.l + b.l; nw = max(a.w,b.w); nh = max(a.h,b.h); }
                        else if (ax == 1) { nl = max(a.l,b.l); nw = a.w + b.w; nh = max(a.h,b.h); }
                        else              { nl = max(a.l,b.l); nw = max(a.w,b.w); nh = a.h + b.h; }
                        if (nl > I.L || nw > I.W || nh > I.H) continue;
                        ll nvol = a.vol + b.vol;
                        if (nvol * 100 < (ll)min_fr * nl * nw * nh) continue;   // % relleno
                        vector<pair<int,int>> need = a.need;
                        for (auto& p : b.need) {
                            bool found = false;
                            for (auto& q : need) if (q.first == p.first) { q.second += p.second; found = true; break; }
                            if (!found) need.push_back(p);
                        }
                        if (!needFeasible(need)) continue;
                        if (add(Block{nl, nw, nh, nvol, move(need)})) added = true;
                    }
            if (!added) break;
        }
    }
    vector<Block> run() {
        genSimple();
        if (min_fr < 100) genGeneral();
        // Orden por volumen utilizado descendente: habilita la poda f(b,r) <= V(b)
        sort(out.begin(), out.end(), [](const Block&a, const Block&b){ return a.vol > b.vol; });
        if ((int)out.size() > max_bl) out.resize(max_bl);
        return move(out);
    }
};

// ----------------------------------------------------------------------------
// 4. ESTADO y colocacion de bloques
// ----------------------------------------------------------------------------
struct Placement { int block; int x, y, z; };

struct State {
    vector<Cuboid> R;        // espacio residual (cover representation)
    vector<int> rem;         // cajas restantes por tipo
    ll loaded = 0;           // volumen cargado
    vector<Placement> plan;  // plan de carga
};

// Distancia Manhattan del cuboide a su esquina ancla del contenedor (K3)
static inline int anchorDist(const Instance& I, const Cuboid& c) {
    return min(c.x1, I.L - c.x2) + min(c.y1, I.W - c.y2) + min(c.z1, I.H - c.z2);
}
// Seleccion del espacio libre: menor distancia, empate por mayor volumen
static int selectSpace(const Instance& I, const vector<Cuboid>& R) {
    int best = 0, bd = INT_MAX; ll bv = -1;
    for (int i = 0; i < (int)R.size(); i++) {
        int d = anchorDist(I, R[i]); ll v = R[i].vol();
        if (d < bd || (d == bd && v > bv)) { bd = d; bv = v; best = i; }
    }
    return best;
}
static inline bool needOK(const Block& b, const vector<int>& rem) {
    for (auto& p : b.need) if (rem[p.first] < p.second) return false;
    return true;
}

// Colocar bloque bi en la esquina ancla del cuboide ri y actualizar R (K1, K5)
static void place(const Instance& I, State& s, int bi, int ri, const vector<Block>& B) {
    const Block& b = B[bi];
    const Cuboid r = s.R[ri];
    int px = (r.x1 <= I.L - r.x2) ? r.x1 : r.x2 - b.l;   // ancla por eje
    int py = (r.y1 <= I.W - r.y2) ? r.y1 : r.y2 - b.w;
    int pz = (r.z1 <= I.H - r.z2) ? r.z1 : r.z2 - b.h;
    Cuboid reg{px, py, pz, px + b.l, py + b.w, pz + b.h};

    // Cuboides residuales (hasta 6 por cuboide solapado)
    vector<Cuboid> nR; nR.reserve(s.R.size() + 8);
    for (const auto& c : s.R) {
        if (!overlapC(c, reg)) { nR.push_back(c); continue; }
        if (reg.x1 > c.x1) nR.push_back({c.x1, c.y1, c.z1, reg.x1, c.y2, c.z2});
        if (reg.x2 < c.x2) nR.push_back({reg.x2, c.y1, c.z1, c.x2, c.y2, c.z2});
        if (reg.y1 > c.y1) nR.push_back({c.x1, c.y1, c.z1, c.x2, reg.y1, c.z2});
        if (reg.y2 < c.y2) nR.push_back({c.x1, reg.y2, c.z1, c.x2, c.y2, c.z2});
        if (reg.z1 > c.z1) nR.push_back({c.x1, c.y1, c.z1, c.x2, c.y2, reg.z1});
        if (reg.z2 < c.z2) nR.push_back({c.x1, c.y1, reg.z2, c.x2, c.y2, c.z2});
    }
    for (auto& p : b.need) s.rem[p.first] -= p.second;
    s.loaded += b.vol;
    s.plan.push_back({bi, px, py, pz});

    // Filtro de cuboides no maximales (duplicados: sobrevive el de menor indice)
    int n = nR.size();
    vector<char> dead(n, 0);
    for (int i = 0; i < n; i++) {
        if (dead[i]) continue;
        for (int j = 0; j < n; j++) {
            if (i == j || dead[j]) continue;
            if (nR[j].containsC(nR[i]) && (!(nR[i] == nR[j]) || j < i)) { dead[i] = 1; break; }
        }
    }
    // Poda: descartar cuboides donde no cabe ninguna caja restante
    s.R.clear();
    for (int i = 0; i < n; i++) {
        if (dead[i]) continue;
        bool fitsAny = false;
        for (int t = 0; t < (int)I.types.size() && !fitsAny; t++) {
            if (s.rem[t] <= 0) continue;
            for (const auto& o : I.orients[t])
                if (o.l <= nR[i].lx() && o.w <= nR[i].ly() && o.h <= nR[i].lz()) { fitsAny = true; break; }
        }
        if (fitsAny) s.R.push_back(nR[i]);
    }
}

// ----------------------------------------------------------------------------
// 5. EVALUACION (K4): f(b,r) = V(b) - Vloss(b,r), con KPA (subset-sum DP)
//    El KPA se computa UNA vez por estado y se reutiliza (igual que el paper).
// ----------------------------------------------------------------------------
static const int MAXB = 4096;   // cota superior de dimension del contenedor

struct KPA { vector<int> bx, by, bz; };  // bX[c] = mayor combinacion lineal <= c

static void addBounded(bitset<MAXB>& bs, int val, int cnt) {
    // knapsack acotado por binary splitting sobre bitset
    ll k = 1;
    while (cnt > 0) {
        ll t = min<ll>(k, cnt);
        ll sh = (ll)val * t;
        if (sh < MAXB) bs |= bs << (size_t)sh;
        cnt -= (int)t; k <<= 1;
    }
}
static vector<int> buildBest(const bitset<MAXB>& bs, int cap) {
    vector<int> best(cap + 1);
    int cur = 0;
    for (int c = 0; c <= cap; c++) { if (bs[c]) cur = c; best[c] = cur; }
    return best;
}
static KPA computeKPA(const Instance& I, const vector<int>& rem) {
    bitset<MAXB> hb, vb;   // dims horizontales (todas) / verticales (con flag)
    hb.set(0); vb.set(0);
    for (int t = 0; t < (int)I.types.size(); t++) {
        if (rem[t] <= 0) continue;
        for (int i = 0; i < 3; i++) {
            addBounded(hb, I.types[t].d[i], rem[t]);
            if (I.types[t].f[i]) addBounded(vb, I.types[t].d[i], rem[t]);
        }
    }
    KPA kp;
    kp.bx = buildBest(hb, I.L);
    kp.by = buildBest(hb, I.W);
    kp.bz = buildBest(vb, I.H);
    return kp;
}
static inline ll fval(const Block& b, const Cuboid& r, const KPA& kp) {
    int dl = r.lx() - b.l, dw = r.ly() - b.w, dh = r.lz() - b.h;
    if (dl < 0 || dw < 0 || dh < 0) return NEG;          // no cabe
    ll lm = kp.bx[dl], wm = kp.by[dw], hm = kp.bz[dh];
    ll vloss = r.vol() - (ll)(b.l + lm) * (b.w + wm) * (b.h + hm);
    return b.vol - vloss;
}

// ----------------------------------------------------------------------------
// 6. GREEDY de evaluacion (determinista, identico al paper)
// ----------------------------------------------------------------------------
struct GOut { ll vol; string sig; State sol; };

static GOut greedy(const Instance& I, State s, const vector<Block>& B,
                   Clock::time_point deadline) {
    KPA kp = computeKPA(I, s.rem);
    int it = 0;
    while (!s.R.empty()) {
        if (((++it) & 7) == 0 && Clock::now() > deadline) break;
        int ri = selectSpace(I, s.R);
        const Cuboid r = s.R[ri];
        ll best = NEG; int bi = -1;
        for (int idx = 0; idx < (int)B.size(); idx++) {
            const Block& b = B[idx];
            if (b.vol <= best) break;                    // poda: f(b,r) <= V(b)
            if (b.l > r.lx() || b.w > r.ly() || b.h > r.lz()) continue;
            if (!needOK(b, s.rem)) continue;
            ll f = fval(b, r, kp);
            if (f > best) { best = f; bi = idx; }
        }
        if (bi < 0) { s.R.erase(s.R.begin() + ri); continue; }   // espacio inutil
        place(I, s, bi, ri, B);
    }
    string sig; sig.reserve(s.rem.size() * 4);
    for (int x : s.rem) { sig += to_string(x); sig += ','; }
    return {s.loaded, move(sig), move(s)};
}

// ----------------------------------------------------------------------------
// 7. OPERADOR PROPUESTO: Expand-RCL
//    K = ceil(gamma * w) candidatos; muestreo de w sin reemplazo con pesos
//    lineales por rango. gamma = 1  ==>  top-w determinista (BSG-CLP original)
// ----------------------------------------------------------------------------
static vector<State> expandRCL(const Instance& I, const State& s0, int wEff,
                               double gamma, mt19937& rng,
                               const vector<Block>& B, Clock::time_point deadline) {
    State cur = s0;
    KPA kp = computeKPA(I, cur.rem);
    while (!cur.R.empty()) {
        if (Clock::now() > deadline) return {};
        int ri = selectSpace(I, cur.R);
        const Cuboid r = cur.R[ri];
        int K = max(1, (int)ceil(gamma * wEff));

        // top-K bloques por f(b,r) con heap minimo + poda por volumen
        priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<>> heap;
        for (int idx = 0; idx < (int)B.size(); idx++) {
            const Block& b = B[idx];
            if ((int)heap.size() == K && b.vol <= heap.top().first) break;
            if (b.l > r.lx() || b.w > r.ly() || b.h > r.lz()) continue;
            if (!needOK(b, cur.rem)) continue;
            ll f = fval(b, r, kp);
            heap.push({f, idx});
            if ((int)heap.size() > K) heap.pop();
        }
        if (heap.empty()) { cur.R.erase(cur.R.begin() + ri); continue; }

        vector<pair<ll,int>> cand;
        while (!heap.empty()) { cand.push_back(heap.top()); heap.pop(); }
        sort(cand.rbegin(), cand.rend());                 // descendente por f
        int m = min<int>(wEff, cand.size());

        // muestreo ponderado sin reemplazo: peso(rango i) = |cand| - i
        vector<int> alive(cand.size());
        iota(alive.begin(), alive.end(), 0);
        vector<int> chosen;
        for (int k = 0; k < m; k++) {
            double total = 0;
            for (int a : alive) total += (double)(cand.size() - a);
            uniform_real_distribution<double> U(0.0, total);
            double x = U(rng); int pick = alive.back(); size_t pos = alive.size() - 1;
            double acc = 0;
            for (size_t q = 0; q < alive.size(); q++) {
                acc += (double)(cand.size() - alive[q]);
                if (x <= acc) { pick = alive[q]; pos = q; break; }
            }
            chosen.push_back(pick);
            alive.erase(alive.begin() + pos);
        }
        vector<State> out; out.reserve(m);
        for (int ci : chosen) {
            State s2 = cur;
            place(I, s2, cand[ci].second, ri, B);
            out.push_back(move(s2));
        }
        return out;
    }
    return {};
}

// ----------------------------------------------------------------------------
// 8. BEAM SEARCH (K6) con remocion de estados similares
// ----------------------------------------------------------------------------
struct BestSol { ll vol = -1; State sol; };

static void beamSearch(const Instance& I, const vector<Block>& B, int w,
                       double gamma, mt19937& rng, Clock::time_point deadline,
                       BestSol& best) {
    State s0;
    s0.R = {{0, 0, 0, I.L, I.W, I.H}};
    s0.rem.resize(I.types.size());
    for (int t = 0; t < (int)I.types.size(); t++) s0.rem[t] = I.types[t].count;

    vector<State> S; S.push_back(move(s0));
    bool root = true;
    while (!S.empty()) {
        if (Clock::now() > deadline) return;
        vector<State> succ;
        int wRoot = (int)min<ll>((ll)w * w, 100000);   // evita overflow en instancias triviales
        for (auto& s : S) {
            auto v = expandRCL(I, s, root ? wRoot : w, gamma, rng, B, deadline);
            for (auto& x : v) succ.push_back(move(x));
        }
        root = false;
        if (succ.empty()) return;

        // Evaluacion greedy de cada sucesor
        struct Ev { ll g; string sig; int idx; ll loaded; };
        vector<Ev> evs; evs.reserve(succ.size());
        for (int i = 0; i < (int)succ.size(); i++) {
            if (Clock::now() > deadline) return;
            GOut go = greedy(I, succ[i], B, deadline);
            if (go.vol > best.vol) { best.vol = go.vol; best.sol = move(go.sol); }
            evs.push_back({go.vol, move(go.sig), i, succ[i].loaded});
        }
        // Remocion de similares: misma firma => sobrevive el de MENOR volumen actual
        unordered_map<string, int> keep;
        for (int e = 0; e < (int)evs.size(); e++) {
            auto it = keep.find(evs[e].sig);
            if (it == keep.end() || evs[it->second].loaded > evs[e].loaded) keep[evs[e].sig] = e;
        }
        vector<int> fil;
        for (auto& kv : keep) fil.push_back(kv.second);
        sort(fil.begin(), fil.end(), [&](int a, int b){ return evs[a].g > evs[b].g; });

        vector<State> nS;
        for (int k = 0; k < (int)fil.size() && k < w; k++)
            nS.push_back(move(succ[evs[fil[k]].idx]));
        S = move(nS);
    }
}

// ----------------------------------------------------------------------------
// 9. BUCLE EXTERNO (doble esfuerzo) y runner
// ----------------------------------------------------------------------------
struct RunResult { double fill; ll vol; int lastw; double secs; BestSol best; };

static RunResult solveOne(const Instance& I, double tlimit, double gamma,
                          unsigned seed, int min_fr, int max_bl, bool verbose) {
    auto start = Clock::now();
    auto deadline = start + chrono::duration_cast<Clock::duration>(chrono::duration<double>(tlimit));
    BlockGenerator gen(I, max_bl, min_fr);
    vector<Block> B = gen.run();
    if (verbose) cerr << "  bloques generados: " << B.size() << "\n";
    mt19937 rng(seed);
    BestSol best;
    ll totalBoxVol = 0;
    for (auto& t : I.types) totalBoxVol += t.vol * t.count;
    const ll optimum = min(totalBoxVol, I.cvol);
    int w = 1, lastw = 1;
    while (Clock::now() < deadline) {
        beamSearch(I, B, w, gamma, rng, deadline, best);
        lastw = w;
        if (best.vol >= optimum) break;            // todas las cajas cargadas: optimo
        int nw = (int)ceil(w * sqrt(2.0));
        w = (nw == w) ? w + 1 : nw;
        w = min(w, 10000);                         // tope de seguridad
        if (verbose) {
            double fill = best.vol < 0 ? 0.0 : 100.0 * best.vol / I.cvol;
            cerr << "  w=" << lastw << "  mejor=" << fixed << setprecision(2) << fill << "%\n";
        }
    }
    double secs = chrono::duration<double>(Clock::now() - start).count();
    double fill = best.vol < 0 ? 0.0 : 100.0 * best.vol / I.cvol;
    return {fill, max<ll>(best.vol, 0), lastw, secs, move(best)};
}

static void dumpSolution(const Instance& I, const vector<Block>& dummy,
                         const BestSol& best, const vector<Block>& B, ostream& os) {
    os << "# contenedor " << I.L << " " << I.W << " " << I.H << "\n";
    os << "# colocaciones: bloque(l w h) pos(x y z) cajas\n";
    for (auto& p : best.sol.plan) {
        const Block& b = B[p.block];
        os << b.l << " " << b.w << " " << b.h << "  @ "
           << p.x << " " << p.y << " " << p.z << "  [";
        for (auto& q : b.need) os << " t" << q.first + 1 << "x" << q.second;
        os << " ]\n";
    }
}

// ----------------------------------------------------------------------------
// 10. SELFTEST de geometria y end-to-end
// ----------------------------------------------------------------------------
static int selftest() {
    // (a) residuales: bloque 5x5x5 en esquina de contenedor 10x10x10
    Instance I; I.L = I.W = I.H = 10; I.cvol = 1000;
    BoxType t{}; t.d[0]=t.d[1]=t.d[2]=5; t.f[0]=t.f[1]=t.f[2]=true; t.count=8; t.vol=125;
    I.types = {t}; I.orients = {genOrients(t)};
    vector<Block> B = {{5,5,5,125,{{0,1}}}};
    State s; s.R = {{0,0,0,10,10,10}}; s.rem = {8};
    place(I, s, 0, 0, B);
    assert(s.R.size() == 3);
    for (auto& c : s.R) assert(c.vol() == 250 || c.vol() == 500);
    assert(s.loaded == 125 && s.rem[0] == 7);
    cerr << "[OK] residuales tras colocar bloque en esquina (3 cuboides maximales)\n";

    // (b) end-to-end: 8 cubos de 5 en contenedor de 10 => 100%
    RunResult r = solveOne(I, 1.0, 1.0, 1, 100, 10000, false);
    assert(fabs(r.fill - 100.0) < 1e-9);
    cerr << "[OK] end-to-end cubo perfecto: relleno = " << r.fill << "%\n";

    // (c) gamma=1 es determinista: dos corridas, mismo resultado
    Instance J; J.L=100; J.W=80; J.H=60; J.cvol=(ll)100*80*60;
    BoxType a{}; a.d[0]=30;a.d[1]=25;a.d[2]=20; a.f[0]=a.f[1]=a.f[2]=true; a.count=12; a.vol=15000;
    BoxType b2{}; b2.d[0]=22;b2.d[1]=18;b2.d[2]=15; b2.f[0]=b2.f[1]=b2.f[2]=true; b2.count=20; b2.vol=5940;
    J.types={a,b2}; J.orients={genOrients(a),genOrients(b2)};
    RunResult r1 = solveOne(J, 1.5, 1.0, 1, 100, 10000, false);
    RunResult r2 = solveOne(J, 1.5, 1.0, 99, 100, 10000, false);
    assert(r1.vol == r2.vol);   // misma solucion con distinta semilla => determinista
    cerr << "[OK] gamma=1 determinista (relleno " << fixed << setprecision(2) << r1.fill << "%)\n";
    cerr << "SELFTEST COMPLETO: todo OK\n";
    return 0;
}

// ----------------------------------------------------------------------------
// main
// ----------------------------------------------------------------------------
static void usage() {
    cerr <<
      "Uso: ./bsg_rcl <archivo_BR> [opciones]\n"
      "  --problems a-b   rango de problemas a resolver (def: 1-1)\n"
      "  --all            todos los problemas del archivo\n"
      "  --time S         limite de tiempo por instancia en segundos (def: 30)\n"
      "  --gamma G        factor RCL; 1.0 = BSG-CLP determinista (def: 1.0)\n"
      "  --seed S         semilla RNG (def: 1)\n"
      "  --min_fr F       %% relleno min. de bloques generales; 100 = solo simples\n"
      "                   (def: 100; usar 98 en BR8-BR15)\n"
      "  --max_bl N       maximo de bloques (def: 10000)\n"
      "  --csv F          anexar resultados a CSV\n"
      "  --sol F          escribir plan de carga de la mejor solucion\n"
      "  --selftest       correr pruebas unitarias y salir\n";
}

int main(int argc, char** argv) {
    string file, csv, solfile;
    double tlimit = 30, gamma = 1.0;
    unsigned seed = 1;
    int min_fr = 100, max_bl = 10000, pa = 1, pb = 1;
    bool all = false, doSelftest = false;

    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        auto next = [&](){ if (i+1 >= argc) { usage(); exit(1); } return string(argv[++i]); };
        if      (a == "--selftest") doSelftest = true;
        else if (a == "--time")     tlimit = stod(next());
        else if (a == "--gamma")    gamma = stod(next());
        else if (a == "--seed")     seed = (unsigned)stoul(next());
        else if (a == "--min_fr")   min_fr = stoi(next());
        else if (a == "--max_bl")   max_bl = stoi(next());
        else if (a == "--csv")      csv = next();
        else if (a == "--sol")      solfile = next();
        else if (a == "--all")      all = true;
        else if (a == "--problems") { string r = next(); size_t d = r.find('-');
            pa = stoi(r.substr(0, d)); pb = (d == string::npos) ? pa : stoi(r.substr(d+1)); }
        else if (a == "--help")     { usage(); return 0; }
        else if (a[0] != '-')       file = a;
        else { cerr << "opcion desconocida: " << a << "\n"; usage(); return 1; }
    }
    if (doSelftest) return selftest();
    if (file.empty()) { usage(); return 1; }

    vector<Instance> inst = parseFile(file);
    if (all) { pa = 1; pb = inst.size(); }
    pa = max(1, pa); pb = min((int)inst.size(), pb);

    bool newCsv = false;
    if (!csv.empty()) { ifstream chk(csv); newCsv = !chk.good() || chk.peek() == EOF; }
    ofstream co;
    if (!csv.empty()) {
        co.open(csv, ios::app);
        if (newCsv) co << "archivo,problema,gamma,semilla,min_fr,max_bl,tlimite,relleno_pct,vol_cargado,vol_contenedor,w_final,segundos\n";
    }

    for (int p = pa; p <= pb; p++) {
        const Instance& I = inst[p - 1];
        cerr << "== " << file << " problema " << p << "  (gamma=" << gamma
             << ", seed=" << seed << ", t=" << tlimit << "s) ==\n";
        // se regeneran bloques dentro de solveOne (1 vez por instancia)
        RunResult r = solveOne(I, tlimit, gamma, seed, min_fr, max_bl, true);
        cout << file << " p" << p << "  gamma=" << fixed << setprecision(2) << gamma
             << " seed=" << seed
             << "  relleno=" << setprecision(3) << r.fill
             << "%  w_final=" << r.lastw << "  t=" << setprecision(1) << r.secs << "s\n";
        if (co.is_open())
            co << file << "," << p << "," << fixed << setprecision(2) << gamma << ","
               << seed << "," << min_fr << ","
               << max_bl << "," << setprecision(1) << tlimit << "," << setprecision(4) << r.fill << ","
               << r.vol << "," << I.cvol << "," << r.lastw << "," << setprecision(2) << r.secs << "\n";
        if (!solfile.empty()) {
            BlockGenerator gen(I, max_bl, min_fr);
            vector<Block> B = gen.run();
            ofstream so(solfile + (pb > pa ? ("." + to_string(p)) : ""));
            dumpSolution(I, B, r.best, B, so);
        }
    }
    return 0;
}
