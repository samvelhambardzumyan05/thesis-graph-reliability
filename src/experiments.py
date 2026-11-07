"""Experiments for Chapter 3. Run: python src/experiments.py  (writes results/results.json and thesis/figures/*.pdf)"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(__file__))
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import connectivity as C
import reliability as R

ROOT = os.path.join(os.path.dirname(__file__), "..")
FIG = os.path.join(ROOT, "thesis", "figures")
out = {}

nets = {
    "Ցիկլ C8": nx.cycle_graph(8),
    "Անիվ W7": nx.wheel_graph(7),
    "Ցանց 3x3": nx.convert_node_labels_to_integers(nx.grid_2d_graph(3, 3)),
    "Լրիվ K5": nx.complete_graph(5),
    "Պետերսեն": nx.petersen_graph(),
}

# 1. structural parameters
rows = []
for name, G in nets.items():
    adj = C.from_networkx(G)
    br, cut = C.bridges_and_articulation_points(adj)
    rows.append(dict(name=name, n=G.number_of_nodes(), m=G.number_of_edges(),
                     delta=min(d for _, d in G.degree()),
                     kappa=C.vertex_connectivity(adj), lam=C.edge_connectivity(adj),
                     bridges=len(br), cut=len(cut)))
out["structure"] = rows

# 2. reliability vs p
ps = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
curves = {}
for name, G in nets.items():
    if G.number_of_edges() <= 16:
        curves[name] = [R.exact_factoring(G, p) for p in ps]
    else:
        curves[name] = [R.monte_carlo(G, p, 20000, seed=7)[0] for p in ps]
out["ps"], out["curves"] = ps, curves
plt.figure(figsize=(6, 4))
for name, ys in curves.items():
    plt.plot(ps, ys, marker="o", label=name)
plt.xlabel("p"); plt.ylabel("R(G)"); plt.grid(alpha=.3); plt.legend()
plt.tight_layout(); plt.savefig(os.path.join(FIG, "rel_vs_p.pdf")); plt.close()

# 3. method timing: enumeration vs factoring vs MC on growing graphs (wheel graphs)
timing = []
for k in range(4, 15):
    G = nx.wheel_graph(k)
    m = G.number_of_edges()
    row = dict(n=k, m=m)
    p = 0.9
    if m <= 16:
        t = time.perf_counter(); ex = R.exact_enumeration(G, p); row["t_enum"] = time.perf_counter() - t
    t = time.perf_counter(); fa = R.exact_factoring(G, p); row["t_fact"] = time.perf_counter() - t
    t = time.perf_counter(); mc, hw = R.monte_carlo(G, p, 20000, seed=3); row["t_mc"] = time.perf_counter() - t
    row.update(R_exact=fa, R_mc=mc, hw=hw)
    timing.append(row)
out["timing"] = timing
plt.figure(figsize=(6, 4))
ms = [r["m"] for r in timing]
em = [(r["m"], r["t_enum"]) for r in timing if "t_enum" in r]
plt.semilogy([a for a, _ in em], [b for _, b in em], marker="o", label="Լրիվ թվարկում")
plt.semilogy(ms, [r["t_fact"] for r in timing], marker="s", label="Գործոնավորում")
plt.semilogy(ms, [r["t_mc"] for r in timing], marker="^", label="Monte Carlo (20000)")
plt.xlabel("Կողերի թիվը m"); plt.ylabel("Ժամանակ, վրկ"); plt.grid(alpha=.3, which="both"); plt.legend()
plt.tight_layout(); plt.savefig(os.path.join(FIG, "timing.pdf")); plt.close()

# 4. MC error vs number of trials (Petersen graph, p=0.8 -> reference from large MC)
G = nets["Ցանց 3x3"]; p = 0.8
ref = R.exact_factoring(G, p)
Ns = [100, 300, 1000, 3000, 10000, 30000, 100000]
errs = []
for N in Ns:
    e = [abs(R.monte_carlo(G, p, N, seed=s)[0] - ref) for s in range(10)]
    errs.append(sum(e) / len(e))
out["mc_ref"], out["mc_N"], out["mc_err"] = ref, Ns, errs
plt.figure(figsize=(6, 4))
plt.loglog(Ns, errs, marker="o", label="միջին բացարձակ սխալ")
plt.loglog(Ns, [errs[0] * (Ns[0] / n) ** 0.5 for n in Ns], "--", label=r"$\propto 1/\sqrt{N}$")
plt.xlabel("Փորձերի թիվը N"); plt.ylabel("Սխալ"); plt.grid(alpha=.3, which="both"); plt.legend()
plt.tight_layout(); plt.savefig(os.path.join(FIG, "mc_error.pdf")); plt.close()

# 5. bridge removal effect: reliability of C8 + chords
G0 = nx.cycle_graph(8)
variants = {"C8": G0.copy()}
G1 = G0.copy(); G1.add_edge(0, 4); variants["C8 + 1 լրաց. կող"] = G1
G2 = G1.copy(); G2.add_edge(2, 6); variants["C8 + 2 լրաց. կող"] = G2
out["chords"] = {k: [R.exact_factoring(g, p) for p in (0.7, 0.9)] for k, g in variants.items()}


# 6. cut-based bounds at p = 0.9
def cut_bounds(G, p):
    """Upper: R <= 1-(1-p)^lambda.  Lower: R >= 1 - sum over vertex bipartitions (1-p)^{|cut|}."""
    nodes = list(G.nodes()); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    E = [(idx[a], idx[b]) for a, b in G.edges()]
    lam = None; s = 0.0
    for mask in range(1, 2 ** (n - 1)):  # vertex n-1 is always on the side "0"
        c = sum(1 for a, b in E if ((mask >> a) & 1) != ((mask >> b) & 1))
        lam = c if lam is None else min(lam, c)
        s += (1 - p) ** c
    return max(0.0, 1 - s), 1 - (1 - p) ** lam
bounds = []
for name, G in nets.items():
    lo, up = cut_bounds(G, 0.9)
    bounds.append(dict(name=name, lower=lo, exact=R.exact_factoring(G, 0.9), upper=up))
out["bounds"] = bounds

json.dump(out, open(os.path.join(ROOT, "results", "results.json"), "w"), ensure_ascii=False, indent=1)
print(json.dumps({k: out[k] for k in ("bounds",)}, ensure_ascii=False, indent=1))
