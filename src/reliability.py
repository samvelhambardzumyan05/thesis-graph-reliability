"""Network reliability (all-terminal): probability that the random subgraph
(each edge works independently with probability p) stays connected.

Three methods:
  * exact_enumeration  - sum over all 2^m edge states
  * exact_factoring    - deletion-contraction  R(G) = p R(G/e) + (1-p) R(G-e)
  * monte_carlo        - simulation with a normal-approximation 95% interval
"""
import itertools
import math
import random


class DSU:
    """Disjoint-set union (union-find) with path compression."""

    def __init__(self, n):
        self.p = list(range(n))
        self.comps = n

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb
            self.comps -= 1


def _index(G):
    nodes = list(G.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    edges = [(idx[u], idx[v]) for u, v in G.edges()]
    return len(nodes), edges


def _connected(n, edges):
    d = DSU(n)
    for a, b in edges:
        d.union(a, b)
    return d.comps == 1


def exact_enumeration(G, p):
    n, edges = _index(G)
    m = len(edges)
    total = 0.0
    for state in itertools.product((0, 1), repeat=m):
        d = DSU(n)
        k = 0
        for s, (a, b) in zip(state, edges):
            if s:
                d.union(a, b)
                k += 1
        if d.comps == 1:
            total += p**k * (1 - p) ** (m - k)
    return total


def exact_factoring(G, p):
    """Deletion-contraction with simple pruning (multigraph edge list, loops dropped)."""
    n, edges = _index(G)

    def rec(n, edges):
        if n == 1:
            return 1.0
        if not _connected(n, edges):  # even with all edges working -> disconnected
            return 0.0
        e = edges[0]
        rest = edges[1:]
        a, b = e
        # deletion: remove e
        r_del = rec(n, rest)
        # contraction: merge b into a, relabel so labels stay 0..n-2
        new = []
        for u, v in rest:
            u = a if u == b else u
            v = a if v == b else v
            if u == v:
                continue  # loop
            new.append((u, v))
        # relabel to keep vertex ids compact
        used = sorted(set(range(n)) - {b})
        mp = {x: i for i, x in enumerate(used)}
        new = [(mp[u], mp[v]) for u, v in new]
        r_con = rec(n - 1, new)
        return p * r_con + (1 - p) * r_del

    return rec(n, edges)


def monte_carlo(G, p, trials=100_000, seed=0):
    """Returns (estimate, half_width_of_95%_CI)."""
    rng = random.Random(seed)
    n, edges = _index(G)
    ok = 0
    for _ in range(trials):
        d = DSU(n)
        for a, b in edges:
            if rng.random() < p:
                d.union(a, b)
        if d.comps == 1:
            ok += 1
    est = ok / trials
    hw = 1.96 * math.sqrt(max(est * (1 - est), 1e-12) / trials)
    return est, hw


# closed forms used for validation --------------------------------------------
def reliability_tree(n, p):
    return p ** (n - 1)


def reliability_cycle(n, p):
    return p**n + n * p ** (n - 1) * (1 - p)


if __name__ == "__main__":
    import networkx as nx

    G = nx.cycle_graph(6)
    G.add_edges_from([(0, 3), (1, 4)])
    for p in (0.5, 0.8, 0.95):
        print(p, exact_enumeration(G, p), exact_factoring(G, p), monte_carlo(G, p, 20000))
