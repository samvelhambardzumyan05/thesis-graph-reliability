"""Network reliability: exact (enumeration) vs Monte Carlo estimate.

All-terminal reliability: probability that the graph with randomly failing
edges (each edge works with probability p) remains connected.
"""
import itertools
import random

import networkx as nx


def exact_reliability(G: nx.Graph, p: float) -> float:
    """Exact all-terminal reliability by enumerating edge states (small graphs only)."""
    edges = list(G.edges())
    total = 0.0
    for state in itertools.product([0, 1], repeat=len(edges)):
        H = nx.Graph()
        H.add_nodes_from(G.nodes())
        H.add_edges_from(e for e, s in zip(edges, state) if s)
        if nx.is_connected(H):
            k = sum(state)
            total += p**k * (1 - p) ** (len(edges) - k)
    return total


def monte_carlo_reliability(G: nx.Graph, p: float, trials: int = 100_000, seed: int = 0) -> float:
    rng = random.Random(seed)
    ok = 0
    for _ in range(trials):
        H = nx.Graph()
        H.add_nodes_from(G.nodes())
        H.add_edges_from(e for e in G.edges() if rng.random() < p)
        ok += nx.is_connected(H)
    return ok / trials


if __name__ == "__main__":
    G = nx.cycle_graph(6)
    G.add_edges_from([(0, 3), (1, 4)])
    for p in (0.5, 0.8, 0.95):
        print(p, exact_reliability(G, p), monte_carlo_reliability(G, p, 20000))
