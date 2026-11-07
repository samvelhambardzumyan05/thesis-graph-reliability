import os, sys, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import networkx as nx
import connectivity as C
import reliability as R


def graphs():
    gs = [nx.path_graph(6), nx.cycle_graph(7), nx.complete_graph(5), nx.petersen_graph(),
          nx.grid_2d_graph(3, 3), nx.barbell_graph(4, 2), nx.star_graph(5), nx.wheel_graph(6)]
    for seed in range(15):
        gs.append(nx.gnp_random_graph(9, 0.3, seed=seed))
    return gs


def test_components():
    for G in graphs():
        adj = C.from_networkx(G)
        exp = sorted(map(sorted, nx.connected_components(G)), key=lambda c: c[0])
        for f in (C.components_bfs, C.components_dfs):
            got = sorted(map(sorted, f(adj)), key=lambda c: c[0])
            assert got == exp


def test_bridges_articulation():
    for G in graphs():
        adj = C.from_networkx(G)
        br, cut = C.bridges_and_articulation_points(adj)
        assert br == {frozenset(e) for e in nx.bridges(G)}
        assert cut == set(nx.articulation_points(G))


def test_connectivity_numbers():
    for G in graphs():
        if not nx.is_connected(G):
            continue
        adj = C.from_networkx(G)
        assert C.edge_connectivity(adj) == nx.edge_connectivity(G)
        assert C.vertex_connectivity(adj) == nx.node_connectivity(G)


def test_whitney_inequality():
    for G in graphs():
        if not nx.is_connected(G):
            continue
        adj = C.from_networkx(G)
        k, l = C.vertex_connectivity(adj), C.edge_connectivity(adj)
        assert k <= l <= min(dict(G.degree()).values())


def test_reliability_closed_forms():
    for p in (0.3, 0.7, 0.95):
        for n in (3, 5, 7):
            assert abs(R.exact_enumeration(nx.path_graph(n), p) - R.reliability_tree(n, p)) < 1e-12
            assert abs(R.exact_enumeration(nx.cycle_graph(n), p) - R.reliability_cycle(n, p)) < 1e-12
            assert abs(R.exact_factoring(nx.cycle_graph(n), p) - R.reliability_cycle(n, p)) < 1e-12


def test_factoring_equals_enumeration():
    for G in [nx.complete_graph(4), nx.wheel_graph(5), nx.grid_2d_graph(2, 4), nx.petersen_graph().subgraph(range(7)).copy()]:
        if not nx.is_connected(G) or G.number_of_edges() > 14:
            continue
        for p in (0.4, 0.9):
            assert abs(R.exact_enumeration(G, p) - R.exact_factoring(G, p)) < 1e-12


def test_monte_carlo_within_ci():
    G = nx.wheel_graph(6)
    p = 0.8
    ex = R.exact_factoring(G, p)
    est, hw = R.monte_carlo(G, p, 40000, seed=1)
    assert abs(est - ex) < 4 * hw


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
