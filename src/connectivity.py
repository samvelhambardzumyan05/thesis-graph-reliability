"""Connectivity algorithms implemented from scratch (no networkx inside the algorithms).

Graphs are given as adjacency dicts: {vertex: set(neighbours)} (simple undirected graphs).
Helper `from_networkx` converts an nx.Graph.
"""
from collections import deque


def from_networkx(G):
    return {v: set(G.neighbors(v)) for v in G.nodes()}


# ---------------------------------------------------------------- components
def components_bfs(adj):
    """Connected components via BFS. O(n + m)."""
    seen, comps = set(), []
    for s in adj:
        if s in seen:
            continue
        comp, queue = {s}, deque([s])
        seen.add(s)
        while queue:
            v = queue.popleft()
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    comp.add(w)
                    queue.append(w)
        comps.append(comp)
    return comps


def components_dfs(adj):
    """Connected components via iterative DFS. O(n + m)."""
    seen, comps = set(), []
    for s in adj:
        if s in seen:
            continue
        comp, stack = set(), [s]
        while stack:
            v = stack.pop()
            if v in seen:
                continue
            seen.add(v)
            comp.add(v)
            stack.extend(w for w in adj[v] if w not in seen)
        comps.append(comp)
    return comps


def is_connected(adj):
    return len(components_bfs(adj)) <= 1


# ------------------------------------------------- bridges / articulation points
def bridges_and_articulation_points(adj):
    """Tarjan's low-link algorithm (recursive DFS). O(n + m).

    disc[v]  - discovery time of v
    low[v]   - smallest discovery time reachable from the DFS subtree of v
               using at most one back edge
    Edge (parent, v) is a bridge  iff low[v] > disc[parent].
    Vertex u is an articulation point iff
       (u is a DFS root with >= 2 children) or (u is not root and low[child] >= disc[u]).
    """
    import sys

    sys.setrecursionlimit(max(10000, 4 * len(adj)))
    disc, low = {}, {}
    bridges, cut = set(), set()
    timer = [0]

    def dfs(u, parent):
        disc[u] = low[u] = timer[0]
        timer[0] += 1
        children = 0
        for v in adj[u]:
            if v == parent:
                continue
            if v in disc:
                low[u] = min(low[u], disc[v])
            else:
                children += 1
                dfs(v, u)
                low[u] = min(low[u], low[v])
                if low[v] > disc[u]:
                    bridges.add(frozenset((u, v)))
                if parent is not None and low[v] >= disc[u]:
                    cut.add(u)
        if parent is None and children >= 2:
            cut.add(u)

    for s in adj:
        if s not in disc:
            dfs(s, None)
    return bridges, cut


# --------------------------------------------------------------- max flow / cuts
def _max_flow_unit(cap, s, t):
    """Edmonds-Karp on a dict-of-dict capacity graph (modified in place). Returns flow value."""
    flow = 0
    while True:
        parent = {s: None}
        queue = deque([s])
        while queue and t not in parent:
            v = queue.popleft()
            for w, c in cap[v].items():
                if c > 0 and w not in parent:
                    parent[w] = v
                    queue.append(w)
        if t not in parent:
            return flow
        # all capacities are integers; bottleneck of the augmenting path
        b, w = float("inf"), t
        while parent[w] is not None:
            b = min(b, cap[parent[w]][w])
            w = parent[w]
        w = t
        while parent[w] is not None:
            cap[parent[w]][w] -= b
            cap[w][parent[w]] += b
            w = parent[w]
        flow += b


def max_edge_disjoint_paths(adj, s, t):
    """Max number of edge-disjoint s-t paths = min s-t edge cut (Menger, edge version)."""
    cap = {v: {w: 1 for w in adj[v]} for v in adj}
    return _max_flow_unit(cap, s, t)


def max_vertex_disjoint_paths(adj, s, t):
    """Max number of internally vertex-disjoint s-t paths (Menger, vertex version).

    Every vertex v (other than s,t) is split into v_in -> v_out with capacity 1.
    """
    INF = len(adj) + 5
    cap = {}

    def add(u, v, c):
        cap.setdefault(u, {})[v] = cap.get(u, {}).get(v, 0) + c
        cap.setdefault(v, {}).setdefault(u, 0)

    for v in adj:
        add((v, "in"), (v, "out"), INF if v in (s, t) else 1)
    for v in adj:
        for w in adj[v]:
            add((v, "out"), (w, "in"), INF)
    return _max_flow_unit(cap, (s, "out"), (t, "in"))


def edge_connectivity(adj):
    """lambda(G): fix one vertex s, minimise max-flow(s, t) over all t != s."""
    vs = list(adj)
    if len(vs) < 2:
        return 0
    s = vs[0]
    return min(max_edge_disjoint_paths(adj, s, t) for t in vs[1:])


def vertex_connectivity(adj):
    """kappa(G) by checking all non-adjacent pairs (small graphs). Complete graph: n-1."""
    vs = list(adj)
    n = len(vs)
    best = n - 1
    for i in range(n):
        for j in range(i + 1, n):
            u, v = vs[i], vs[j]
            if v not in adj[u]:
                best = min(best, max_vertex_disjoint_paths(adj, u, v))
    return best
