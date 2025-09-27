import networkx as nx
import matplotlib.pyplot as plt

class DBG:
    def __init__(self, k: int):
        self.k = k
        self.G = nx.DiGraph()
        self.contigs = []

    def build_graph_weighted(self, kmer_counts: dict):
        for kmer, c in kmer_counts.items():
            u = kmer[:-1]; v = kmer[1:]
            prev = self.G.get_edge_data(u, v, default={'w': 0})['w']
            self.G.add_edge(u, v, w=prev + c)
        return self.G
    
    def print_g(self):
        edges_with_data = list(self.G.edges(data=True))
        print("Edges with data:", edges_with_data)

    def wccs(self):
        """Return a list of subgraphs, one per weakly connected component (island)."""
        return [self.G.subgraph(nodes).copy()
                for nodes in nx.weakly_connected_components(self.G)]

    def euler_contig_or_none(self, H, label="WCC"):
        """
        If H (a subgraph/WCC) is Eulerian (weighted), emit its contig and print a note.
        Otherwise print that it's not Eulerian and return None.
        """
        # nodes with any incident weight
        nodes = [n for n in H.nodes if H.in_degree(n, weight='w') + H.out_degree(n, weight='w') > 0]
        if not nodes or not nx.is_weakly_connected(H.subgraph(nodes)):
            print(f"{label} is NOT Eulerian.")
            return None

        # weighted degree balances
        diffs = {n: H.out_degree(n, weight='w') - H.in_degree(n, weight='w') for n in nodes}
        plus  = [n for n, d in diffs.items() if d == 1]
        minus = [n for n, d in diffs.items() if d == -1]
        balanced = all(d == 0 for d in diffs.values())

        if balanced:
            status = "circuit"
            start = next((n for n in nodes if H.out_degree(n, weight='w') > 0), None)
        elif len(plus) == 1 and len(minus) == 1:
            status = "path"
            start = plus[0]
        else:
            print(f"{label} is NOT Eulerian.")
            return None

        # Hierholzer-style traversal consuming 'w'
        H = H.copy()
        edges = []
        stack = [start]
        while stack:
            v = stack[-1]
            advanced = False
            for w in sorted(H.successors(v)):  # sorted = stable/deterministic
                if H.has_edge(v, w) and H[v][w].get('w', 0) > 0:
                    H[v][w]['w'] -= 1
                    if H[v][w]['w'] == 0:
                        H.remove_edge(v, w)
                    stack.append(w)
                    edges.append((v, w))
                    advanced = True
                    break
            if not advanced:
                stack.pop()

        if not edges:
            print(f"{label} is NOT Eulerian.")
            return None

        seq = self._contig_from_edge_list(edges)
        print(f"{label} is Eulerian ({status}). Emitting contig of len={len(seq)}.")
        return seq
    
    def unitigs_weighted(self, H=None):
        """
        Return unitigs (maximal non-branching paths) in H using weighted degrees.
        Break whenever in_w != 1 or out_w != 1 (i.e., at repeats/ambiguity).
        Does not consume/remove edges; uses a seen set to avoid re-emitting edges.
        """
        H = self.G if H is None else H
        seen = set()     # track directed edges (u,v) we already emitted
        outs = []

        def in_w(n):  return H.in_degree(n, weight='w')
        def out_w(n): return H.out_degree(n, weight='w')
        def is_branch(n):  # branching if not exactly one in and one out (by weight)
            return in_w(n) != 1 or out_w(n) != 1

        # Start unitigs from every outgoing edge of a branch node
        for u in H.nodes:
            if out_w(u) == 0:
                continue
            if is_branch(u):
                for _, v in H.out_edges(u):
                    if (u, v) in seen or H[u][v].get('w', 0) <= 0:
                        continue
                    path = [(u, v)]
                    seen.add((u, v))
                    cur = v
                    # extend through non-branching nodes with a unique outgoing edge
                    while not is_branch(cur) and out_w(cur) == 1:
                        nxt = next(iter(H.successors(cur)))
                        if (cur, nxt) in seen or H[cur][nxt].get('w', 0) <= 0:
                            break
                        seen.add((cur, nxt))
                        path.append((cur, nxt))
                        cur = nxt
                    seq = self._contig_from_edge_list(path)
                    if seq:
                        outs.append(seq)

        # Handle pure simple cycles: every node on the cycle has in_w = out_w = 1
        for u, v, d in H.edges(data=True):
            if (u, v) in seen or d.get('w', 0) <= 0:
                continue
            if in_w(u) == 1 and out_w(u) == 1:
                path = [(u, v)]
                seen.add((u, v))
                cur = v
                ok = True
                while cur != u:
                    if in_w(cur) != 1 or out_w(cur) != 1:
                        ok = False
                        break
                    nxt = next(iter(H.successors(cur)))
                    if (cur, nxt) in seen or H[cur][nxt].get('w', 0) <= 0:
                        break
                    seen.add((cur, nxt))
                    path.append((cur, nxt))
                    cur = nxt
                if ok:
                    seq = self._contig_from_edge_list(path)
                    if seq:
                        outs.append(seq)

        outs.sort(key=len, reverse=True)
        return outs

    def _contig_from_edge_list(self, edge_path):
        if not edge_path:
            return ""
        # Normalize tuples to (u, v)
        first = edge_path[0]
        if len(first) == 3:
            u0, _, _ = first
        else:
            u0, _ = first
        seq = u0
        for e in edge_path:
            if len(e) == 3:
                _, v, _ = e
            else:
                _, v = e
            seq += v[-1]
        return seq


    def show(self):
        pos = nx.spring_layout(self.G, seed=0)
        nx.draw(self.G, pos, with_labels=True, node_size=800, font_size=8)
        plt.show()
    