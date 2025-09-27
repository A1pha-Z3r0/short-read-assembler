# pylint: disable=trailing-whitespace
"""
This script has two helper functions to produce final contigs and
to write contigs to output fasta
"""

from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import networkx as nx
from utils.dbg import DBG

def merge_by_unique_k_overlap(seq_list : list, k : int):
    """
    Merge sequences by exact k-length suffix/prefix overlap, but ONLY when the
    next choice is unique. Ambiguities stop the chain. Returns merged contigs.
    """

    n = len(seq_list)
    if n == 0:
        return []

    # index k-prefix -> list of seq indices that start with that prefix
    pref_index = defaultdict(list)
    for i, s in enumerate(seq_list):
        if len(s) >= k:
            pref_index[s[:k]].append(i)

    succ = {}                 # i -> j (unique successor)
    indeg = [0] * n
    ambiguous = set()         # i where multiple successors exist (or conflicting succ set)

    # decide successors based on unique match of suffix(seq[i], k) -> prefix(seq[j], k)
    for i, s in enumerate(seq_list):
        if len(s) < k:
            continue
        cand = [j for j in pref_index.get(s[-k:], []) if j != i]
        if len(cand) == 1:
            j = cand[0]
            if i in succ:
                # conflicting successor discovered → mark ambiguous, retract
                ambiguous.add(i)
                del succ[i]
            else:
                succ[i] = j
                indeg[j] += 1
        elif len(cand) > 1:
            ambiguous.add(i)

    used = [False] * n
    contigs = []

    def stitch(chain):
        out = seq_list[chain[0]]
        for idx in chain[1:]:
            out += seq_list[idx][k:]
        return out

    # start chains at heads: indeg==0, not ambiguous, not used
    for i in range(n):
        if len(seq_list[i]) < k:
            continue
        if indeg[i] == 0 and i not in ambiguous and not used[i]:
            chain = [i]; used[i] = True
            cur = i
            while cur in succ and succ[cur] not in ambiguous and not used[succ[cur]]:
                cur = succ[cur]
                used[cur] = True
                chain.append(cur)
            contigs.append(stitch(chain))

    # anything not used becomes its own contig
    for i in range(n):
        if not used[i] and seq_list[i]:
            contigs.append(seq_list[i])

    contigs = [c for c in contigs if c]
    contigs.sort(key=len, reverse=True)
    return contigs



def kmer_counts_from_sequence(t_counts : dict,
                              seq: str,
                              k: int):
    """Mutates a dict in place to get counts of each k-mer

    Args:
        t_counts (dict): a dict containing key = kmer; value = count
        seq (str): the sequence read that is being processed
        k (int): the length of k-mer we want 

    Returns:
        t_counts (dict): a dict containing key = kmer; value = count
    """
    n = len(seq)
    for i in range(0, n - k + 1):
        kmer = seq[i:i+k]
        t_counts[kmer] = t_counts.get(kmer, 0) + 1
    return t_counts

def read_assembly(inp_file_path: str,
                out_file_path: str,
                k: int):
    records = []
    total_counts = {}
    
    # Detect FASTA by attempting to parse a first record
    try:
        first = next(SeqIO.parse(inp_file_path, "fasta-blast"), None)
        is_fasta = first is not None
    except Exception:
        is_fasta = False


    if is_fasta:
        for rec in SeqIO.parse(inp_file_path, "fasta-blast"):
            kmer_counts_from_sequence(t_counts = total_counts, seq = rec.seq, k=k)

    else:
        with open(inp_file_path, "r" , encoding="UTF-8") as file:
            for _, line in enumerate(file):
                seq = line.strip()
                if not seq:
                    continue
                kmer_counts_from_sequence(t_counts = total_counts, seq = seq, k=k)


        dbg = DBG(k=k)
        dbg.build_graph_weighted(total_counts)
        #dbg.print_g()

        wccs = dbg.wccs()
        for i, h_subgraph in enumerate(wccs, start=1):
            print(f"WCC {i}: |V|={h_subgraph.number_of_nodes()}, |E|={h_subgraph.number_of_edges()}")

        contigs = []
        for i, nodes in enumerate(nx.weakly_connected_components(dbg.G),
                         start=1):
            h_subgraph = dbg.G.subgraph(nodes).copy()

            seq = dbg.euler_contig_or_none(h_subgraph, label=f"WCC {i}")  # Eulerian? emit full contig
            if seq is not None:
                contigs.append(seq)
            else:
                uts = dbg.unitigs_weighted(h_subgraph)                     # not Eulerian → unitigs
                #print(f"WCC {i}: {len(uts)} unitig(s) before merge")
                merged = merge_by_unique_k_overlap(uts, k-1)        # <-- merge unitigs by unique k-overlap
                #print(f"WCC {i}: {len(merged)} contig(s) after merge")
                contigs.extend(merged)

        contigs.sort(key=len, reverse=True)
        #print(contigs)

        if not contigs:
            raise ValueError("Assembly produced no contigs")

        contigs.sort(key=len, reverse=True)

        records = [
            SeqRecord(Seq(c), id=f"Contig{i+1}", description=f"len={len(c)};k={k}")
            for i, c in enumerate(contigs)
        ]
        SeqIO.write(records, out_file_path, "fasta")

    return contigs
