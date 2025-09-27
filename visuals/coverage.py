import os
import matplotlib.pyplot as plt


def coverage_vector_unique(contig: str, reads: list[str], k: int) -> list[int]:
    """
    Per-base coverage using unique placements only.
    1) Exact match(s): if exactly one occurrence in contig, count it; if 0 or >1, don't.
    2) Else, k-mer–anchored placement:
       - For each k-mer in the read, find where it occurs in the contig.
       - Each hit implies an alignment start = contig_pos - kmer_index_in_read.
       - If exactly one distinct start is implied across all anchors, count it.
       - Otherwise (0 or >1 starts), skip the read.
    Returns coverage vector (len(contig)).
    """
    n = len(contig)
    cov = [0] * n
    if n == 0:
        return cov

    for r in reads:
        if not r:
            continue
        # 1) Exact placement(s)
        starts = []
        scan = 0
        while True:
            pos = contig.find(r, scan)
            if pos == -1:
                break
            starts.append(pos)
            scan = pos + 1  # allow overlapping exact hits
        if len(starts) == 1:
            s = starts[0]
            e = min(s + len(r), n)
            for i in range(s, e):
                cov[i] += 1
            continue
        if len(starts) > 1:
            # ambiguous exact matches → skip (conservative)
            continue

        # 2) k-mer–anchored placement
        kk = min(k, len(r))
        candidate_starts = set()
        for i in range(0, len(r) - kk + 1):
            anchor = r[i:i+kk]
            seek = 0
            while True:
                j = contig.find(anchor, seek)
                if j == -1:
                    break
                start_pos = j - i  # where r[0] would align on contig
                # Keep candidates that overlap the contig at least 1 base
                if start_pos < n and start_pos + len(r) > 0:
                    candidate_starts.add(start_pos)
                seek = j + 1

        if len(candidate_starts) == 1:
            s = next(iter(candidate_starts))
            # clamp to contig bounds (count only the overlapping slice)
            left = max(s, 0)
            right = min(s + len(r), n)
            if left < right:
                for i in range(left, right):
                    cov[i] += 1
        # else: 0 or >1 starts → skip (conservative)

    return cov


def save_coverage_hist(contig_id: str, coverage: list[int], bins: int = 30,
                       out_dir: str = "coverage_plots"):
    os.makedirs(out_dir, exist_ok=True)
    if not coverage:
        return None
    plt.figure(figsize=(4.5, 3))
    plt.hist(coverage, bins=bins)
    plt.xlabel("Per-base coverage")
    plt.ylabel("# positions")
    plt.title(f"{contig_id} coverage histogram")
    plt.tight_layout()
    out_path = os.path.join(out_dir, f"{contig_id}_coverage_hist.png")
    plt.savefig(out_path, dpi=200)
    plt.close()
    return out_path