"""Local sequence alignment, used to map an epitope onto a paralog.

Note on the substitution matrix: this uses a similarity matrix derived from
the physicochemical property tables in aminoacids.py (hydrophobicity, volume,
charge), not BLOSUM62. That is a deliberate choice -- the properties are the
thing we actually care about when asking "would a binder shaped for surface A
also fit surface B", and deriving the matrix from data already in the repo
keeps it auditable. It is not a drop-in replacement for BLOSUM62 for
phylogenetic work.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .aminoacids import AA20, EISENBERG, FORMAL_CHARGE, RESIDUE_VOLUME

GAP_OPEN = -11.0
GAP_EXTEND = -1.0


@lru_cache(maxsize=1)
def property_matrix() -> dict[tuple[str, str], float]:
    """Similarity in 'BLOSUM-like' units: ~+6 identical, ~-4 very dissimilar."""
    # Normalising constants so each property contributes comparably.
    h_span = max(EISENBERG.values()) - min(EISENBERG.values())
    v_span = max(RESIDUE_VOLUME.values()) - min(RESIDUE_VOLUME.values())

    mat: dict[tuple[str, str], float] = {}
    for a in AA20:
        for b in AA20:
            if a == b:
                mat[(a, b)] = 6.0
                continue
            dh = abs(EISENBERG[a] - EISENBERG[b]) / h_span
            dv = abs(RESIDUE_VOLUME[a] - RESIDUE_VOLUME[b]) / v_span
            dq = abs(FORMAL_CHARGE[a] - FORMAL_CHARGE[b]) / 2.0
            dissimilarity = 0.45 * dh + 0.25 * dv + 0.30 * dq
            mat[(a, b)] = 4.0 - 10.0 * dissimilarity
    return mat


def sub_score(a: str, b: str) -> float:
    if a not in AA20 or b not in AA20:
        return -4.0
    return property_matrix()[(a, b)]


@dataclass
class Alignment:
    score: float
    query_aligned: str
    target_aligned: str
    query_start: int
    target_start: int
    # query index -> target index, for positions aligned without a gap
    mapping: dict[int, int]

    @property
    def identity(self) -> float:
        pairs = [(q, t) for q, t in zip(self.query_aligned, self.target_aligned)
                 if q != "-" and t != "-"]
        if not pairs:
            return 0.0
        return sum(1 for q, t in pairs if q == t) / len(pairs)

    @property
    def length(self) -> int:
        return len(self.query_aligned)


def smith_waterman(query: str, target: str,
                   gap_open: float = GAP_OPEN,
                   gap_extend: float = GAP_EXTEND) -> Alignment:
    """Local alignment with affine gaps.

    Plain dynamic programming. The sequences here are a few hundred residues
    at most, so the O(nm) cost is irrelevant and clarity wins.
    """
    n, m = len(query), len(target)
    if n == 0 or m == 0:
        return Alignment(0.0, "", "", 0, 0, {})

    NEG = float("-inf")
    # M: ends in a match. X: gap in target. Y: gap in query.
    M = [[0.0] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    ptr = [[0] * (m + 1) for _ in range(n + 1)]  # 0 stop, 1 M, 2 X, 3 Y

    best, bi, bj = 0.0, 0, 0
    for i in range(1, n + 1):
        qi = query[i - 1]
        for j in range(1, m + 1):
            s = sub_score(qi, target[j - 1])

            prev = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1], 0.0)
            M[i][j] = prev + s
            if M[i][j] < 0:
                M[i][j] = 0.0

            X[i][j] = max(M[i - 1][j] + gap_open, X[i - 1][j] + gap_extend)
            Y[i][j] = max(M[i][j - 1] + gap_open, Y[i][j - 1] + gap_extend)

            cell = max(M[i][j], X[i][j], Y[i][j])
            if cell <= 0:
                ptr[i][j] = 0
            elif cell == M[i][j]:
                ptr[i][j] = 1
            elif cell == X[i][j]:
                ptr[i][j] = 2
            else:
                ptr[i][j] = 3

            if cell > best:
                best, bi, bj = cell, i, j

    # Traceback
    qa: list[str] = []
    ta: list[str] = []
    mapping: dict[int, int] = {}
    i, j = bi, bj
    while i > 0 and j > 0 and ptr[i][j] != 0:
        move = ptr[i][j]
        if move == 1:
            qa.append(query[i - 1])
            ta.append(target[j - 1])
            mapping[i - 1] = j - 1
            i, j = i - 1, j - 1
        elif move == 2:
            qa.append(query[i - 1])
            ta.append("-")
            i -= 1
        else:
            qa.append("-")
            ta.append(target[j - 1])
            j -= 1

    qa.reverse()
    ta.reverse()
    return Alignment(
        score=best,
        query_aligned="".join(qa),
        target_aligned="".join(ta),
        query_start=i,
        target_start=j,
        mapping=mapping,
    )


def map_positions(query: str, target: str,
                  positions: list[int]) -> dict[int, str | None]:
    """For each query index, the aligned residue in target (None if gapped)."""
    aln = smith_waterman(query, target)
    out: dict[int, str | None] = {}
    for p in positions:
        tj = aln.mapping.get(p)
        out[p] = target[tj] if tj is not None else None
    return out
