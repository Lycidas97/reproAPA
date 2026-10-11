"""Physical-site changes; structure bounds are not sampling intervals."""
import hashlib
import numpy as np

STATES = ['S', 'R', 'E', 'P']
THRESHOLDS = [.05, .10, .15, .20]


def site_id(species, chrom, strand, start, end):
    return '|'.join(map(str, [species, chrom, strand, start, end]))


def factor_id(kind, names):
    return kind + ':' + hashlib.sha256(';'.join(sorted(names)).encode()).hexdigest()[:16]


def masks(start, end, strand, pas):
    if strand == '+':
        return np.array([[end <= p['member_cut_site_min'], end <= p['member_cut_site_max']] for p in pas], dtype=float)
    return np.array([[start >= p['member_cut_site_max']-1, start >= p['member_cut_site_min']-1] for p in pas], dtype=float)


def change(usage, reference, mask):
    contributions = (np.array(usage)-np.array(reference))[:, None] * mask
    return np.array([contributions.min(axis=1).sum(), contributions.max(axis=1).sum()])


def classify(lo, hi, threshold):
    if hi <= -threshold + 1e-12:
        return 'loss'
    if lo >= threshold - 1e-12:
        return 'gain'
    if lo > -threshold and hi < threshold:
        return 'below_threshold'
    return 'cut_uncertain'


def ambiguous_overlaps(intervals):
    """Return all coordinates in non-identical, overlapping components."""
    result = set()
    block = []
    edge = -1
    for span in sorted(set(intervals)):
        if span[0] >= edge:
            if len(block) > 1:
                result.update(block)
            block = []
        block.append(span)
        edge = max(edge, span[1]) if len(block) > 1 else span[1]
    if len(block) > 1:
        result.update(block)
    return result
