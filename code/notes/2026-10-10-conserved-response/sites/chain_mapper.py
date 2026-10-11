"""Stream existing UCSC chains without files or a best-hit choice."""
from bisect import bisect_left
from collections import defaultdict
import gzip


def mapped_block_query(query, target_start, query_start, query_size, query_strand):
    lo, hi = query['start0'] - target_start, query['end0'] - target_start
    if query_strand == '+':
        return query_start + lo, query_start + hi, query['strand']
    return query_size - query_start - hi, query_size - query_start - lo, '-' if query['strand'] == '+' else '+'


def map_chain(lines, queries):
    """Each entire source interval must lie in one gap-free alignment block."""
    by_chrom = defaultdict(list)
    result = {q['query_id']: dict(hits=[], partial_overlap=False) for q in queries}
    for query in queries:
        assert query['start0'] < query['end0']
        by_chrom[query['chrom']].append(query)
    by_chrom = {chrom: sorted(values, key=lambda q: q['start0']) for chrom, values in by_chrom.items()}
    starts = {chrom: [q['start0'] for q in values] for chrom, values in by_chrom.items()}
    maximum_length = max([q['end0'] - q['start0'] for q in queries] + [1])
    header, pending = None, None

    def emit(block):
        if block is None or header['tName'] not in by_chrom:
            return
        tpos, qpos, size = block
        chrom = header['tName']
        lo = bisect_left(starts[chrom], tpos - maximum_length + 1)
        hi = bisect_left(starts[chrom], tpos + size)
        for query in by_chrom[chrom][lo:hi]:
            if query['end0'] <= tpos:
                continue
            item = result[query['query_id']]
            item['partial_overlap'] = True
            if query['start0'] >= tpos and query['end0'] <= tpos + size:
                a, b, strand = mapped_block_query(query, tpos, qpos, header['qSize'], header['qStrand'])
                assert b - a == query['end0'] - query['start0']
                item['hits'].append(dict(chrom=header['qName'], strand=strand, start0=a, end0=b, chain_id=header['id']))

    for line in lines:
        fields = line.split()
        if not fields:
            emit(pending); pending = None
            continue
        if fields[0] == 'chain':
            emit(pending); pending = None
            assert len(fields) == 13 and fields[4] == '+'
            header = dict(tName=fields[2], qName=fields[7], qSize=int(fields[8]), qStrand=fields[9], id=fields[12])
            tpos, qpos = int(fields[5]), int(fields[10])
            continue
        size = int(fields[0])
        if pending is not None and pending[0] + pending[2] == tpos and pending[1] + pending[2] == qpos:
            pending = (pending[0], pending[1], pending[2] + size)
        else:
            emit(pending); pending = (tpos, qpos, size)
        tpos += size; qpos += size
        if len(fields) == 3:
            dt, dq = int(fields[1]), int(fields[2])
            if dt or dq:
                emit(pending); pending = None
            tpos += dt; qpos += dq
        else:
            assert len(fields) == 1
            emit(pending); pending = None
    emit(pending)
    for item in result.values():
        item['status'] = 'unique' if len(item['hits']) == 1 else 'ambiguous' if len(item['hits']) > 1 else 'gap_or_partial' if item['partial_overlap'] else 'unmapped'
    return result


if __name__ == '__main__':
    import hashlib
    import json
    import sys
    from pathlib import Path
    for species, key in [('human', 'chain_hg38_mm39'), ('mouse', 'chain_mm39_hg38')]:
        resource = REMOTE_RESOURCES[key]
        path = Path(resource['path'])
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1048576), b''):
                digest.update(block)
        assert digest.hexdigest() == resource['sha256']
        print(json.dumps(dict(resource=key, path=str(path), sha256=digest.hexdigest(), bytes=path.stat().st_size)), file=sys.stderr, flush=True)
        queries = [q for q in REMOTE_QUERIES if q['species'] == species]
        with gzip.open(path, 'rt') as stream:
            output = map_chain(stream, queries)
        for query in queries:
            print(json.dumps(dict(query_id=query['query_id'], species=species, **output[query['query_id']])), flush=True)
