"""Bounded, endpoint-locked display smoothing in final-figure millimetres."""
import numpy as np
from matplotlib.path import Path


def point_to_polyline(points, polyline):
    if len(points) > 512:
        return np.concatenate([point_to_polyline(part, polyline)
                               for part in np.array_split(points, int(np.ceil(len(points)/512)))])
    start, delta = polyline[:-1], np.diff(polyline, axis=0)
    denominator = np.sum(delta*delta, axis=1)
    keep = denominator > 1e-20
    start, delta, denominator = start[keep], delta[keep], denominator[keep]
    if not len(start):
        return np.linalg.norm(points-polyline[0], axis=1)
    rel = points[:, None, :] - start
    fraction = np.clip(np.sum(rel*delta, axis=2)/denominator, 0, 1)
    return np.linalg.norm(rel-fraction[:, :, None]*delta, axis=2).min(axis=1)


def simplify(points, tolerance):
    if len(points) <= 2:
        return points.copy()
    distance = point_to_polyline(points, points[[0, -1]])
    index = int(distance.argmax())
    if distance[index] <= tolerance:
        return points[[0, -1]].copy()
    return np.vstack([simplify(points[:index+1], tolerance)[:-1],
                      simplify(points[index:], tolerance)])


def rounded_path(points, cut_mm, step_mm=.04):
    vertices, codes, sampled = [points[0]], [Path.MOVETO], [points[:1]]
    current = points[0]

    def line_to(end):
        nonlocal current
        n = max(2, int(np.ceil(np.linalg.norm(end-current)/step_mm))+1)
        sampled.append(np.linspace(current, end, n)[1:])
        vertices.append(end)
        codes.append(Path.LINETO)
        current = end

    for before, corner, after in zip(points[:-2], points[1:-1], points[2:]):
        incoming, outgoing = corner-before, after-corner
        a, b = np.linalg.norm(incoming), np.linalg.norm(outgoing)
        cut = min(cut_mm, .45*a, .45*b)
        enter, leave = corner-cut*incoming/a, corner+cut*outgoing/b
        line_to(enter)
        n = max(3, int(np.ceil(2*cut/step_mm))+1)
        t = np.linspace(0, 1, n)[:, None]
        sampled.append(((1-t)**2*enter + 2*(1-t)*t*corner + t**2*leave)[1:])
        vertices.extend([corner, leave])
        codes.extend([Path.CURVE3, Path.CURVE3])
        current = leave
    line_to(points[-1])
    return Path(np.asarray(vertices), codes), np.vstack(sampled)


def smooth_chain(points_mm, tolerance_mm=.55, cut_mm=2.5, max_deviation_mm=1.5):
    """Remove sub-mm zigzags and round corners; never move either chain end.

    Caller splits at every topological branch, so locking chain ends also locks
    all roots, forks and terminals. Original timed nodes remain in the OT table;
    the resulting curve is explicitly a display-only geometry.
    """
    points = np.array(points_mm, dtype=float, copy=True)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 2 or not np.isfinite(points).all():
        raise ValueError('A smoothing chain needs at least two finite 2-D points')
    points = points[np.r_[True, np.linalg.norm(np.diff(points, axis=0), axis=1) > 1e-10]]
    if len(points) < 2:
        raise ValueError('A smoothing chain has zero length')
    reduced = simplify(points, tolerance_mm)
    # A dense sampling plus a 0.04-mm allowance bounds curve approximation and
    # distance sampling error; reduce rounding automatically if the cap fails.
    _, original_sample = rounded_path(points, 0.)
    for _ in range(12):
        path, sample = rounded_path(reduced, cut_mm)
        deviation = max(float(point_to_polyline(sample, points).max()),
                        float(point_to_polyline(original_sample, sample).max()))
        if deviation + .04 <= max_deviation_mm:
            break
        cut_mm *= .5
    else:
        raise ValueError('Smoothing exceeds the permitted geometric deviation')
    if not np.array_equal(path.vertices[[0, -1]], points[[0, -1]]):
        raise ValueError('Smoothing moved a root, fork or terminal')
    return path, {
        'original_vertices': len(points), 'retained_waypoints': len(reduced),
        'quadratic_corners': max(0, len(reduced)-2),
        'simplification_tolerance_mm': tolerance_mm, 'rounding_cut_mm': cut_mm,
        'max_sampled_deviation_mm': deviation, 'deviation_upper_bound_mm': deviation+.04,
        'locked_endpoint_shift_mm': 0.,
    }
