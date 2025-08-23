import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from dataclasses import dataclass

# --- Core Self-Implemented Geometric Primitives ---

def is_inside_polygon(point, polygon_vertices):
    x, y = point
    n = len(polygon_vertices)
    inside = False
    p1x, p1y = polygon_vertices[0]
    for i in range(n + 1):
        p2x, p2y = polygon_vertices[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside

def polygon_area(vertices):
    if vertices.shape[0] < 3: return 0.0
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

def polygon_centroid(vertices):
    if vertices.shape[0] < 3: return np.mean(vertices, axis=0)
    x, y = vertices[:, 0], vertices[:, 1]
    cross_prod_term = x * np.roll(y, -1) - np.roll(x, -1) * y
    area_val = 0.5 * np.sum(cross_prod_term)
    if np.abs(area_val) < 1e-9: return np.mean(vertices, axis=0)
    cx = (1 / (6 * area_val)) * np.sum((x + np.roll(x, -1)) * cross_prod_term)
    cy = (1 / (6 * area_val)) * np.sum((y + np.roll(y, -1)) * cross_prod_term)
    return np.array([cx, cy])

def sutherland_hodgman_clip(subject_polygon, clip_polygon):
    def is_inside_edge(p, p1, p2):
        return (p2[0] - p1[0]) * (p[1] - p1[1]) - (p2[1] - p1[1]) * (p[0] - p1[0]) >= 0

    def intersection(s1, s2, p1, p2):
        # *** FIX 2: Replace np.cross to remove DeprecationWarning ***
        dc = p1 - p2
        dp = s1 - s2
        n1 = s1[0]*s2[1] - s1[1]*s2[0] # cross(s1, s2)
        n2 = p1[0]*p2[1] - p1[1]*p2[0] # cross(p1, p2)
        denominator = dp[0]*dc[1] - dp[1]*dc[0] # cross(dp, dc)
        if abs(denominator) < 1e-9: return None
        return (n1 * dc - n2 * dp) / denominator

    output_list = subject_polygon
    clip_len = len(clip_polygon)
    
    for i in range(clip_len):
        clip_p1 = clip_polygon[i]
        clip_p2 = clip_polygon[(i + 1) % clip_len]
        
        input_list = output_list
        output_list = []
        
        # This check works because input_list is always a NumPy array now
        if input_list.size == 0:
            return np.array([])

        s = input_list[-1]
        for j in range(len(input_list)):
            e = input_list[j]
            s_inside = is_inside_edge(s, clip_p1, clip_p2)
            e_inside = is_inside_edge(e, clip_p1, clip_p2)

            if e_inside:
                if not s_inside:
                    inter = intersection(s, e, clip_p1, clip_p2)
                    if inter is not None: output_list.append(inter)
                output_list.append(e)
            elif s_inside:
                inter = intersection(s, e, clip_p1, clip_p2)
                if inter is not None: output_list.append(inter)
            s = e
        
        # *** FIX 1: Convert list back to NumPy array for the next iteration ***
        output_list = np.array(output_list) if len(output_list) > 0 else np.array([])
            
    return output_list

# --- Hex Grid Logic ---

def axial_to_cartesian(q, r, size):
    x = size * (3.0/2.0 * q)
    y = size * (np.sqrt(3)/2.0 * q + np.sqrt(3) * r)
    return np.array([x, y])

def create_pointy_hexagon_vertices(center, size):
    # Flat-top orientation to match axial_to_cartesian mapping
    corners = [[center[0] + size * np.cos(i * np.pi/3), center[1] + size * np.sin(i * np.pi/3)] for i in range(6)]
    return np.array(corners)

NEIGHBOR_MAP = [(1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1), (1, 0)]

def fmt_qr(k):
    return f"({int(k[0])},{int(k[1])})"

# --- Main Rasterizer Class ---

@dataclass
class Element:
    geoms: list
    area: float
    com: np.ndarray
    is_base: bool
    id: str
    com_orig: np.ndarray | None = None
    merged_piece_info: list | None = None

@dataclass
class Frag:
    qr: tuple
    i: int
    geom: np.ndarray
    area: float
    com: np.ndarray

class HexRasterizerNP:
    def __init__(self, polygon_vertices, hex_size, verbosity=1, update_com=True):
        self.polygon_verts = np.array(polygon_vertices)
        self.hex_size = float(hex_size)
        self.full_hex_area = 1.5 * np.sqrt(3) * self.hex_size**2
        self.verbosity = verbosity
        self.update_com = bool(update_com)  # True: dynamic COM updates on merge; False: keep base COM static

    def rasterize(self, area_keep_ratio=0.7):
        if self.verbosity >= 1: print("Stage 1: Generating hexagonal grid...")
        hex_centers, hex_verts_map = self._generate_grid()

        # Two-phase pipeline
        if self.verbosity >= 1: print("Stage 2: First pass (collect bases and fragments)...")
        elements, fragments = self._first_pass_collect(hex_centers, hex_verts_map, area_keep_ratio)

        if self.verbosity >= 1: print("Stage 3: Second pass (merge fragments into nearest bases)...")
        self._merge_fragments(fragments, hex_centers, hex_verts_map, elements)
        
        if self.verbosity >= 1:
            original_area = polygon_area(self.polygon_verts)
            total_element_area = sum(el.area for el in elements.values())
            print("\nSanity Check:")
            print(f"Original Polygon Area: {original_area:.4f}")
            print(f"Total Element Area:    {total_element_area:.4f}")
            if original_area > 1e-9:
                 print(f"Difference (%):        {100 * abs(original_area - total_element_area) / original_area:.4f}%")
        
        return elements

    def _generate_grid(self):
        min_coords = np.min(self.polygon_verts, axis=0)
        max_coords = np.max(self.polygon_verts, axis=0)
        minx, miny, maxx, maxy = min_coords[0], min_coords[1], max_coords[0], max_coords[1]
        
        q_min = int(np.floor((2/3 * minx) / self.hex_size)) - 2
        q_max = int(np.ceil((2/3 * maxx) / self.hex_size)) + 2
        r_min = int(np.floor((-1/3 * maxx + np.sqrt(3)/3 * miny) / self.hex_size)) - 2
        r_max = int(np.ceil((-1/3 * minx + np.sqrt(3)/3 * maxy) / self.hex_size)) + 2

        q_grid, r_grid = np.meshgrid(np.arange(q_min, q_max), np.arange(r_min, r_max))
        hex_indices = np.stack([q_grid.ravel(), r_grid.ravel()], axis=-1)
        
        hex_centers, hex_verts_map = {}, {}
        for q, r in hex_indices:
            qr = (q, r)
            center = axial_to_cartesian(q, r, self.hex_size)
            if (center[0] + self.hex_size > minx and center[0] - self.hex_size < maxx and
                center[1] + self.hex_size > miny and center[1] - self.hex_size < maxy):
                hex_centers[qr] = center
                hex_verts_map[qr] = create_pointy_hexagon_vertices(center, self.hex_size)
        return hex_centers, hex_verts_map

    def _classify_hexagons(self, hex_verts_map):
        full_indices, partial_indices = [], []
        for qr, vertices in hex_verts_map.items():
            if all(is_inside_polygon(v, self.polygon_verts) for v in vertices):
                full_indices.append(qr)
            else:
                clipped = sutherland_hodgman_clip(self.polygon_verts, vertices)
                if clipped.size > 0 and polygon_area(clipped) > 1e-9:
                    partial_indices.append(qr)
        return full_indices, partial_indices
    
    def _create_initial_elements(self, full_indices, hex_verts_map, hex_centers):
        # Initialize base (full) hex elements; keep original COM for plotting/debug
        return {
            qr: Element(geoms=[hex_verts_map[qr]], area=self.full_hex_area, com=hex_centers[qr].copy(), is_base=True, id=f"{qr[0]},{qr[1]}", com_orig=hex_centers[qr].copy(), merged_piece_info=None)
            for qr in full_indices
        }

    # --- Two-phase: First pass ---
    def _first_pass_collect(self, hex_centers, hex_verts_map, area_keep_ratio):
        elements: dict[tuple, Element] = {}
        fragments: list[Frag] = []
        for qr, vertices in hex_verts_map.items():
            q, r = qr
            center = hex_centers[qr]
            inside_all = all(is_inside_polygon(v, self.polygon_verts) for v in vertices)
            if inside_all:
                # Full hex -> base element
                elements[qr] = Element(geoms=[vertices], area=self.full_hex_area, com=center.copy(), is_base=True, id=f"{q},{r}", com_orig=center.copy(), merged_piece_info=None)
                if self.verbosity >= 2:
                    print(f"BASE_FULL: base={q},{r} A_full={self.full_hex_area:.6f} COM={center}")
                continue

            # Clip hex by polygon
            clipped_hex = sutherland_hodgman_clip(self.polygon_verts, vertices)
            if clipped_hex.size > 0 and clipped_hex.shape[0] >= 3:
                a_clip = polygon_area(clipped_hex)
                ratio = a_clip / self.full_hex_area
            else:
                a_clip = 0.0
                ratio = 0.0

            if ratio >= area_keep_ratio and a_clip > 1e-12:
                # Large partial -> base element with corrected COM
                com_clip = polygon_centroid(clipped_hex)
                elements[qr] = Element(geoms=[clipped_hex], area=a_clip, com=com_clip.copy(), is_base=True, id=f"{q},{r}", com_orig=com_clip.copy(), merged_piece_info=None)
                if self.verbosity >= 2:
                    print(f"BASE_PART: base={q},{r} A_clip={a_clip:.6f} ratio={ratio:.3f} COM={com_clip}")
                continue

            # Otherwise decompose into 6 triangles and collect clipped fragments
            for i in range(6):
                p1, p2 = vertices[i], vertices[(i + 1) % 6]
                triangle_verts = np.array([center, p1, p2])
                clipped_piece = sutherland_hodgman_clip(self.polygon_verts, triangle_verts)
                if clipped_piece.shape[0] < 3:
                    if self.verbosity >= 4: print(f"OUT: frag={q},{r},{i} -> clipped 0")
                    continue
                piece_area = polygon_area(clipped_piece)
                if piece_area < 1e-12:
                    if self.verbosity >= 4: print(f"SMALL: frag={q},{r},{i} -> clipped ~0")
                    continue
                piece_com = polygon_centroid(clipped_piece)
                fragments.append(Frag(qr=qr, i=i, geom=clipped_piece, area=piece_area, com=piece_com))
        return elements, fragments

    # --- Two-phase: Second pass ---
    def _merge_fragments(self, fragments, hex_centers, hex_verts_map, elements):
        for frag in fragments:
            qr = frag.qr
            q, r = qr
            i = frag.i
            clipped_piece = frag.geom
            piece_area = frag.area
            piece_com = frag.com

            # Candidate bases
            def base_keys(): return [k for k, v in elements.items() if isinstance(k, tuple) and v.is_base]
            candidates = base_keys()
            if len(candidates) == 0:
                # Seed zero-area base to enable merging when absolutely no bases exist yet
                if (qr not in elements) or (not elements[qr].is_base):
                    elements[qr] = Element(geoms=[], area=0.0, com=hex_centers[qr].copy(), is_base=True, id=f"{q},{r}", com_orig=hex_centers[qr].copy(), merged_piece_info=None)
            candidates = base_keys()

            self._dbg_nei_and_candidates(q, r, i, candidates, piece_com, elements)

            # Adjacency-first selection: find neighbor across the shared edge (i,i+1)
            verts = hex_verts_map[qr]
            p1, p2 = verts[i], verts[(i + 1) % 6]
            edge_mid = 0.5 * (p1 + p2)
            nei = [(q + dq, r + dr) for (dq, dr) in NEIGHBOR_MAP]
            nei_centers = [axial_to_cartesian(nk[0], nk[1], self.hex_size) for nk in nei]
            d_mid = np.linalg.norm(np.array(nei_centers) - edge_mid, axis=1)
            adj_key = nei[int(np.argmin(d_mid))]
            if (adj_key in elements) and elements[adj_key].is_base:
                best_key = adj_key
            else:
                # restrict search to 1-ring neighbor bases if available
                nei_bases = [k for k in nei if (k in elements) and elements[k].is_base]
                if len(nei_bases) > 0:
                    nb_coms = np.array([elements[k].com for k in nei_bases])
                    nb_dists = np.linalg.norm(nb_coms - piece_com, axis=1)
                    best_key = nei_bases[int(np.argmin(nb_dists))]
                else:
                    # fall back to global candidates
                    base_coms = np.array([elements[k].com for k in candidates])
                    dists = np.linalg.norm(base_coms - piece_com, axis=1)
                    best_key = candidates[int(np.argmin(dists))]
            if self.verbosity > 2:
                print(f"ADJ:  frag={q},{r},{i} -> edge->{fmt_qr(adj_key)} {'OK' if (adj_key == best_key) else '-'}; sel={fmt_qr(best_key)} mode={'adj' if (best_key == adj_key) else ('nei' if (best_key in nei) else 'global')}")
            neighbor_el = elements[best_key]
            old_area = neighbor_el.area
            old_com = neighbor_el.com
            new_area = old_area + piece_area
            if self.update_com:
                neighbor_el.com = (old_com * old_area + piece_com * piece_area) / new_area
            neighbor_el.area = new_area
            neighbor_el.geoms.append(clipped_piece)
            frag_id = f"{q},{r},{i}"
            # Keep lightweight diagnostic trail only if verbosity requests it
            if self.verbosity > 1:
                if neighbor_el.merged_piece_info is None: neighbor_el.merged_piece_info = []
                neighbor_el.merged_piece_info.append({'id': frag_id, 'com': piece_com})
                print(f"MERG: frag={frag_id} A_clip={piece_area:.6f} -> base{fmt_qr(best_key)} A:{old_area:.6f}->{new_area:.6f} {'COM:'+str(old_com)+'->'+str(neighbor_el.com) if self.update_com else 'COM static '+str(old_com)}")

    def _dbg_nei_and_candidates(self, q, r, i, candidates, piece_com, elements):
        if self.verbosity <= 3: return
        base_coms = np.array([elements[k].com for k in candidates]) if len(candidates) else np.zeros((0,2))
        dists = np.linalg.norm(base_coms - piece_com, axis=1) if len(candidates) else np.array([])
        nei = [(q + dq, r + dr) for (dq, dr) in NEIGHBOR_MAP]
        nei_info = [f"{fmt_qr(nk)}:{'B' if ((nk in elements) and elements[nk].is_base) else '-'}" for nk in nei]
        print(f"NEI:  frag={q},{r},{i} -> 1-ring {' '.join(nei_info)}")
        if len(candidates):
            cand_pairs = sorted(list(zip(candidates, dists)), key=lambda x: x[1])
            preview = ", ".join([f"{fmt_qr(k)} d={d:.3f}" for k, d in cand_pairs[:12]])
            print(f"CAND: frag={q},{r},{i} -> candidates (nearest first): {preview}")

    

# --- Visualization ---
def _fnv1a32_mix(nums):
    h = 2166136261  # FNV-1a 32-bit offset basis
    for n in nums:
        h ^= (n & 0xFFFFFFFF)
        h = (h * 16777619) & 0xFFFFFFFF
    return h

def _wang_hash_32(x):
    x = (x ^ 61) ^ (x >> 16)
    x = (x + (x << 3)) & 0xFFFFFFFF
    x = x ^ (x >> 4)
    x = (x * 0x27d4eb2d) & 0xFFFFFFFF
    x = x ^ (x >> 15)
    return x & 0xFFFFFFFF

def _id_to_color(id_value, cmap=plt.cm.hsv):
    parts = str(id_value).split(',')
    # fold (q,r,i) -> 3 ints; pad with zeros if missing
    try:
        nums = [int(p) for p in parts if p != '']
    except ValueError:
        nums = [0, 0, 0]
    while len(nums) < 3: nums.append(0)
    seed = _fnv1a32_mix(nums[:3])
    h = _wang_hash_32(seed)
    v = h / 0xFFFFFFFF  # scalar in [0,1]
    return cmap(v)

def plot_rasterization_np(polygon_verts, elements, hex_size, draw_edges=True, plot_com_shift=False, plot_merge_lines=False, plot_base_labels=False, plot_frag_labels=False):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.add_patch(MplPolygon(polygon_verts, fill=False, edgecolor='k', linewidth=2.5, zorder=10, label='Original Polygon Boundary'))
    
    # Distinct coloring per element using stable hash of element ID
    cmap = plt.cm.hsv
    added_base_orig = False
    added_base_merged = False
    added_boundary = False
    for idx, (key, element) in enumerate(elements.items()):
        lbl = getattr(element, 'id', str(key))
        color = _id_to_color(lbl, cmap)
        ec = 'k' if draw_edges else 'none'
        lw = 0.5 if draw_edges else 0.0
        for geom_verts in element.geoms:
            ax.add_patch(MplPolygon(geom_verts, facecolor=color, alpha=0.6, edgecolor=ec, linewidth=lw))
        com = element.com
        if element.is_base:
            com0 = element.com_orig if (element.com_orig is not None) else com
            # original COM (blue x), updated COM (red o)
            ax.plot(com0[0], com0[1], 'bx', markersize=6, zorder=6, label=('Base COM (orig)' if not added_base_orig else None))
            added_base_orig = True
            ax.plot(com[0], com[1], 'ro', markersize=4, zorder=6, label=('Base COM (merged)' if not added_base_merged else None))
            added_base_merged = True
            if plot_com_shift and (np.linalg.norm(com - com0) > 1e-12):
                ax.plot([com0[0], com[0]], [com0[1], com[1]], color='r', linestyle='--', linewidth=1.0, alpha=0.9)
            if plot_merge_lines and element.merged_piece_info:
                for info in element.merged_piece_info:
                    pc = info['com']
                    ax.plot([com[0], pc[0]], [com[1], pc[1]], color='gray', linestyle='-', linewidth=0.7, alpha=0.7)
                    # mark merged fragment COMs and optionally label
                    dot_col = _id_to_color(info.get('id', ''), cmap)
                    ax.plot(pc[0], pc[1], marker='.', color=dot_col, markersize=3, alpha=0.9)
                    if plot_frag_labels and ('id' in info):
                        ax.text(pc[0], pc[1], f"F {info['id']}", color='black', fontsize=7, ha='left', va='bottom')
            if plot_base_labels:
                lbl = getattr(element, 'id', str(key))
                ax.text(com[0], com[1], f"B {lbl}", color='black', fontsize=8, ha='center', va='center', zorder=100)
        else:
            ax.plot(com[0], com[1], 'go', markersize=4, zorder=6, label=('Boundary element COM' if not added_boundary else None))
            added_boundary = True
            if plot_frag_labels:
                ax.text(com[0], com[1], f"F {getattr(element, 'id', str(key))}", color='black', fontsize=8, ha='center', va='center', zorder=100)
            if plot_merge_lines:
                # Draw line to nearest base element for debugging relation
                base_keys = [k for k, v in elements.items() if isinstance(k, tuple) and v.is_base]
                if len(base_keys) > 0:
                    base_coms = np.array([elements[k].com for k in base_keys])
                    dists = np.linalg.norm(base_coms - com, axis=1)
                    best_key = base_keys[int(np.argmin(dists))]
                    bc = elements[best_key].com
                    ax.plot([bc[0], com[0]], [bc[1], com[1]], color='gray', linestyle='--', linewidth=0.7, alpha=0.7)

    min_coords, max_coords = np.min(polygon_verts, axis=0) - hex_size, np.max(polygon_verts, axis=0) + hex_size
    ax.set_xlim(min_coords[0], max_coords[0])
    ax.set_ylim(min_coords[1], max_coords[1])
    ax.set_aspect('equal', 'box')
    ax.set_title("Hex Rasterization with Boundary Merging (debug)")
    ax.legend(loc='upper right')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

# --- Main Execution ---
if __name__ == '__main__':
    # run like this:
    #   python PolygonRasterization.py --size 1.0 --areaFactor 0.4 --staticCOM --plotCOMShift --plotMergeLines --plotBaseLabels --plotFragLabels --noEdges
    #   python PolygonRasterization.py --size 1.0 --areaFactor 0.7 --plotMergeLines --plotBaseLabels --plotFragLabels --verbosity 4
    import argparse
    parser = argparse.ArgumentParser(description='Hexagonal rasterization of polygon with area-aware boundary handling')
    parser.add_argument('--size', type=float, default=1.0, help='Hexagon size (grid step)')
    parser.add_argument('--areaFactor', type=float, default=0.7, help='Keep clipped hex if A_clip/A_full >= areaFactor; otherwise split triangles')
    #parser.add_argument('--mergeFactor', type=float, default=0.3, help='Triangle piece merge threshold as fraction of full hex area')
    parser.add_argument('--plotCOMShift', action='store_true', help='Plot line between original base COM and merged COM')
    parser.add_argument('--plotMergeLines', action='store_true', help='Plot lines from merged base COM to merged piece COMs')
    parser.add_argument('--noEdges', action='store_true', help='If set, do not draw polygon edges, only filled areas')
    parser.add_argument('--plotBaseLabels', action='store_true', help='Plot labels for base/final elements')
    parser.add_argument('--plotFragLabels', action='store_true', help='Plot labels for fragment pieces (standalone or merged)')
    parser.add_argument('--verbosity', type=int, default=1, help='Verbosity level; >3 prints candidate neighbors and distances for merges')
    parser.add_argument('--staticCOM', action='store_true', help='Keep base COM static during merges (do not update COM)')
    args = parser.parse_args()

    #subject_polygon_verts =  np.array( [(0, 8), (-2, 2), (-8, 2), (-4, -2), (-5, -8), (0, -5), (5, -8), (4, -2), (8, 2), (2, 2)])
    subject_polygon_verts = np.array([(0, 8),(3, 8),  (8, 2),  (2, 2)])
    
    HEX_SIZE = args.size

    rasterizer = HexRasterizerNP(subject_polygon_verts, HEX_SIZE, verbosity=args.verbosity, update_com=(not args.staticCOM))
    final_elements = rasterizer.rasterize(area_keep_ratio=args.areaFactor)
    
    print(f"\nRasterization complete. Generated {len(final_elements)} final elements.")
    plot_rasterization_np(subject_polygon_verts, final_elements, HEX_SIZE, draw_edges=(not args.noEdges), plot_com_shift=args.plotCOMShift, plot_merge_lines=args.plotMergeLines, plot_base_labels=args.plotBaseLabels, plot_frag_labels=args.plotFragLabels)