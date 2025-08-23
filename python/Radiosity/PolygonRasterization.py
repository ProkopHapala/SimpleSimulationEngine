import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import argparse

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

# --- Main Rasterizer Class ---

class HexRasterizerNP:
    def __init__(self, polygon_vertices, hex_size):
        self.polygon_verts = np.array(polygon_vertices)
        self.hex_size = float(hex_size)
        self.full_hex_area = 1.5 * np.sqrt(3) * self.hex_size**2

    def rasterize(self, merge_threshold_ratio=0.25, area_keep_ratio=0.7):
        print("Stage 1: Generating hexagonal grid...")
        hex_centers, hex_verts_map = self._generate_grid()

        print("Stage 2: Classifying hexagons...")
        full_indices, partial_indices = self._classify_hexagons(hex_verts_map)

        print("Stage 3: Creating initial elements from full hexagons...")
        # *** FIX 3: Pass hex_centers to the creation function ***
        elements = self._create_initial_elements(full_indices, hex_verts_map, hex_centers)

        print("Stage 4: Processing partial hexagons (divide-and-conquer + keep-if-large)...")
        self._process_partial_hexagons_granular(partial_indices, hex_centers, hex_verts_map, elements, merge_threshold_ratio, area_keep_ratio)
        
        original_area = polygon_area(self.polygon_verts)
        total_element_area = sum(el['area'] for el in elements.values())
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
            qr: {
                'geoms': [hex_verts_map[qr]],
                'area': self.full_hex_area,
                'com': hex_centers[qr].copy(),
                'com_orig': hex_centers[qr].copy(),
                'is_base': True,
                'merged_pieces': 0,
            } for qr in full_indices
        }

    def _process_partial_hexagons_granular(self, partial_indices, hex_centers, hex_verts_map, elements, merge_threshold_ratio, area_keep_ratio):
        merge_area_threshold = self.full_hex_area * merge_threshold_ratio
        new_element_id_counter = -1
        
        print(f"DEBUG: merge_threshold_area = {merge_area_threshold:.6f} ({merge_threshold_ratio*100:.1f}% of full hex)")
        print(f"DEBUG: area_keep_ratio     = {area_keep_ratio:.3f} (keep clipped hex if A_clip/A_full >= ratio)")
        
        def base_keys(): return [k for k, v in elements.items() if isinstance(k, tuple) and v.get('is_base', False)]
        
        for qr in partial_indices:
            center, vertices, (q, r) = hex_centers[qr], hex_verts_map[qr], qr

            clipped_hex = sutherland_hodgman_clip(self.polygon_verts, vertices)
            if clipped_hex.size > 0 and clipped_hex.shape[0] >= 3:
                a_clip = polygon_area(clipped_hex)
                ratio = a_clip / self.full_hex_area
            else:
                a_clip = 0.0
                ratio = 0.0

            if ratio >= area_keep_ratio and a_clip > 1e-12:
                com_clip = polygon_centroid(clipped_hex)
                elements[qr] = {
                    'geoms': [clipped_hex],
                    'area': a_clip,
                    'com': com_clip,
                    'com_orig': com_clip.copy(),
                    'is_base': True,
                    'merged_pieces': 0,
                }
                print(f"KEEP_HEX: qr={qr} A_clip={a_clip:.6f} ratio={ratio:.3f} -> kept as base element with corrected COM {com_clip}")
                continue

            for i in range(6):
                p1, p2 = vertices[i], vertices[(i + 1) % 6]
                triangle_verts = np.array([center, p1, p2])
                tri_area = polygon_area(triangle_verts)
                clipped_piece = sutherland_hodgman_clip(self.polygon_verts, triangle_verts)

                if clipped_piece.shape[0] < 3:
                    # DEBUG: triangle fully outside
                    print(f"OUT: qr={qr} tri[{i}] A_tri={tri_area:.6f} -> clipped 0")
                    continue
                piece_area = polygon_area(clipped_piece)
                if piece_area < 1e-12:
                    print(f"SMALL: qr={qr} tri[{i}] A_tri={tri_area:.6f} -> clipped ~0")
                    continue
                
                piece_com = polygon_centroid(clipped_piece)
                
                if piece_area < merge_area_threshold:
                    # Merge into nearest base (full) neighbor if exists
                    candidates = base_keys()
                    if len(candidates) == 0:
                        # No base elements at all (unlikely), create standalone
                        elements[new_element_id_counter] = {'geoms': [clipped_piece], 'area': piece_area, 'com': piece_com, 'is_base': False}
                        print(f"NEW:  qr={qr} tri[{i}] A_tri={tri_area:.6f} A_clip={piece_area:.6f} -> elem_id={new_element_id_counter} com={piece_com}")
                        new_element_id_counter -= 1
                        continue
                    
                    # Find nearest base element by COM distance
                    base_coms = np.array([elements[k]['com'] for k in candidates])
                    dists = np.linalg.norm(base_coms - piece_com, axis=1)
                    best_idx = int(np.argmin(dists))
                    best_key = candidates[best_idx]
                    neighbor_el = elements[best_key]
                    old_area = neighbor_el['area']
                    old_com = neighbor_el['com']
                    new_area = old_area + piece_area
                    new_com = (old_com * old_area + piece_com * piece_area) / new_area
                    neighbor_el['area'] = new_area
                    neighbor_el['com'] = new_com
                    neighbor_el['geoms'].append(clipped_piece)
                    neighbor_el['merged_pieces'] = neighbor_el.get('merged_pieces', 0) + 1
                    print(f"MERG: qr={qr} tri[{i}] A_tri={tri_area:.6f} A_clip={piece_area:.6f} -> base{best_key} A:{old_area:.6f}->{new_area:.6f} COM:{old_com}->{new_com}")
                else:
                    # Create an independent element
                    elements[new_element_id_counter] = {'geoms': [clipped_piece], 'area': piece_area, 'com': piece_com, 'is_base': False}
                    print(f"KEEP: qr={qr} tri[{i}] A_tri={tri_area:.6f} A_clip={piece_area:.6f} -> elem_id={new_element_id_counter} com={piece_com}")
                    new_element_id_counter -= 1

# --- Visualization ---
def plot_rasterization_np(polygon_verts, elements, hex_size):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.add_patch(MplPolygon(polygon_verts, fill=False, edgecolor='k', linewidth=2.5, zorder=10, label='Original Polygon Boundary'))
    
    # Distinct coloring per element using evenly spaced hues (HSV)
    N = max(len(elements), 1)
    cmap = plt.cm.hsv
    added_base_orig = False
    added_base_merged = False
    added_boundary = False
    for idx, (key, element) in enumerate(elements.items()):
        color = cmap(idx / N)
        for geom_verts in element['geoms']:
            ax.add_patch(MplPolygon(geom_verts, facecolor=color, alpha=0.6, edgecolor='k', linewidth=0.5))
        com = element['com']
        if element.get('is_base', False):
            com0 = element.get('com_orig', com)
            # original COM (blue x), updated COM (red o)
            ax.plot(com0[0], com0[1], 'bx', markersize=6, zorder=6, label=('Base COM (orig)' if not added_base_orig else None))
            added_base_orig = True
            ax.plot(com[0], com[1], 'ro', markersize=4, zorder=6, label=('Base COM (merged)' if not added_base_merged else None))
            added_base_merged = True
        else:
            ax.plot(com[0], com[1], 'go', markersize=4, zorder=6, label=('Boundary element COM' if not added_boundary else None))
            added_boundary = True

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
    import argparse
    parser = argparse.ArgumentParser(description='Hexagonal rasterization of polygon with area-aware boundary handling')
    parser.add_argument('--size', type=float, default=1.0, help='Hexagon size (grid step)')
    parser.add_argument('--areaFactor', type=float, default=0.7, help='Keep clipped hex if A_clip/A_full >= areaFactor; otherwise split triangles')
    parser.add_argument('--mergeFactor', type=float, default=0.3, help='Triangle piece merge threshold as fraction of full hex area')
    args = parser.parse_args()

    #subject_polygon_verts =  np.array( [(0, 8), (-2, 2), (-8, 2), (-4, -2), (-5, -8), (0, -5), (5, -8), (4, -2), (8, 2), (2, 2)])
    subject_polygon_verts = np.array([(0, 8),(3, 8),  (8, 2),  (2, 2)])
    
    HEX_SIZE = args.size

    rasterizer = HexRasterizerNP(subject_polygon_verts, HEX_SIZE)
    final_elements = rasterizer.rasterize(merge_threshold_ratio=args.mergeFactor, area_keep_ratio=args.areaFactor)
    
    print(f"\nRasterization complete. Generated {len(final_elements)} final elements.")
    plot_rasterization_np(subject_polygon_verts, final_elements, HEX_SIZE)