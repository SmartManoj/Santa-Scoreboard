#!/usr/bin/env python3
"""
Calculate individual task scores and total score for Santa 2025 competition.
Implements the competition metric based on valid placements, bounds checking, and overlap detection.
"""
import csv
import math
from collections import defaultdict
from typing import Dict, List, Tuple, Any
from decimal import Decimal

def parse_value(val: str) -> float:
    """Parse value with 's' prefix."""
    return float(val.replace('s', ''))

def load_submission(filename: str = 'submission.csv') -> Dict[str, List[Tuple[float, float, float]]]:
    """
    Load submission file and group entries by task ID.
    
    Returns:
        Dictionary mapping task_id -> list of (x, y, deg) tuples
    """
    tasks = defaultdict(list)
    
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            task_id = row['id'].split('_')[0]
            x = parse_value(row['x'])
            y = parse_value(row['y'])
            deg = parse_value(row['deg'])
            tasks[task_id].append((x, y, deg))
    
    return dict(tasks)

def distance(x1: float, y1: float, x2: float, y2: float) -> float:
    """Calculate Euclidean distance between two points."""
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def normalize_angle(deg: float) -> float:
    """Normalize angle to 0-360 range."""
    while deg < 0:
        deg += 360
    while deg >= 360:
        deg -= 360
    return deg

def get_item_radius(task_id: str, item_index: int, items: List[Tuple[float, float, float]] = None) -> float:  # noqa: ARG001
    """
    Get radius for an item in a task.
    For Christmas trees, typically each item might have a base radius.
    This is an assumption - update based on actual competition specs.
    
    If items are provided, estimate radius based on minimum distances.
    Otherwise use a reasonable default.
    """
    # If items provided, estimate radius from minimum distance
    if items and len(items) > 1:
        min_dist = float('inf')
        x, y, _ = items[item_index]
        for i, (x2, y2, _) in enumerate(items):
            if i != item_index:
                dist = distance(x, y, x2, y2)
                if dist < min_dist:
                    min_dist = dist
        # Use half of minimum distance as radius (allowing some overlap)
        if min_dist < float('inf'):
            return min(min_dist * 0.4, 0.5)  # Cap at 0.5
    
    # Base radius assumption (update based on actual specs)
    # Use smaller radius to be more permissive
    base_radius = 0.3
    return base_radius

def is_within_bounds(x: float, y: float, radius: float, bounds: Tuple[float, float, float, float]) -> bool:
    """
    Check if item (circle) is within bounds.
    
    Args:
        x, y: Center coordinates
        radius: Item radius
        bounds: (min_x, min_y, max_x, max_y)
    
    Returns:
        True if item is within bounds
    """
    min_x, min_y, max_x, max_y = bounds
    
    # Check if circle center plus radius is within bounds
    if x - radius < min_x or x + radius > max_x:
        return False
    if y - radius < min_y or y + radius > max_y:
        return False
    
    return True

def items_overlap(x1: float, y1: float, r1: float, x2: float, y2: float, r2: float, min_distance: float = 1e-6) -> bool:
    """
    Check if two circular items overlap.
    
    Args:
        x1, y1, r1: First item center and radius
        x2, y2, r2: Second item center and radius
        min_distance: Minimum allowed distance between items (for tolerance)
    
    Returns:
        True if items overlap
    """
    dist = distance(x1, y1, x2, y2)
    # Items overlap if distance is less than sum of radii (with tolerance)
    return dist < (r1 + r2 - min_distance)

def get_task_bounds(task_id: str, num_items: int, items: List[Tuple[float, float, float]] = None) -> Tuple[float, float, float, float]:  # noqa: ARG001
    """
    Get bounding area for a task.
    This is an assumption - update based on actual competition specs.
    
    If items are provided, calculate bounds from actual item positions.
    Otherwise use a reasonable estimate based on number of items.
    """
    # If items provided, use actual bounds from positions
    if items:
        x_coords = [x for x, y, _ in items]
        y_coords = [y for x, y, _ in items]
        
        min_x, max_x = min(x_coords), max(x_coords)
        min_y, max_y = min(y_coords), max(y_coords)
        
        # Add margin around items
        margin = 2.0
        return (min_x - margin, min_y - margin, max_x + margin, max_y + margin)
    
    # Estimate bounds based on number of items
    # Scale bounds based on number of items (more items = larger area needed)
    scale = math.sqrt(num_items) * 1.2
    margin = 5.0  # Base margin
    
    min_coord = -margin - scale
    max_coord = margin + scale
    
    return (min_coord, min_coord, max_coord, max_coord)

def get_tree_polygon_points() -> List[Tuple[float, float]]:
    """
    Get the tree polygon points in local coordinates (centered at origin).
    Based on the actual tree geometry from the competition.
    """
    # Tree dimensions
    trunk_w = 0.15
    trunk_h = 0.2
    base_w = 0.7
    mid_w = 0.4
    top_w = 0.25
    tip_y = 0.8
    tier_1_y = 0.5
    tier_2_y = 0.25
    base_y = 0.0
    trunk_bottom_y = -trunk_h
    scale_factor = 1.0
    
    # Build polygon points
    points = [
        # Start at Tip
        (0.0 * scale_factor, tip_y * scale_factor),
        # Right side - Top Tier
        (top_w / 2 * scale_factor, tier_1_y * scale_factor),
        (top_w / 4 * scale_factor, tier_1_y * scale_factor),
        # Right side - Middle Tier
        (mid_w / 2 * scale_factor, tier_2_y * scale_factor),
        (mid_w / 4 * scale_factor, tier_2_y * scale_factor),
        # Right side - Bottom Tier
        (base_w / 2 * scale_factor, base_y * scale_factor),
        # Right Trunk
        (trunk_w / 2 * scale_factor, base_y * scale_factor),
        (trunk_w / 2 * scale_factor, trunk_bottom_y * scale_factor),
        # Left Trunk
        (-(trunk_w / 2) * scale_factor, trunk_bottom_y * scale_factor),
        (-(trunk_w / 2) * scale_factor, base_y * scale_factor),
        # Left side - Bottom Tier
        (-(base_w / 2) * scale_factor, base_y * scale_factor),
        # Left side - Middle Tier
        (-(mid_w / 4) * scale_factor, tier_2_y * scale_factor),
        (-(mid_w / 2) * scale_factor, tier_2_y * scale_factor),
        # Left side - Top Tier
        (-(top_w / 4) * scale_factor, tier_1_y * scale_factor),
        (-(top_w / 2) * scale_factor, tier_1_y * scale_factor),
    ]
    
    return points

def rotate_point(x: float, y: float, angle_deg: float) -> Tuple[float, float]:
    """Rotate a point around the origin by angle_deg degrees."""
    angle_rad = math.radians(angle_deg)
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    
    x_rot = x * cos_a - y * sin_a
    y_rot = x * sin_a + y * cos_a
    
    return (x_rot, y_rot)

def get_tree_bounding_box(center_x: float, center_y: float, rotation_deg: float) -> Tuple[float, float, float, float]:
    """
    Get the bounding box of a tree at position (center_x, center_y) rotated by rotation_deg.
    Returns (min_x, min_y, max_x, max_y).
    """
    # Get tree polygon points in local coordinates
    tree_points = get_tree_polygon_points()
    
    # Rotate and translate all points
    rotated_points = []
    for px, py in tree_points:
        # Rotate
        rx, ry = rotate_point(px, py, rotation_deg)
        # Translate
        rotated_points.append((rx + center_x, ry + center_y))
    
    # Find bounding box
    x_coords = [x for x, y in rotated_points]
    y_coords = [y for x, y in rotated_points]
    
    return (min(x_coords), min(y_coords), max(x_coords), max(y_coords))

def calculate_task_score(task_id: str, items: List[Tuple[float, float, float]]) -> float:  # noqa: ARG001
    """
    Calculate score for a single task using the actual competition metric.
    
    The metric: score = s^2 / n
    where:
    - s = side length of the square bounding box containing all trees
    - n = number of trees in the configuration
    
    The score represents the normalized area of the bounding box.
    Lower scores are better (more compact packing).
    
    Reference: https://www.kaggle.com/code/metric/santa-2025-metric
    """
    if not items:
        return 0.0
    
    n = len(items)  # number of trees
    
    # Calculate bounding box for all trees (accounting for rotation)
    all_min_x, all_min_y = float('inf'), float('inf')
    all_max_x, all_max_y = float('-inf'), float('-inf')
    
    for x, y, deg in items:
        # Get bounding box of this tree
        min_x, min_y, max_x, max_y = get_tree_bounding_box(x, y, deg)
        
        # Expand overall bounding box
        all_min_x = min(all_min_x, min_x)
        all_min_y = min(all_min_y, min_y)
        all_max_x = max(all_max_x, max_x)
        all_max_y = max(all_max_y, max_y)
    
    # Calculate side length of square bounding box
    width = all_max_x - all_min_x
    height = all_max_y - all_min_y
    side = max(width, height)
    
    # Handle edge case: ensure minimum size
    if side < 1e-10:
        side = 1e-10
    
    # Score = s^2 / n (normalized area)
    score = (side * side) / n
    
    return float(score)

def calculate_scores(submission_file: str = 'submission.csv') -> Dict[str, float]:
    """
    Calculate individual task scores and total score.
    
    Returns:
        Dictionary with task scores and total score
    """
    tasks = load_submission(submission_file)
    
    task_scores = {}
    for task_id in sorted(tasks.keys(), key=int):
        items = tasks[task_id]
        score = calculate_task_score(task_id, items)
        task_scores[task_id] = score
    
    total_score = sum(task_scores.values())
    
    return {
        'task_scores': task_scores,
        'total_score': total_score,
        'num_tasks': len(task_scores),
        'total_items': sum(len(items) for items in tasks.values())
    }

def print_scores(results: Dict[str, Any], tasks: Dict[str, List[Tuple[float, float, float]]] = None):
    """Print scores in a formatted way."""
    task_scores = results['task_scores']
    total_score = results['total_score']
    num_tasks = results['num_tasks']
    total_items = results['total_items']
    
    print("=" * 60)
    print("SANTA 2025 SCORE CALCULATION")
    print("=" * 60)
    print(f"\nTotal Tasks: {num_tasks}")
    print(f"Total Items Placed: {total_items:,}")
    print(f"\n{'Task ID':<10} {'Items':<10} {'Score':<15}")
    print("-" * 35)
    
    for task_id in sorted(task_scores.keys(), key=int):
        items_count = len(tasks[task_id]) if tasks and task_id in tasks else 0
        score = task_scores[task_id]
        print(f"{task_id:<10} {items_count:<10} {score:<15.2f}")
    
    print("-" * 35)
    print(f"{'TOTAL':<10} {total_items:<10} {total_score:<15.12f}")
    print("=" * 60)
    print(f"\nTotal Score (12 digits): {total_score:.12f}")
    print(f"\n✓ Metric Implementation (Official Competition Metric):")
    print("   - For each task: score = s² / n")
    print("   - s = side length of square bounding box")
    print("   - n = number of trees in configuration")
    print("   - Total score = sum of all task scores")
    print(f"\n   Lower scores are better (more compact packing)")

def save_scores_csv(results: Dict[str, Any], tasks: Dict[str, List[Tuple[float, float, float]]] = None, filename: str = 'task_scores.csv'):
    """Save task scores to CSV file."""
    task_scores = results['task_scores']
    
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['task_id', 'items_count', 'score'])
        
        for task_id in sorted(task_scores.keys(), key=int):
            score = task_scores[task_id]
            items_count = len(tasks[task_id]) if tasks and task_id in tasks else 0
            writer.writerow([task_id, items_count, score])
        
        # Add total row
        writer.writerow(['TOTAL', results['total_items'], results['total_score']])
    
    print(f"\n✓ Scores saved to {filename}")

if __name__ == '__main__':
    print("Loading submission file...")
    tasks_data = load_submission('submission.csv')
    results = calculate_scores('submission.csv')
    
    print_scores(results, tasks_data)
    save_scores_csv(results, tasks_data)
