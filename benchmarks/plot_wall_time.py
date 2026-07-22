#!/usr/bin/env python3

"""Render a dependency-free SVG wall-time plot from a benchmark summary."""

import argparse
import csv
import math
from pathlib import Path


WIDTH = 900
HEIGHT = 520
LEFT = 86
RIGHT = 80
TOP = 58
BOTTOM = 78


def read_points(path):
    points = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row["run"].startswith("threads-"):
                continue
            points.append(
                {
                    "threads": int(row["run"].removeprefix("threads-")),
                    "seconds": float(row["real_seconds"]),
                }
            )
    return sorted(points, key=lambda point: point["threads"])


def render(points):
    if len(points) < 2:
        raise ValueError("At least two thread measurements are required")

    plot_width = WIDTH - LEFT - RIGHT
    plot_height = HEIGHT - TOP - BOTTOM
    x_min = math.log2(points[0]["threads"])
    x_max = math.log2(points[-1]["threads"])
    minute_values = [point["seconds"] / 60 for point in points]
    y_min = math.floor(min(minute_values) / 5) * 5
    y_max = math.ceil(max(minute_values) / 5) * 5

    def x_position(threads):
        return LEFT + (math.log2(threads) - x_min) / (x_max - x_min) * plot_width

    def y_position(minutes):
        return TOP + (y_max - minutes) / (y_max - y_min) * plot_height

    best = min(points, key=lambda point: point["seconds"])
    first = points[0]
    last = points[-1]
    description = (
        f"Wall time falls from {first['seconds'] / 60:.2f} minutes with "
        f"{first['threads']} requested thread to {best['seconds'] / 60:.2f} "
        f"minutes at {best['threads']} threads. The final measurement at "
        f"{last['threads']} threads is {last['seconds'] / 60:.2f} minutes."
    )
    path_data = " ".join(
        ("M" if index == 0 else "L")
        + f" {x_position(point['threads']):.1f} {y_position(point['seconds'] / 60):.1f}"
        for index, point in enumerate(points)
    )

    output = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}" role="img" aria-labelledby="title description">',
        '<title id="title">TelSeq Parallel wall time versus requested threads</title>',
        f'<desc id="description">{description}</desc>',
        "<style>",
        ".text{fill:#24292f;font:14px -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif}",
        ".muted{fill:#57606a;font:13px -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif}",
        ".title{fill:#24292f;font:600 20px -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif}",
        ".grid{stroke:#d0d7de;stroke-width:1}.axis{stroke:#57606a;stroke-width:1.5}",
        ".series{fill:none;stroke:#0969da;stroke-width:3}.point{fill:#0969da;stroke:#fff;stroke-width:2}",
        ".best{fill:#1a7f37;stroke:#fff;stroke-width:3}.guide{stroke:#1a7f37;stroke-width:1.5;stroke-dasharray:6 5}",
        "@media (prefers-color-scheme:dark){.text,.title{fill:#f0f6fc}.muted{fill:#8c959f}.grid{stroke:#30363d}.axis{stroke:#8c959f}.series{stroke:#58a6ff}.point{fill:#58a6ff;stroke:#0d1117}.best{fill:#3fb950;stroke:#0d1117}.guide{stroke:#3fb950}}",
        "</style>",
        '<text class="title" x="86" y="31">Real WGS wall time by requested thread count</text>',
        '<text class="muted" x="86" y="50">Lower is better · x-axis uses a log₂ scale</text>',
    ]

    for value in range(int(y_min), int(y_max) + 1, 5):
        y = y_position(value)
        output.append(f'<line class="grid" x1="{LEFT}" y1="{y:.1f}" x2="{WIDTH - RIGHT}" y2="{y:.1f}"/>')
        output.append(f'<text class="muted" x="{LEFT - 12}" y="{y + 5:.1f}" text-anchor="end">{value}</text>')

    output.extend(
        [
            f'<line class="axis" x1="{LEFT}" y1="{TOP}" x2="{LEFT}" y2="{HEIGHT - BOTTOM}"/>',
            f'<line class="axis" x1="{LEFT}" y1="{HEIGHT - BOTTOM}" x2="{WIDTH - RIGHT}" y2="{HEIGHT - BOTTOM}"/>',
            f'<line class="guide" x1="{LEFT}" y1="{y_position(best["seconds"] / 60):.1f}" x2="{WIDTH - RIGHT}" y2="{y_position(best["seconds"] / 60):.1f}"/>',
            f'<path class="series" d="{path_data}"/>',
        ]
    )

    for point in points:
        threads = point["threads"]
        minutes = point["seconds"] / 60
        x = x_position(threads)
        y = y_position(minutes)
        point_class = "best" if point is best else "point"
        label_y = y - 15 if threads != 1 else y + 26
        output.append(f'<circle class="{point_class}" cx="{x:.1f}" cy="{y:.1f}" r="6"/>')
        output.append(f'<text class="text" x="{x:.1f}" y="{label_y:.1f}" text-anchor="middle">{minutes:.2f}</text>')
        output.append(f'<text class="muted" x="{x:.1f}" y="{HEIGHT - BOTTOM + 25}" text-anchor="middle">{threads}</text>')

    best_x = x_position(best["threads"])
    best_y = y_position(best["seconds"] / 60)
    annotation_x = best_x - 18 if best_x > WIDTH - RIGHT - 150 else best_x + 18
    annotation_anchor = "end" if best_x > WIDTH - RIGHT - 150 else "start"
    annotation_y = (
        best_y - 28
        if best_y > HEIGHT - BOTTOM - 45
        else best_y + 29
    )
    output.extend(
        [
            f'<text class="text" x="{annotation_x:.1f}" y="{annotation_y:.1f}" text-anchor="{annotation_anchor}">fastest observed</text>',
            f'<text class="text" x="{(LEFT + WIDTH - RIGHT) / 2:.1f}" y="{HEIGHT - 20}" text-anchor="middle">Requested threads</text>',
            f'<text class="text" transform="translate(24 {(TOP + HEIGHT - BOTTOM) / 2:.1f}) rotate(-90)" text-anchor="middle">Wall time (minutes)</text>',
            "</svg>",
        ]
    )
    return "\n".join(output) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    points = read_points(args.summary)
    args.output.write_text(render(points), encoding="utf-8")


if __name__ == "__main__":
    main()
