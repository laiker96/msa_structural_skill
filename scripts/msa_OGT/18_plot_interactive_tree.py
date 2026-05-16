#!/usr/bin/env python3
"""Step 18: plot the final OGT-aware tree as a rectangular interactive HTML tree."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

cfg = import_module("00_config")

DEFAULT_TREE = cfg.INTER_DIR / "regime_clades" / "regime_clades_rooted_with_outgroups.nwk"
DEFAULT_METADATA = cfg.INTER_DIR / "regime_clades" / "regime_clades_tip_metadata.tsv"
DEFAULT_CLADE_TABLE = cfg.INTER_DIR / "regime_clades" / "regime_clades.tsv"
DEFAULT_OUT_DIR = cfg.INTER_DIR / "interactive_tree"
DEFAULT_HELIX_CALLS = cfg.INTER_DIR / "helix_extension_classifier" / "helix_extension_calls.tsv"
HELIX_PARTIAL_CALLS = {"partial_extension"}
HELIX_FULL_RUN_MIN = 6
HELIX_PARTIAL_RUN_MIN = 3

REGIME_COLORS = {
    "query": "#111827",
    "psychro": "#2563eb",
    "meso": "#16a34a",
    "thermo": "#dc2626",
    "unknown": "#9ca3af",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    return {row["id"]: row for row in read_tsv(path) if row.get("id")}


def helix_signal_for_row(row: dict[str, str]) -> str:
    signal = row.get("extension_signal", "")
    if signal:
        return signal
    call_text = str(row.get("call", "") or "")
    try:
        extension_run_len = int(row.get("extension_run_len", "") or 0)
    except (TypeError, ValueError):
        extension_run_len = 0
    if extension_run_len >= HELIX_FULL_RUN_MIN:
        return "positive"
    if extension_run_len >= HELIX_PARTIAL_RUN_MIN or call_text in HELIX_PARTIAL_CALLS:
        return "partial"
    if call_text in {"missing_structure", "structure_error", "missing_core_region"}:
        return "missing"
    return "negative"


def helix_signal_for_call(call: object) -> str:
    call_text = str(call or "")
    if call_text == "extended_helix":
        return "positive"
    if call_text in HELIX_PARTIAL_CALLS:
        return "partial"
    if call_text in {"missing_structure", "structure_error", "missing_core_region"}:
        return "missing"
    return "negative"


def load_helix_calls(path: Path) -> dict[str, dict[str, str]]:
    calls = {}
    for row in read_tsv(path):
        sid = row.get("id", "")
        if not sid:
            continue
        row = dict(row)
        row["extension_signal"] = helix_signal_for_row(row)
        calls[sid] = row
    return calls


def tip_regime(tip_id: str, row: dict[str, str] | None) -> str:
    if tip_id == "photoHymenobact":
        return "query"
    if not row:
        return "unknown"
    return row.get("threshold_regime") or row.get("regime") or "unknown"


def clade_id(row: dict[str, str] | None) -> str:
    if not row:
        return ""
    return row.get("regime_clade_id") or row.get("regime_clade_context_id") or ""


def ensure_node_names(tree) -> None:
    used = {clade.name for clade in tree.find_clades() if clade.name}
    idx = 1
    for clade in tree.get_nonterminals(order="level"):
        if clade.name:
            continue
        while f"N{idx}" in used:
            idx += 1
        clade.name = f"N{idx}"
        used.add(clade.name)


def cumulative_depths(tree) -> dict[object, float]:
    depths = {}

    def walk(clade, depth: float):
        depths[clade] = depth
        for child in clade.clades:
            walk(child, depth + float(child.branch_length or 0.0))

    walk(tree.root, 0.0)
    return depths


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    pos = fraction * (len(ordered) - 1)
    lower = math.floor(pos)
    upper = math.ceil(pos)
    if lower == upper:
        return ordered[lower]
    weight = pos - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def branch_break_limit(tree, args) -> float | None:
    if args.no_branch_breaks:
        return None
    if args.branch_break_threshold != "auto":
        try:
            threshold = float(args.branch_break_threshold)
        except ValueError as exc:
            raise SystemExit("ERROR: --branch-break-threshold must be 'auto' or a positive number") from exc
        if threshold <= 0:
            raise SystemExit("ERROR: --branch-break-threshold must be positive")
        return threshold

    lengths = [
        float(clade.branch_length)
        for clade in tree.find_clades(order="preorder")
        if clade.branch_length is not None and float(clade.branch_length) > 0
    ]
    if len(lengths) < 8:
        return None
    median = percentile(lengths, 0.5)
    if median <= 0 or max(lengths) <= median * 20:
        return None
    return max(median * 10, percentile(lengths, 0.90))


def display_depths(tree, args) -> tuple[dict[object, float], dict[int, dict[str, float]], float | None]:
    threshold = branch_break_limit(tree, args)
    depths = {}
    breaks = {}

    def displayed_length(length: float) -> float:
        if threshold is None or length <= threshold:
            return length
        return threshold * (1.0 + math.log1p((length - threshold) / threshold))

    def walk(clade, depth: float):
        depths[clade] = depth
        for child in clade.clades:
            raw = max(0.0, float(child.branch_length or 0.0))
            shown = displayed_length(raw)
            if threshold is not None and raw > threshold:
                breaks[id(child)] = {"raw": raw, "shown": shown}
            walk(child, depth + shown)

    walk(tree.root, 0.0)
    return depths, breaks, threshold


def terminal_order(tree) -> dict[object, float]:
    y = {terminal: float(idx) for idx, terminal in enumerate(tree.get_terminals())}

    def assign(clade):
        if clade in y:
            return y[clade]
        child_y = [assign(child) for child in clade.clades]
        y[clade] = sum(child_y) / len(child_y)
        return y[clade]

    assign(tree.root)
    return y


def regime_counts_for_clade(clade, metadata: dict[str, dict[str, str]]) -> Counter:
    counts = Counter()
    for terminal in clade.get_terminals():
        counts[tip_regime(terminal.name, metadata.get(terminal.name))] += 1
    return counts


def helix_counts_for_clade(clade, helix_calls: dict[str, dict[str, str]]) -> Counter:
    counts = Counter()
    for terminal in clade.get_terminals():
        row = helix_calls.get(terminal.name)
        signal = row.get("extension_signal", "") if row else "missing"
        counts[signal or "missing"] += 1
    return counts


def clade_tip_ids(metadata: dict[str, dict[str, str]]) -> dict[str, list[str]]:
    out = defaultdict(list)
    for tip_id, row in metadata.items():
        cid = clade_id(row)
        if cid:
            out[cid].append(tip_id)
    return out


def shorten(value: str, n: int) -> str:
    return value if len(value) <= n else value[: n - 1] + "..."


def afdb_code_from_helix_row(tip_id: str, helix_row: dict[str, str]) -> str:
    if tip_id == "photoHymenobact":
        return "ANKros"
    structure_path = helix_row.get("structure_path", "")
    if structure_path:
        stem = Path(structure_path).name
        if stem.startswith("AF-"):
            return stem.removesuffix(".pdb").removesuffix(".cif").removesuffix(".mmcif")
    accession = helix_row.get("accession", "")
    return f"AF-{accession}-F1" if accession else ""


def format_support(value) -> str:
    if value is None:
        return ""
    try:
        support = float(value)
    except (TypeError, ValueError):
        return ""
    if math.isclose(support, round(support), abs_tol=1e-9):
        return str(int(round(support)))
    return f"{support:.3g}"


HELIX_MARKERS = {
    "positive": {"shape": "star", "color": "#0f766e", "label": "extension positive"},
    "partial": {"shape": "triangle", "color": "#d97706", "label": "partial extension"},
    "negative": {"shape": "square", "color": "#475569", "label": "negative"},
    "missing": {"shape": "circle", "color": "#cbd5e1", "label": "missing"},
}


def polygon_points(points: list[tuple[float, float]]) -> str:
    return " ".join(f"{x:.3f},{y:.3f}" for x, y in points)


def star_points(cx: float, cy: float, outer: float = 7.2, inner: float = 3.25) -> str:
    points = []
    for idx in range(10):
        angle = -math.pi / 2 + idx * math.pi / 5
        radius = outer if idx % 2 == 0 else inner
        points.append((cx + math.cos(angle) * radius, cy + math.sin(angle) * radius))
    return polygon_points(points)


def triangle_points(cx: float, cy: float, radius: float = 6.4) -> str:
    return polygon_points(
        [
            (cx, cy - radius),
            (cx + radius * 0.9, cy + radius * 0.65),
            (cx - radius * 0.9, cy + radius * 0.65),
        ]
    )


def helix_marker_svg(
    signal: str,
    cx: float,
    cy: float,
    attrs: str = "",
    css_class: str | None = None,
    fill: str | None = None,
) -> str:
    marker = HELIX_MARKERS.get(signal, HELIX_MARKERS["missing"])
    color = html.escape(fill or marker["color"])
    css_class = css_class or f"helix-marker helix-{html.escape(signal)}"
    if marker["shape"] == "star":
        return f'<polygon class="{css_class}" {attrs} points="{star_points(cx, cy)}" fill="{color}"/>'
    if marker["shape"] == "triangle":
        return f'<polygon class="{css_class}" {attrs} points="{triangle_points(cx, cy)}" fill="{color}"/>'
    if marker["shape"] == "square":
        size = 10.2
        return f'<rect class="{css_class}" {attrs} x="{cx - size / 2:.3f}" y="{cy - size / 2:.3f}" width="{size:.3f}" height="{size:.3f}" fill="{color}"/>'
    return f'<circle class="{css_class}" {attrs} cx="{cx:.3f}" cy="{cy:.3f}" r="5.4" fill="{color}"/>'


def build_clade_payload(tree, metadata: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    tree_tip_names = {terminal.name for terminal in tree.get_terminals()}
    out = []
    for cid, tips in sorted(clade_tip_ids(metadata).items()):
        keep = sorted(set(tips) & tree_tip_names)
        if not keep:
            continue
        try:
            mrca = tree.common_ancestor(keep)
        except Exception:
            continue
        out.append({
            "clade": cid,
            "node": mrca.name,
            "n_tips": len(keep),
            "regime": cid.split("_", 1)[0] if "_" in cid else "",
        })
    return out


def clade_geometry_summary(row: dict[str, str] | None) -> str:
    if not row:
        return ""
    parts = []
    for domain in ("antenna", "linker", "catalytic"):
        mean_rmsd = row.get(f"{domain}_mean_ca_rmsd_a", "")
        n_pairs = row.get(f"{domain}_geometry_n_pairs", "")
        if mean_rmsd:
            parts.append(f"{domain} mean RMSD {mean_rmsd} A ({n_pairs or 0} pairs)")
    return "; ".join(parts)


def build_html(
    args,
    tree,
    metadata: dict[str, dict[str, str]],
    clade_rows: list[dict[str, str]],
    helix_calls: dict[str, dict[str, str]],
) -> str:
    depths, branch_breaks, branch_break_threshold = display_depths(tree, args)
    y_by_clade = terminal_order(tree)
    terminals = tree.get_terminals()
    nonterminals = tree.get_nonterminals()
    max_depth = max(depths.values()) or 1.0
    clade_payload = build_clade_payload(tree, metadata)
    clade_row_by_id = {row.get("clade_id", ""): row for row in clade_rows if row.get("clade_id")}
    for item in clade_payload:
        row = clade_row_by_id.get(str(item["clade"]), {})
        item["geometry"] = clade_geometry_summary(row)
        item["n_labelled"] = row.get("n_labelled", "")
        item["n_unlabelled"] = row.get("n_unlabelled", "")
        item["target_fraction"] = row.get("target_fraction", "")
        tips = [terminal for terminal in tree.get_terminals() if clade_id(metadata.get(terminal.name)) == item["clade"]]
        hcounts = Counter()
        for terminal in tips:
            hrow = helix_calls.get(terminal.name, {})
            hcounts[hrow.get("extension_signal") or "missing"] += 1
        callable_n = hcounts.get("positive", 0) + hcounts.get("partial", 0) + hcounts.get("negative", 0)
        item["helix_positive"] = hcounts.get("positive", 0)
        item["helix_partial"] = hcounts.get("partial", 0)
        item["helix_negative"] = hcounts.get("negative", 0)
        item["helix_missing"] = hcounts.get("missing", 0)
        item["helix_positive_fraction_callable"] = (
            f"{hcounts.get('positive', 0) / callable_n:.3f}" if callable_n else ""
        )
    clade_node = {item["clade"]: item["node"] for item in clade_payload}
    clades_by_node = defaultdict(list)
    for item in clade_payload:
        clades_by_node[item["node"]].append(item["clade"])

    margin_left = 36
    margin_right = int(args.label_width)
    margin_top = 44
    tip_gap = float(args.tip_gap)
    width = int(args.width)
    tree_width = width - margin_left - margin_right
    height = int(margin_top * 2 + max(1, len(terminals) - 1) * tip_gap)

    def x(clade) -> float:
        return margin_left + depths[clade] / max_depth * tree_width

    def y(clade) -> float:
        return margin_top + y_by_clade[clade] * tip_gap

    branches = []
    stems = []
    break_marks = []
    nodes = []
    labels = []
    support_labels = []
    terminal_payload = []

    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            px, py = x(parent), y(parent)
            cx, cy = x(child), y(child)
            branches.append(f'<path class="branch" d="M {px:.3f} {cy:.3f} L {cx:.3f} {cy:.3f}"/>')
            stems.append(f'<path class="branch" d="M {px:.3f} {py:.3f} L {px:.3f} {cy:.3f}"/>')
            if id(child) in branch_breaks:
                break_marks.append(
                    f'<text class="branch-break" x="{((px + cx) / 2):.3f}" y="{cy:.3f}" '
                    f'data-title="branch length: {branch_breaks[id(child)]["raw"]:.6g}; '
                    f'displayed compressed">{html.escape("//")}</text>'
                )
            support = format_support(child.confidence)
            if support and not child.is_terminal():
                support_labels.append(
                    f'<text class="support-label" data-node="{html.escape(child.name or "")}" '
                    f'x="{(px + cx) / 2:.3f}" y="{cy - 4:.3f}">{html.escape(support)}</text>'
                )

    for clade in nonterminals:
        counts = regime_counts_for_clade(clade, metadata)
        helix_counts = helix_counts_for_clade(clade, helix_calls)
        helix_callable = (
            helix_counts.get("positive", 0)
            + helix_counts.get("partial", 0)
            + helix_counts.get("negative", 0)
        )
        dominant = counts.most_common(1)[0][0] if counts else "unknown"
        clade_tags = sorted(clades_by_node.get(clade.name, []))
        support = format_support(clade.confidence)
        tooltip = [
            f"node: {clade.name}",
            f"support: {support or 'none'}",
            f"descendant tips: {sum(counts.values())}",
            "regimes: " + ", ".join(f"{key}={counts[key]}" for key in sorted(counts)),
            "helix extension: "
            + ", ".join(f"{key}={helix_counts[key]}" for key in sorted(helix_counts))
            + (f" (positive/callable={helix_counts.get('positive', 0)}/{helix_callable})" if helix_callable else ""),
        ]
        if clade_tags:
            tooltip.append("called clades: " + ", ".join(clade_tags))
            for tag in clade_tags:
                geometry = clade_geometry_summary(clade_row_by_id.get(tag))
                if geometry:
                    tooltip.append(f"{tag} geometry: {geometry}")
        title = html.escape("\n".join(tooltip), quote=True)
        radius = 4.8 if clade_tags else 2.5
        cls = "internal-node clade-node" if clade_tags else "internal-node"
        nodes.append(
            f'<circle class="{cls}" data-node="{html.escape(clade.name)}" '
            f'data-regime="{html.escape(dominant)}" data-clades="{html.escape(" ".join(clade_tags))}" '
            f'data-title="{title}" cx="{x(clade):.3f}" cy="{y(clade):.3f}" r="{radius:.1f}"/>'
        )
        labels.append(
            f'<text class="node-label" data-node="{html.escape(clade.name)}" '
            f'x="{x(clade) + 7:.3f}" y="{y(clade) - 5:.3f}">{html.escape(clade.name)}</text>'
        )
        if clade_tags:
            labels.append(
                f'<text class="clade-label" data-node="{html.escape(clade.name)}" '
                f'x="{x(clade) + 7:.3f}" y="{y(clade) + 8:.3f}">{html.escape(",".join(clade_tags[:3]))}</text>'
            )

    for terminal in terminals:
        row = metadata.get(terminal.name)
        regime = tip_regime(terminal.name, row)
        cid = clade_id(row)
        organism = row.get("organism", "") if row else ""
        ogt = row.get("ogt", "") if row else ""
        helix_row = helix_calls.get(terminal.name, {})
        helix_call = helix_row.get("call", "missing_structure" if terminal.name != "photoHymenobact" else "")
        helix_signal = helix_row.get("extension_signal") or helix_signal_for_call(helix_call)
        extension_run_len = helix_row.get("extension_run_len", "")
        extension_qpos_span = helix_row.get("extension_qpos_span", "")
        extension_sequence = helix_row.get("extension_sequence", "")
        afdb_code = afdb_code_from_helix_row(terminal.name, helix_row)
        display = terminal.name if args.tip_label == "id" else organism or terminal.name
        if terminal.name == "photoHymenobact":
            display = "ANKros (photoHymenobact)"
        elif afdb_code:
            display = f"{display} [{afdb_code}]"
        label = shorten(display, args.label_chars)
        tooltip = [
            f"id: {terminal.name}",
            f"AFDB code: {afdb_code or 'none'}",
            f"regime: {regime}",
            f"clade: {cid or 'none'}",
        ]
        if ogt:
            tooltip.append(f"OGT: {ogt}")
        if organism:
            tooltip.append(f"organism: {organism}")
        if helix_call:
            tooltip.append(f"helix extension call: {helix_call}")
            tooltip.append(f"helix extension signal: {helix_signal}")
        if extension_run_len:
            tooltip.append(f"extension run: {extension_run_len} residues")
        if extension_qpos_span:
            tooltip.append(f"extension qpos: {extension_qpos_span}")
        if extension_sequence:
            tooltip.append(f"extension sequence: {extension_sequence}")
        title = html.escape("\n".join(tooltip), quote=True)
        cls = f"tip-node tip-{html.escape(regime)} helix-marker helix-{html.escape(helix_signal)}"
        if terminal.name == "photoHymenobact":
            cls += " key-tip"
        marker_attrs = (
            f'data-tip="{html.escape(terminal.name)}" data-regime="{html.escape(regime)}" '
            f'data-clade="{html.escape(cid)}" data-helix="{html.escape(helix_signal)}" '
            f'data-title="{title}"'
        )
        nodes.append(
            helix_marker_svg(
                helix_signal,
                x(terminal),
                y(terminal),
                marker_attrs,
                css_class=cls,
                fill=REGIME_COLORS.get(regime, REGIME_COLORS["unknown"]),
            )
        )
        label_cls = "tip-label key-tip-label" if terminal.name == "photoHymenobact" else "tip-label"
        labels.append(
            f'<text class="{label_cls}" data-tip="{html.escape(terminal.name)}" '
            f'data-regime="{html.escape(regime)}" data-clade="{html.escape(cid)}" data-helix="{html.escape(helix_signal)}" '
            f'x="{x(terminal) + 10:.3f}" y="{y(terminal) + 3.4:.3f}">{html.escape(label)}</text>'
        )
        terminal_payload.append({
            "id": terminal.name,
            "regime": regime,
            "clade": cid,
            "organism": organism,
            "afdb_code": afdb_code,
            "ogt": ogt,
            "helix_signal": helix_signal,
            "helix_call": helix_call,
        })

    regime_counts = Counter(item["regime"] for item in terminal_payload)
    helix_signal_counts = Counter(item["helix_signal"] for item in terminal_payload)
    branch_break_note = ""
    if branch_breaks:
        threshold_text = "auto" if branch_break_threshold is None else f"{branch_break_threshold:.6g}"
        branch_break_note = (
            f"<div class=\"small\">Branch breaks: {len(branch_breaks)} long branches "
            f"compressed for display; // marks branches above {html.escape(threshold_text)}.</div>"
        )
    clade_options = ["<option value=''>Select clade</option>"] + [
        (
            f'<option value="{html.escape(item["clade"])}" title="{html.escape(str(item.get("geometry", "")), quote=True)}">'
            f'{html.escape(item["clade"])} ({item["n_tips"]}; labelled={html.escape(str(item.get("n_labelled", "")))}'
            f'; unknown={html.escape(str(item.get("n_unlabelled", "")))}'
            f'; helix+={html.escape(str(item.get("helix_positive", "")))})'
            f'{(" - " + html.escape(str(item.get("geometry", "")))) if item.get("geometry") else ""}</option>'
        )
        for item in clade_payload
    ]
    jump_options = ["<option value='photoHymenobact'>ANKros query tip</option>"]
    jump_options.extend(
        f'<option value="{html.escape(clade.name)}">{html.escape(clade.name)}</option>'
        for clade in nonterminals
    )
    jump_options.extend(
        f'<option value="{html.escape(item["id"])}">{html.escape(shorten(((item["organism"] or item["id"]) + (" [" + item["afdb_code"] + "]" if item["afdb_code"] else "")), 120))}</option>'
        for item in sorted(terminal_payload, key=lambda row: (row["organism"] or row["id"]).lower())
    )
    legend_rows = []
    for regime in ("query", "psychro", "meso", "thermo", "unknown"):
        legend_rows.append(
            f'<label><input class="regime-filter" type="checkbox" value="{regime}" checked> '
            f'<span class="swatch" style="background:{REGIME_COLORS[regime]}"></span>{regime}</label>'
            f'<span>{regime_counts.get(regime, 0)}</span>'
        )
    helix_legend_rows = []
    for signal in ("positive", "partial", "negative", "missing"):
        marker = HELIX_MARKERS[signal]
        helix_legend_rows.append(
            f'<label><input class="helix-filter" type="checkbox" value="{signal}" checked> '
            f'<span class="shape-swatch shape-{html.escape(marker["shape"])}" '
            f'style="--marker-color:#111827"></span>{html.escape(marker["label"])}</label>'
            f'<span>{helix_signal_counts.get(signal, 0)}</span>'
        )

    css = f"""
    html, body {{ height:100%; margin:0; overflow:hidden; }}
    body {{ font:13px/1.35 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; color:#111827; background:#f8fafc; }}
    .layout {{ display:grid; grid-template-columns:330px 1fr; height:100vh; overflow:hidden; }}
    aside {{ position:sticky; top:0; height:100vh; overflow:auto; box-sizing:border-box; padding:18px; border-right:1px solid #d1d5db; background:white; }}
    main {{ min-height:0; overflow:hidden; }}
    h1 {{ font-size:17px; margin:0 0 10px; }}
    h2 {{ font-size:12px; margin:18px 0 8px; color:#374151; text-transform:uppercase; letter-spacing:.04em; }}
    input, select, button {{ font:inherit; }}
    input[type=search], select {{ box-sizing:border-box; width:100%; padding:7px 9px; border:1px solid #cbd5e1; border-radius:6px; }}
    button {{ border:1px solid #cbd5e1; border-radius:6px; background:#f9fafb; padding:6px 8px; cursor:pointer; }}
    .row {{ display:flex; gap:8px; align-items:center; flex-wrap:wrap; margin:8px 0; }}
    .small {{ color:#4b5563; font-size:12px; }}
    .legend {{ display:grid; grid-template-columns:1fr auto; gap:5px 10px; }}
    .swatch {{ display:inline-block; width:10px; height:10px; border-radius:50%; margin-right:6px; vertical-align:-1px; }}
    .tree-wrap {{ width:100%; height:100vh; overflow:auto; background:white; scroll-behavior:auto; overscroll-behavior:contain; }}
    svg {{ display:block; background:white; user-select:none; }}
    .branch {{ fill:none; stroke:#9ca3af; stroke-width:1.05; vector-effect:non-scaling-stroke; }}
    .branch-break {{ font-size:12px; font-weight:800; fill:#374151; text-anchor:middle; paint-order:stroke; stroke:white; stroke-width:4px; cursor:default; }}
    .tip-node {{ stroke:white; stroke-width:.7; vector-effect:non-scaling-stroke; }}
    .tip-query {{ fill:{REGIME_COLORS['query']}; }}
    .tip-psychro {{ fill:{REGIME_COLORS['psychro']}; }}
    .tip-meso {{ fill:{REGIME_COLORS['meso']}; }}
    .tip-thermo {{ fill:{REGIME_COLORS['thermo']}; }}
    .tip-unknown {{ fill:{REGIME_COLORS['unknown']}; }}
    .helix-marker {{ stroke:white; stroke-width:1.4; vector-effect:non-scaling-stroke; cursor:default; }}
    .shape-swatch {{ display:inline-block; position:relative; width:14px; height:14px; margin-right:6px; vertical-align:-2px; }}
    .shape-swatch::before {{ content:""; position:absolute; inset:2px; background:var(--marker-color); }}
    .shape-circle::before {{ border-radius:50%; }}
    .shape-square::before {{ border-radius:2px; }}
    .shape-triangle::before {{ inset:1px 0 0 0; width:0; height:0; background:transparent; border-left:7px solid transparent; border-right:7px solid transparent; border-bottom:12px solid var(--marker-color); }}
    .shape-star::before {{ content:"★"; inset:-3px 0 0 0; background:transparent; color:var(--marker-color); font-size:17px; line-height:14px; text-align:center; }}
    .internal-node {{ fill:white; stroke:#6b7280; stroke-width:1.1; vector-effect:non-scaling-stroke; }}
    .clade-node {{ fill:#fef3c7; stroke:#b45309; stroke-width:1.8; }}
    .key-tip {{ stroke:#111827; stroke-width:2.4; }}
    text {{ font-size:10px; fill:#111827; dominant-baseline:middle; }}
    .tip-label {{ display:none; }}
    .labels-on .tip-label {{ display:block; }}
    .key-tip-label {{ display:block; font-weight:800; font-size:12px; paint-order:stroke; stroke:white; stroke-width:4px; }}
    .node-label {{ display:none; font-weight:700; font-size:9px; paint-order:stroke; stroke:white; stroke-width:3px; }}
    .clade-label {{ display:none; font-weight:800; font-size:10px; fill:#92400e; paint-order:stroke; stroke:white; stroke-width:4px; }}
    .clade-labels-on .clade-label {{ display:block; }}
    .support-label {{ display:none; font-weight:700; font-size:9px; fill:#374151; text-anchor:middle; paint-order:stroke; stroke:white; stroke-width:3px; }}
    .support-on .support-label {{ display:block; }}
    .hidden-regime {{ opacity:.08; }}
    .hidden-helix {{ opacity:.06; }}
    .highlight {{ opacity:1 !important; stroke:#f59e0b !important; stroke-width:3 !important; }}
    .highlight-tip {{ opacity:1 !important; stroke:#111827 !important; stroke-width:2.2 !important; }}
    .search-label {{ display:block !important; font-weight:800; paint-order:stroke; stroke:white; stroke-width:4px; }}
    #tooltip {{ position:fixed; display:none; max-width:380px; white-space:pre-wrap; padding:8px 10px; background:#111827; color:white; border-radius:6px; pointer-events:none; z-index:3; box-shadow:0 10px 30px rgba(15,23,42,.24); }}
    """
    js = """
    const data = {terminals: TERMINAL_PAYLOAD, clades: CLADE_PAYLOAD};
    const svg = document.getElementById('treeSvg');
    const treeWrap = document.getElementById('treeWrap');
    const tooltip = document.getElementById('tooltip');
    const search = document.getElementById('search');
    const jump = document.getElementById('jump');
    const cladeSelect = document.getElementById('cladeSelect');
    let viewBox = {x:0, y:0, w:SVG_WIDTH, h:SVG_HEIGHT};
    function clampViewBox(){
      viewBox.w = Math.min(Math.max(viewBox.w, 80), SVG_WIDTH);
      viewBox.h = Math.min(Math.max(viewBox.h, 80), SVG_HEIGHT);
      viewBox.x = Math.min(Math.max(viewBox.x, 0), Math.max(0, SVG_WIDTH - viewBox.w));
      viewBox.y = Math.min(Math.max(viewBox.y, 0), Math.max(0, SVG_HEIGHT - viewBox.h));
    }
    function setViewBox(){ clampViewBox(); svg.setAttribute('viewBox', `${viewBox.x} ${viewBox.y} ${viewBox.w} ${viewBox.h}`); }
    function showTip(evt, el){ const text = el.dataset.title || ''; if(!text) return; tooltip.textContent = text; tooltip.style.display='block'; tooltip.style.left=Math.min(window.innerWidth-420, evt.clientX+14)+'px'; tooltip.style.top=Math.min(window.innerHeight-160, evt.clientY+14)+'px'; }
    function hideTip(){ tooltip.style.display='none'; }
    document.querySelectorAll('[data-title]').forEach(el => { el.addEventListener('mousemove', evt => showTip(evt, el)); el.addEventListener('mouseleave', hideTip); });
    function applyFilters(){ const enabledRegimes = new Set(Array.from(document.querySelectorAll('.regime-filter:checked')).map(x => x.value)); const enabledHelix = new Set(Array.from(document.querySelectorAll('.helix-filter:checked')).map(x => x.value)); document.querySelectorAll('[data-regime]').forEach(el => el.classList.toggle('hidden-regime', !enabledRegimes.has(el.dataset.regime || 'unknown'))); document.querySelectorAll('[data-helix]').forEach(el => el.classList.toggle('hidden-helix', !enabledHelix.has(el.dataset.helix || 'missing'))); }
    document.querySelectorAll('.regime-filter').forEach(el => el.addEventListener('change', applyFilters));
    document.querySelectorAll('.helix-filter').forEach(el => el.addEventListener('change', applyFilters));
    function clearHighlights(){ document.querySelectorAll('.highlight,.highlight-tip,.search-label').forEach(el => el.classList.remove('highlight','highlight-tip','search-label')); }
    function centerOn(el){ if(!el) return; const b=el.getBBox(); const cx=b.x+b.width/2; const cy=b.y+b.height/2; viewBox.w=Math.max(SVG_WIDTH*.22,280); viewBox.h=Math.max(SVG_HEIGHT*.22,220); viewBox.x=cx-viewBox.w/2; viewBox.y=cy-viewBox.h/2; setViewBox(); }
    function markNode(node){ document.querySelectorAll(`[data-node="${CSS.escape(node)}"]`).forEach(el => el.classList.add('highlight')); document.querySelectorAll(`.node-label[data-node="${CSS.escape(node)}"],.clade-label[data-node="${CSS.escape(node)}"]`).forEach(el => el.classList.add('search-label')); }
    function highlightNode(node, center=false){ clearHighlights(); markNode(node); if(center) centerOn(document.querySelector(`[data-node="${CSS.escape(node)}"]`)); }
    function highlightTip(tip, center=false){ clearHighlights(); document.querySelectorAll(`[data-tip="${CSS.escape(tip)}"]`).forEach(el => { if(el.classList.contains('tip-node')) el.classList.add('highlight-tip'); if(el.classList.contains('tip-label')) el.classList.add('search-label'); }); if(center) centerOn(document.querySelector(`.tip-node[data-tip="${CSS.escape(tip)}"]`)); }
    function highlightClade(clade, center=false){ clearHighlights(); document.querySelectorAll(`[data-clade="${CSS.escape(clade)}"]`).forEach(el => { if(el.classList.contains('tip-node')) el.classList.add('highlight-tip'); if(el.classList.contains('tip-label')) el.classList.add('search-label'); }); const c=data.clades.find(x => x.clade === clade); if(c) markNode(c.node); if(center) centerOn(c ? document.querySelector(`[data-node="${CSS.escape(c.node)}"]`) : document.querySelector(`[data-clade="${CSS.escape(clade)}"]`)); }
    function firstHit(q){ const s=q.trim().toLowerCase(); if(!s) return null; if(/^n\\d+$/.test(s)) return {type:'node', id:s.toUpperCase()}; const t=data.terminals.find(x => x.id.toLowerCase()===s || `${x.id} ${x.organism} ${x.afdb_code} ${x.clade} ${x.ogt} ${x.helix_signal} ${x.helix_call}`.toLowerCase().includes(s)); if(t) return {type:'tip', id:t.id}; const c=data.clades.find(x => x.clade.toLowerCase()===s || x.clade.toLowerCase().includes(s)); if(c) return {type:'clade', id:c.clade}; return null; }
    function go(value){ const hit=firstHit(value); if(!hit) return; if(hit.type==='node') highlightNode(hit.id,true); if(hit.type==='tip') highlightTip(hit.id,true); if(hit.type==='clade') highlightClade(hit.id,true); }
    search.addEventListener('keydown', evt => { if(evt.key === 'Enter') go(search.value); });
    search.addEventListener('input', () => { clearHighlights(); const q=search.value.trim().toLowerCase(); if(!q) return; let hits=0; data.terminals.forEach(t => { if(`${t.id} ${t.organism} ${t.afdb_code} ${t.clade} ${t.ogt} ${t.helix_signal} ${t.helix_call}`.toLowerCase().includes(q)){ hits++; document.querySelectorAll(`[data-tip="${CSS.escape(t.id)}"]`).forEach(el => { if(el.classList.contains('tip-node')) el.classList.add('highlight-tip'); if(el.classList.contains('tip-label')) el.classList.add('search-label'); }); }}); document.getElementById('searchHits').textContent = `${hits} tip hits`; });
    document.getElementById('go').addEventListener('click', () => go(jump.value || search.value));
    jump.addEventListener('change', () => go(jump.value));
    cladeSelect.addEventListener('change', () => highlightClade(cladeSelect.value, true));
    document.getElementById('labels').addEventListener('change', evt => svg.classList.toggle('labels-on', evt.target.checked));
    document.getElementById('support').addEventListener('change', evt => svg.classList.toggle('support-on', evt.target.checked));
    document.getElementById('cladeLabels').addEventListener('change', evt => svg.classList.toggle('clade-labels-on', evt.target.checked));
    function scrollTreeStart(){ document.documentElement.scrollTop=0; document.body.scrollTop=0; window.scrollTo(0,0); treeWrap.scrollLeft=0; treeWrap.scrollTop=0; }
    document.getElementById('reset').addEventListener('click', () => { clearHighlights(); search.value=''; jump.value=''; cladeSelect.value=''; viewBox={x:0,y:0,w:SVG_WIDTH,h:SVG_HEIGHT}; setViewBox(); scrollTreeStart(); });
    svg.addEventListener('wheel', evt => { if(!evt.ctrlKey) return; evt.preventDefault(); const scale = evt.deltaY < 0 ? 0.9 : 1.1; const pt=svg.createSVGPoint(); pt.x=evt.clientX; pt.y=evt.clientY; const cursor=pt.matrixTransform(svg.getScreenCTM().inverse()); viewBox.x=cursor.x-(cursor.x-viewBox.x)*scale; viewBox.y=cursor.y-(cursor.y-viewBox.y)*scale; viewBox.w*=scale; viewBox.h*=scale; setViewBox(); }, {passive:false});
    let dragging=null; svg.addEventListener('mousedown', evt => dragging={x:evt.clientX,y:evt.clientY,vb:{...viewBox}}); window.addEventListener('mouseup',()=>dragging=null); window.addEventListener('mousemove', evt => { if(!dragging) return; const sx=viewBox.w/svg.clientWidth; const sy=viewBox.h/svg.clientHeight; viewBox.x=dragging.vb.x-(evt.clientX-dragging.x)*sx; viewBox.y=dragging.vb.y-(evt.clientY-dragging.y)*sy; setViewBox(); });
    if('scrollRestoration' in history) history.scrollRestoration = 'manual';
    window.addEventListener('pageshow', () => setTimeout(scrollTreeStart, 0));
    setViewBox(); applyFilters(); requestAnimationFrame(scrollTreeStart); setTimeout(scrollTreeStart, 50);
    """
    js = (
        js.replace("TERMINAL_PAYLOAD", json.dumps(terminal_payload))
        .replace("CLADE_PAYLOAD", json.dumps(clade_payload))
        .replace("SVG_WIDTH", str(width))
        .replace("SVG_HEIGHT", str(height))
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OGT-Aware MSA Tree</title>
<style>{css}</style>
</head>
<body>
<div class="layout">
<aside>
  <h1>OGT-Aware MSA Tree</h1>
  <div class="small">Tree: {html.escape(str(args.tree))}</div>
  <div class="small">{len(terminals)} tips; {len(nonterminals)} named internal nodes; {len(clade_payload)} called clades</div>
  {branch_break_note}
  <h2>Search</h2>
  <input id="search" type="search" placeholder="Tip, organism, clade, node N1">
  <div class="small" id="searchHits"></div>
  <h2>Go To</h2>
  <div class="row"><input id="jump" list="jumpOptions" type="search" placeholder="Taxon, node, or clade"><button id="go">Go</button></div>
  <datalist id="jumpOptions">{''.join(jump_options)}</datalist>
  <h2>Regime Filter</h2>
  <div class="legend">{''.join(legend_rows)}</div>
  <h2>Helix Extension</h2>
  <div class="legend">{''.join(helix_legend_rows)}</div>
  <h2>Clades</h2>
  <select id="cladeSelect">{''.join(clade_options)}</select>
  <h2>Controls</h2>
  <label class="small"><input id="labels" type="checkbox"> Show all tip labels</label>
  <label class="small"><input id="cladeLabels" type="checkbox"> Show clade labels</label>
  <label class="small"><input id="support" type="checkbox"> Show branch support values</label>
  <div class="row"><button id="reset">Reset</button></div>
  <div class="small">Scroll moves the tree. Ctrl+scroll zooms. Drag pans. Hover tips and nodes for metadata.</div>
</aside>
<main>
  <div class="tree-wrap" id="treeWrap">
    <svg id="treeSvg" width="{width}" height="{height}">
      {''.join(stems)}
      {''.join(branches)}
      {''.join(break_marks)}
      {''.join(nodes)}
      {''.join(support_labels)}
      {''.join(labels)}
    </svg>
  </div>
</main>
</div>
<div id="tooltip"></div>
<script>{js}</script>
</body>
</html>
"""


def write_node_summary(
    path: Path,
    tree,
    metadata: dict[str, dict[str, str]],
    helix_calls: dict[str, dict[str, str]],
) -> None:
    fields = [
        "node",
        "n_tips",
        "query",
        "psychro",
        "meso",
        "thermo",
        "unknown",
        "helix_positive",
        "helix_partial",
        "helix_negative",
        "helix_missing",
        "helix_positive_fraction_callable",
    ]
    rows = []
    for clade in tree.get_nonterminals(order="level"):
        counts = regime_counts_for_clade(clade, metadata)
        helix_counts = helix_counts_for_clade(clade, helix_calls)
        helix_callable = (
            helix_counts.get("positive", 0)
            + helix_counts.get("partial", 0)
            + helix_counts.get("negative", 0)
        )
        row = {"node": clade.name, "n_tips": sum(counts.values())}
        row.update({regime: counts.get(regime, 0) for regime in fields[2:]})
        row.update(
            {
                "helix_positive": helix_counts.get("positive", 0),
                "helix_partial": helix_counts.get("partial", 0),
                "helix_negative": helix_counts.get("negative", 0),
                "helix_missing": helix_counts.get("missing", 0),
                "helix_positive_fraction_callable": (
                    f"{helix_counts.get('positive', 0) / helix_callable:.3f}" if helix_callable else ""
                ),
            }
        )
        rows.append(row)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", default=str(DEFAULT_TREE))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--clades", default=str(DEFAULT_CLADE_TABLE))
    parser.add_argument("--helix-calls", default=str(DEFAULT_HELIX_CALLS))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--width", type=int, default=2800)
    parser.add_argument("--label-width", type=int, default=560)
    parser.add_argument("--tip-gap", type=float, default=18.0)
    parser.add_argument("--tip-label", choices=("id", "organism"), default="organism")
    parser.add_argument("--label-chars", type=int, default=76)
    parser.add_argument(
        "--branch-break-threshold",
        default="auto",
        help="Branch length above which the SVG display is compressed and marked with //. Use 'auto' or a positive number.",
    )
    parser.add_argument(
        "--no-branch-breaks",
        action="store_true",
        help="Draw branch lengths linearly without visual compression.",
    )
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    ensure_node_names(tree)
    metadata = load_metadata(Path(args.metadata))
    clades = read_tsv(Path(args.clades))
    helix_calls = load_helix_calls(Path(args.helix_calls))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    html_path = out_dir / "ogt_msa_tree.html"
    summary_path = out_dir / "ogt_msa_tree_node_summary.tsv"
    named_tree_path = out_dir / "ogt_msa_tree_named_nodes.nwk"

    html_path.write_text(build_html(args, tree, metadata, clades, helix_calls))
    write_node_summary(summary_path, tree, metadata, helix_calls)
    Phylo.write(tree, str(named_tree_path), "newick")

    print(f"Saved: {html_path}")
    print(f"Saved: {summary_path}")
    print(f"Saved: {named_tree_path}")
    print(f"Tips: {len(tree.get_terminals())}; internal nodes: {len(tree.get_nonterminals())}")


if __name__ == "__main__":
    main()
