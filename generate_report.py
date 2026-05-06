#!/usr/bin/env python3
"""
Build a self-contained HTML report directly from per-record stats.yaml files.

By default this script is intended to be executed from the directory that
contains `.looper.yaml`. It reads `output_dir` from that config and writes the
report to `<output_dir>/report.html`.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import re
import shutil
from pathlib import Path

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - depends on runtime environment
    yaml = None


STATUS_ORDER = ("running", "completed", "failed", "waiting", "partial")
STATUS_PATTERN_TEMPLATE = r"^{pipeline}_(?P<record>.+)_(?P<status>{statuses})\.flag$"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the NGS.analysis HTML report from a looper project directory."
    )
    parser.add_argument(
        "--looper-config",
        default=None,
        help="Optional path to .looper.yaml. Defaults to ./.looper.yaml.",
    )
    parser.add_argument(
        "--results-dir",
        default=None,
        help="Optional override for the results directory. By default it is read from .looper.yaml.",
    )
    parser.add_argument(
        "--pipeline-name",
        default=None,
        help="Optional pipeline name override. By default it is inferred from stats.yaml.",
    )
    return parser.parse_args()


def _parse_scalar(value: str) -> object:
    if value == "{}":
        return {}
    if value in {"[]"}:
        return []
    if value in {"null", "Null", "NULL"}:
        return None
    if value in {"true", "True"}:
        return True
    if value in {"false", "False"}:
        return False
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    try:
        if "." in value:
            return float(value)
        return int(value)
    except ValueError:
        return value


def simple_yaml_load(text: str) -> dict:
    root: dict = {}
    stack: list[tuple[int, dict]] = [(-1, root)]
    for raw_line in text.splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        line = raw_line.strip()
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        while stack and indent <= stack[-1][0]:
            stack.pop()
        current = stack[-1][1]
        if value == "":
            current[key] = {}
            stack.append((indent, current[key]))
        else:
            current[key] = _parse_scalar(value)
    return root


def load_yaml_file(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    if yaml is not None:
        return yaml.safe_load(text) or {}
    return simple_yaml_load(text)


def load_looper_config(looper_config_path: Path) -> dict:
    if not looper_config_path.exists():
        raise SystemExit(f"Looper config not found: {looper_config_path}")
    return load_yaml_file(looper_config_path)


def resolve_results_dir(looper_dir: Path, args: argparse.Namespace) -> Path:
    if args.results_dir:
        return (looper_dir / args.results_dir).resolve() if not Path(args.results_dir).is_absolute() else Path(args.results_dir).resolve()
    config_path = looper_dir / ".looper.yaml" if args.looper_config is None else Path(args.looper_config).resolve()
    config = load_looper_config(config_path)
    output_dir = config.get("output_dir")
    if not output_dir:
        raise SystemExit("Could not determine output_dir from .looper.yaml")
    output_path = Path(output_dir)
    if output_path.is_absolute():
        return output_path.resolve()
    return (config_path.parent / output_path).resolve()


def infer_pipeline_name(results_dir: Path, override: str | None) -> str:
    if override:
        return override
    for stats_path in sorted(results_dir.rglob("stats.yaml")):
        data = load_yaml_file(stats_path)
        if data:
            return next(iter(data.keys()))
    raise SystemExit(f"Could not infer pipeline name from stats files in {results_dir}")


def discover_statuses(flags_dir: Path, pipeline_name: str) -> dict[str, str]:
    statuses: dict[str, str] = {}
    if not flags_dir.exists():
        return statuses
    status_group = "|".join(re.escape(status) for status in STATUS_ORDER)
    pattern = re.compile(
        STATUS_PATTERN_TEMPLATE.format(
            pipeline=re.escape(pipeline_name), statuses=status_group
        )
    )
    for flag_path in flags_dir.iterdir():
        match = pattern.match(flag_path.name)
        if match:
            statuses[match.group("record")] = match.group("status")
    return statuses


def read_stats_file(path: Path, pipeline_name: str) -> tuple[dict[str, dict], dict[str, dict]]:
    data = load_yaml_file(path)
    pipeline_data = data.get(pipeline_name, {}) or {}
    return (
        pipeline_data.get("sample", {}) or {},
        pipeline_data.get("project", {}) or {},
    )


def collect_records(results_dir: Path, pipeline_name: str) -> tuple[dict[str, dict], dict[str, dict]]:
    sample_records: dict[str, dict] = {}
    project_records: dict[str, dict] = {}
    for stats_path in sorted(results_dir.rglob("stats.yaml")):
        samples, projects = read_stats_file(stats_path, pipeline_name)
        sample_records.update(samples)
        project_records.update(projects)
    return sample_records, project_records


def split_record_values(record: dict) -> tuple[dict[str, object], dict[str, dict]]:
    scalars: dict[str, object] = {}
    objects: dict[str, dict] = {}
    for key, value in (record or {}).items():
        if key == "meta":
            continue
        if isinstance(value, dict) and "path" in value:
            objects[key] = value
        else:
            scalars[key] = value
    return scalars, objects


def is_number(value: object) -> bool:
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return True
    if isinstance(value, str):
        try:
            float(value)
            return True
        except ValueError:
            return False
    return False


def numeric_value(value: object) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    return float(str(value))


def human_title(name: str) -> str:
    return name.replace("_", " ")


def _suffix_relative_to_results(target: str | Path, results_dir: Path) -> Path | None:
    target_str = str(target).replace("\\", "/")
    results_name = results_dir.name
    marker = f"/{results_name}/"
    if marker in target_str:
        suffix = target_str.split(marker, 1)[1]
        return Path(*[part for part in suffix.split("/") if part])
    return None


def rel_href(from_dir: Path, target: str | Path, results_dir: Path | None = None) -> str:
    target_path = Path(target)
    try:
        return os.path.relpath(target_path, from_dir).replace("\\", "/")
    except ValueError:
        if results_dir is not None:
            suffix = _suffix_relative_to_results(target, results_dir)
            if suffix is not None:
                return os.path.relpath(results_dir / suffix, from_dir).replace("\\", "/")
        return str(target_path).replace("\\", "/")


def infer_thumbnail_path(path_str: str | None, thumbnail_path: str | None) -> str | None:
    if thumbnail_path:
        return thumbnail_path
    if not path_str:
        return None
    path = Path(path_str)
    if path.suffix.lower() == ".pdf":
        candidate = path.with_suffix(".png")
        return str(candidate)
    return None


def object_entries(
    objects: dict[str, dict], page_dir: Path, results_dir: Path
) -> tuple[list[dict], list[dict]]:
    figures: list[dict] = []
    links: list[dict] = []
    for key, value in objects.items():
        path = value.get("path")
        if not path:
            continue
        entry = {
            "key": key,
            "title": value.get("title") or human_title(key),
            "path": rel_href(page_dir, path, results_dir),
        }
        thumbnail_path = infer_thumbnail_path(path, value.get("thumbnail_path"))
        if thumbnail_path:
            entry["thumbnail_path"] = rel_href(page_dir, thumbnail_path, results_dir)
            figures.append(entry)
        else:
            links.append(entry)
    return figures, links


def collect_sample_object_catalog(
    sample_objects_by_sample: dict[str, dict[str, dict]], page_dir: Path, results_dir: Path
) -> tuple[dict[str, dict], list[str]]:
    catalog: dict[str, dict] = {}
    ordered_keys: list[str] = []
    for sample_name, sample_objects in sample_objects_by_sample.items():
        figures, links = object_entries(sample_objects, page_dir, results_dir)
        for item in figures + links:
            key = item["key"]
            if key not in catalog:
                catalog[key] = {"title": item["title"], "items": []}
                ordered_keys.append(key)
            catalog[key]["items"].append(
                {
                    "sample": sample_name,
                    "title": item["title"],
                    "path": item["path"],
                    "thumbnail_path": item.get("thumbnail_path"),
                }
            )
    return catalog, ordered_keys


def record_dir(results_dir: Path, record_name: str) -> Path:
    candidates = (
        results_dir / record_name,
        results_dir / "samples" / record_name,
        results_dir / "project" if record_name == "project" else None,
    )
    for candidate in candidates:
        if candidate is not None and candidate.exists():
            return candidate
    for stats_path in sorted(results_dir.rglob("stats.yaml")):
        if stats_path.parent.name == record_name:
            return stats_path.parent
    return results_dir / record_name


def find_log_link(results_dir: Path, record_name: str, pipeline_name: str, page_dir: Path) -> str | None:
    log_path = record_dir(results_dir, record_name) / f"{pipeline_name}_log.md"
    if log_path.exists():
        return rel_href(page_dir, log_path, results_dir)
    return None


def collect_sample_columns(sample_scalars: dict[str, dict[str, object]]) -> list[str]:
    columns: list[str] = []
    for record in sample_scalars.values():
        for key in record:
            if key not in columns:
                columns.append(key)
    preferred = ["Success", "Time"]
    ordered = [key for key in preferred if key in columns]
    ordered.extend(key for key in columns if key not in ordered)
    return ordered


def collect_numeric_columns(sample_scalars: dict[str, dict[str, object]]) -> list[str]:
    numeric_columns: list[str] = []
    all_columns = collect_sample_columns(sample_scalars)
    for key in all_columns:
        values = [record.get(key) for record in sample_scalars.values() if key in record]
        if values and all(is_number(value) for value in values):
            numeric_columns.append(key)
    return numeric_columns


def render_project_section(results_dir: Path, pipeline_name: str, project_records: dict[str, dict]) -> str:
    if not project_records:
        return ""
    chunks: list[str] = []
    for record_name, record in project_records.items():
        page_dir = results_dir
        scalars, objects = split_record_values(record)
        figures, links = object_entries(objects, page_dir, results_dir)
        log_href = find_log_link(results_dir, record_name, pipeline_name, page_dir)
        scalar_rows = [
            f"<tr><th>{html.escape(human_title(key))}</th><td>{html.escape(str(value))}</td></tr>"
            for key, value in scalars.items()
        ]
        figure_html = "".join(
            (
                '<figure class="card figure-card">'
                f'<a href="{html.escape(item["path"])}"><img src="{html.escape(item["thumbnail_path"])}" alt="" loading="lazy"></a>'
                f'<figcaption>{html.escape(item["title"])}</figcaption>'
                "</figure>"
            )
            for item in figures
        )
        link_html = "".join(
            f'<li><a href="{html.escape(item["path"])}">{html.escape(item["title"])}</a></li>'
            for item in links
        )
        log_html = f'<p><a href="{html.escape(log_href)}">Project log</a></p>' if log_href else ""
        chunks.append(
            "<section class=\"section\">"
            f"<h2>Project: {html.escape(record_name)}</h2>"
            f"{log_html}"
            + (f'<table class="kv-table">{"".join(scalar_rows)}</table>' if scalar_rows else "")
            + (f'<div class="figure-grid">{figure_html}</div>' if figure_html else "")
            + (
                f'<div class="card"><h3>Project files</h3><ul class="link-list">{link_html}</ul></div>'
                if link_html
                else ""
            )
            + "</section>"
        )
    return "".join(chunks)


def render_sample_page(
    results_dir: Path,
    report_dir: Path,
    sample_name: str,
    pipeline_name: str,
    sample_scalars: dict[str, object],
    sample_objects: dict[str, dict],
    previous_name: str | None,
    next_name: str | None,
    has_sample_object_page: bool,
) -> str:
    page_dir = report_dir
    figures, links = object_entries(sample_objects, page_dir, results_dir)
    log_href = find_log_link(results_dir, sample_name, pipeline_name, page_dir)
    scalar_rows = "".join(
        f"<tr><th>{html.escape(human_title(key))}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in sample_scalars.items()
    )
    nav_links = ['<a href="../report.html">Main report</a>']
    if previous_name:
        nav_links.append(f'<a href="{html.escape(previous_name)}.html">Previous sample</a>')
    if next_name:
        nav_links.append(f'<a href="{html.escape(next_name)}.html">Next sample</a>')
    if has_sample_object_page:
        nav_links.append('<a href="sample_objects.html">Compare sample outputs</a>')
    sample_links = "".join(nav_links)
    figure_html = "".join(
        (
            '<figure class="card figure-card">'
            f'<a class="figure-thumb" href="{html.escape(item["path"])}"><img src="{html.escape(item["thumbnail_path"])}" alt="" loading="lazy"></a>'
            f'<figcaption>{html.escape(item["title"])}</figcaption>'
            "</figure>"
        )
        for item in figures
    )
    link_html = "".join(
        f'<li><a href="{html.escape(item["path"])}">{html.escape(item["title"])}</a></li>'
        for item in links
    )
    log_html = f'<p><a href="{html.escape(log_href)}">Sample log</a></p>' if log_href else ""
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(sample_name)} report</title>
  <link rel="stylesheet" href="report.css">
</head>
<body>
  <header class="topbar">
    <div class="topbar-inner">
      <a class="brand" href="../report.html">{html.escape(pipeline_name)}</a>
      <nav class="nav-inline">{sample_links}</nav>
    </div>
  </header>
  <main class="layout">
    <section class="section">
      <h1>{html.escape(sample_name)}</h1>
      {log_html}
      <table class="kv-table">{scalar_rows}</table>
    </section>
    <section class="section">
      <h2>Sample Figures</h2>
      {f'<div class="figure-grid">{figure_html}</div>' if figure_html else '<p>No sample-specific figures were reported.</p>'}
    </section>
    <section class="section">
      <h2>Sample Files</h2>
      {f'<ul class="link-list">{link_html}</ul>' if link_html else '<p>No sample-specific files were reported.</p>'}
    </section>
  </main>
</body>
</html>
"""


def render_sample_object_summary_page(
    results_dir: Path,
    report_dir: Path,
    pipeline_name: str,
    sample_object_catalog: dict[str, dict],
    object_order: list[str],
) -> str:
    options = "".join(
        f'<option value="{html.escape(key)}">{html.escape(sample_object_catalog[key]["title"])}</option>'
        for key in object_order
    )
    catalog_payload = json.dumps(
        {
            key: {
                "title": sample_object_catalog[key]["title"],
                "items": sample_object_catalog[key]["items"],
            }
            for key in object_order
        }
    ).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(pipeline_name)} sample comparisons</title>
  <link rel="stylesheet" href="report.css">
</head>
<body>
  <header class="topbar">
    <div class="topbar-inner">
      <a class="brand" href="../report.html">{html.escape(pipeline_name)}</a>
      <nav class="nav-inline"><a href="../report.html">Main report</a></nav>
    </div>
  </header>
  <main class="layout">
    <section class="section">
      <h1>Sample-specific output comparison</h1>
      <p class="muted">Compare one sample-level object across all samples.</p>
      <div class="compare-toolbar">
        <label for="sample-object-select">Sample object</label>
        <select id="sample-object-select">
          {options}
        </select>
      </div>
      <div id="sample-object-grid" class="compare-grid"></div>
    </section>
  </main>
  <script>
    const sampleObjectCatalog = {catalog_payload};
    const sampleObjectSelect = document.getElementById("sample-object-select");
    const sampleObjectGrid = document.getElementById("sample-object-grid");

    function escapeHtml(value) {{
      return String(value ?? "")
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;")
        .replace(/'/g, "&#39;");
    }}

    function safeUrl(value) {{
      return escapeHtml(value || "#");
    }}

    function renderSampleObject(objectKey) {{
      const payload = sampleObjectCatalog[objectKey];
      if (!payload) {{
        sampleObjectGrid.innerHTML = "<p>No sample objects available.</p>";
        return;
      }}
      const hasImages = payload.items.some((item) => item.thumbnail_path);
      sampleObjectGrid.className = hasImages ? "compare-grid compare-grid-images" : "compare-grid compare-grid-files";
      sampleObjectGrid.innerHTML = payload.items.map((item) => {{
        const sample = escapeHtml(item.sample);
        const title = escapeHtml(payload.title);
        const path = safeUrl(item.path);
        if (!item.thumbnail_path) {{
          return `
            <article class="compare-file-row">
              <div>
                <div class="figure-sample-name">${{sample}}</div>
                <div class="muted">${{title}}</div>
              </div>
              <a class="file-button" href="${{path}}">Open</a>
            </article>
          `;
        }}
        const thumbnail = safeUrl(item.thumbnail_path);
        return `
          <figure class="card figure-card compare-card">
            <figcaption>
              <span class="figure-sample-name">${{sample}}</span>
              <span class="muted">${{title}}</span>
            </figcaption>
            <a class="figure-thumb compare-thumb" href="${{path}}">
              <img src="${{thumbnail}}" alt="" loading="lazy">
            </a>
          </figure>
        `;
      }}).join("");
    }}

    if (sampleObjectSelect) {{
      sampleObjectSelect.addEventListener("change", () => renderSampleObject(sampleObjectSelect.value));
      renderSampleObject(sampleObjectSelect.value);
    }}
  </script>
</body>
</html>
"""


def render_main_page(
    results_dir: Path,
    report_dir: Path,
    pipeline_name: str,
    sample_scalars: dict[str, dict[str, object]],
    project_records: dict[str, dict],
    statuses: dict[str, str],
    has_sample_object_page: bool,
) -> str:
    sample_names = list(sample_scalars.keys())
    sample_count = len(sample_names)
    columns = [key for key in collect_sample_columns(sample_scalars) if key != "meta"]
    numeric_columns = collect_numeric_columns(sample_scalars)
    headers = "".join(f"<th>{html.escape(human_title(column))}</th>" for column in columns)
    rows = []
    for sample_name in sample_names:
        row = [
            f'<td><a href="report/{html.escape(sample_name)}.html">{html.escape(sample_name)}</a></td>',
            f"<td>{html.escape(statuses.get(sample_name, 'completed' if sample_scalars[sample_name].get('Success') else 'unknown'))}</td>",
        ]
        log_href = find_log_link(results_dir, sample_name, pipeline_name, results_dir)
        log_cell = f'<a href="{html.escape(log_href)}">log</a>' if log_href else ""
        row.append(f"<td>{log_cell}</td>")
        for column in columns:
            value = sample_scalars[sample_name].get(column, "")
            row.append(f"<td>{html.escape(str(value))}</td>")
        rows.append("<tr>" + "".join(row) + "</tr>")

    plot_data = {
        sample_name: {
            column: numeric_value(values[column])
            for column in numeric_columns
            if column in values and is_number(values[column])
        }
        for sample_name, values in sample_scalars.items()
    }
    project_html = render_project_section(results_dir, pipeline_name, project_records)
    sample_objects_link = (
        '<a href="report/sample_objects.html">Compare sample outputs</a>'
        if has_sample_object_page
        else ""
    )
    metric_options = "".join(
        (
            f'<button type="button" class="metric-item{" active" if index == 0 else ""}" '
            f'data-metric="{html.escape(column)}">{html.escape(human_title(column))}</button>'
        )
        for index, column in enumerate(numeric_columns)
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(pipeline_name)} report</title>
  <link rel="stylesheet" href="report/report.css">
</head>
<body>
  <header class="topbar">
    <div class="topbar-inner">
      <span class="brand">{html.escape(pipeline_name)}</span>
      <nav class="nav-inline">{sample_objects_link}</nav>
    </div>
  </header>
  <main class="layout">
    <section class="section">
      <h1>Results summary</h1>
      <p class="muted">Generated directly from per-record <code>stats.yaml</code> files in the results directory.</p>
      <div class="table-wrap">
        <table class="summary-table">
          <thead>
            <tr>
              <th>Sample</th>
              <th>Status</th>
              <th>Log</th>
              {headers}
            </tr>
          </thead>
          <tbody>
            {''.join(rows)}
          </tbody>
        </table>
      </div>
    </section>
    <section class="section">
      <h2>Numeric sample plots</h2>
      <p class="muted">Only values reported as numeric across samples are offered here.</p>
      <div class="plot-panel">
        <aside class="metric-sidebar">
          <div class="metric-sidebar-title">Metrics</div>
          <div id="metric-list" class="metric-list">
            {metric_options}
          </div>
        </aside>
        <div class="plot-area">
          <div class="plot-caption">
            <span id="plot-title" class="plot-title"></span>
            <span id="plot-meta" class="plot-meta"></span>
          </div>
          <div class="plot-scroller">
            <div id="plot-root" class="plot-root"></div>
          </div>
        </div>
      </div>
    </section>
    {project_html}
  </main>
  <script>
    const plotData = {json.dumps(plot_data)};
    const metricList = document.getElementById("metric-list");
    const plotRoot = document.getElementById("plot-root");
    const plotTitle = document.getElementById("plot-title");
    const plotMeta = document.getElementById("plot-meta");
    const metricButtons = Array.from(document.querySelectorAll(".metric-item"));

    function autoBarWidth(sampleCount) {{
      if (sampleCount <= 3) return 72;
      if (sampleCount <= 6) return 56;
      if (sampleCount <= 12) return 38;
      if (sampleCount <= 20) return 28;
      return 22;
    }}

    function displayLabel(sample) {{
      return sample.length > 18 ? `${{sample.slice(0, 17)}}...` : sample;
    }}

    function renderBars(metric) {{
      if (!metric) {{
        plotRoot.innerHTML = "<p>No numeric sample metrics were reported.</p>";
        plotTitle.textContent = "";
        plotMeta.textContent = "";
        return;
      }}
      const rows = Object.entries(plotData)
        .filter(([, values]) => values[metric] !== undefined)
        .map(([sample, values]) => [sample, values[metric]]);
      if (!rows.length) {{
        plotRoot.innerHTML = "<p>No numeric values available for this metric.</p>";
        plotTitle.textContent = metric.replaceAll("_", " ");
        plotMeta.textContent = "";
        return;
      }}
      const barWidth = autoBarWidth(rows.length);
      const maxValue = Math.max(...rows.map(([, value]) => value), 1);
      const chartHeight = 260;
      const chartTop = 24;
      const axisLeft = 78;
      const labelZone = 96;
      const chartBottom = chartTop + chartHeight;
      const innerGap = rows.length <= 3 ? 28 : 14;
      const minPlotWidth = 560;
      const barsSpan = rows.length * barWidth + Math.max(0, rows.length - 1) * innerGap;
      const plotWidth = Math.max(minPlotWidth, barsSpan + 60);
      const svgWidth = axisLeft + plotWidth + 24;
      const svgHeight = chartBottom + labelZone;
      const startX = axisLeft + Math.max(20, (plotWidth - barsSpan) / 2);
      const gridTicks = [0, 0.25, 0.5, 0.75, 1];

      const gridHtml = gridTicks.map((fraction) => {{
        const value = maxValue * (1 - fraction);
        const y = chartTop + chartHeight * fraction;
        const label = Number.isInteger(value) ? value : value.toFixed(2);
        return `
          <line x1="${{axisLeft}}" y1="${{y}}" x2="${{axisLeft + plotWidth}}" y2="${{y}}" class="svg-grid-line" />
          <text x="${{axisLeft - 10}}" y="${{y + 5}}" text-anchor="end" class="svg-axis-text">${{label}}</text>
        `;
      }}).join("");

      const barsHtml = rows.map(([sample, value], index) => {{
        const x = startX + index * (barWidth + innerGap);
        const height = Math.max(2, (value / maxValue) * chartHeight);
        const y = chartBottom - height;
        const xCenter = x + barWidth / 2;
        const shownLabel = displayLabel(sample);
        const shownValue = Number.isInteger(value) ? value : value.toFixed(2);
        return `
          <g>
            <title>${{sample}}: ${{shownValue}}</title>
            <text x="${{xCenter}}" y="${{y - 8}}" text-anchor="middle" class="svg-value-text">${{shownValue}}</text>
            <rect x="${{x}}" y="${{y}}" width="${{barWidth}}" height="${{height}}" rx="8" ry="8" class="svg-bar" />
            <text x="${{xCenter}}" y="${{chartBottom + 14}}" text-anchor="end" transform="rotate(-45 ${{xCenter}} ${{chartBottom + 14}})" class="svg-label-text">${{shownLabel}}</text>
          </g>
        `;
      }}).join("");

      plotRoot.innerHTML = `
        <svg class="plot-svg" viewBox="0 0 ${{svgWidth}} ${{svgHeight}}" width="${{svgWidth}}" height="${{svgHeight}}" aria-label="${{metric}}">
          <defs>
            <linearGradient id="barGradient" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%" stop-color="#67a7f2" />
              <stop offset="100%" stop-color="#2f76d2" />
            </linearGradient>
          </defs>
          <line x1="${{axisLeft}}" y1="${{chartTop}}" x2="${{axisLeft}}" y2="${{chartBottom}}" class="svg-axis-line" />
          <line x1="${{axisLeft}}" y1="${{chartBottom}}" x2="${{axisLeft + plotWidth}}" y2="${{chartBottom}}" class="svg-axis-line" />
          ${{gridHtml}}
          ${{barsHtml}}
        </svg>
      `;
      plotTitle.textContent = metric.replaceAll("_", " ");
      plotMeta.textContent = `${{rows.length}} sample${{rows.length === 1 ? "" : "s"}}`;
    }}

    if (metricButtons.length) {{
      metricButtons.forEach((button) => {{
        button.addEventListener("click", () => {{
          metricButtons.forEach((item) => item.classList.remove("active"));
          button.classList.add("active");
          renderBars(button.dataset.metric);
        }});
      }});
      renderBars(metricButtons[0].dataset.metric);
    }} else {{
      renderBars(null);
    }}
  </script>
</body>
</html>
"""


def stylesheet() -> str:
    return """body {
  margin: 0;
  font-family: "Segoe UI", Arial, sans-serif;
  background: #f6f8fb;
  color: #1f2937;
}

.topbar {
  background: #18222f;
  color: #fff;
  padding: 16px 24px;
}

.muted-light {
  color: #d7e3f4;
}

.topbar-inner {
  display: flex;
  gap: 24px;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
}

.brand {
  color: #fff;
  text-decoration: none;
  font-size: 1.2rem;
  font-weight: 700;
}

.nav-inline {
  display: flex;
  gap: 14px;
  flex-wrap: wrap;
}

.nav-inline a {
  color: #d7e3f4;
  text-decoration: none;
}

.layout {
  max-width: 1400px;
  margin: 0 auto;
  padding: 24px;
}

.section {
  background: #fff;
  border-radius: 14px;
  padding: 24px;
  margin-bottom: 24px;
  box-shadow: 0 8px 24px rgba(15, 23, 42, 0.08);
}

.muted {
  color: #5b6575;
}

.table-wrap {
  overflow-x: auto;
}

.summary-table, .kv-table {
  width: 100%;
  border-collapse: collapse;
}

.summary-table th, .summary-table td, .kv-table th, .kv-table td {
  border-bottom: 1px solid #e5e7eb;
  padding: 10px 12px;
  text-align: left;
  vertical-align: top;
}

.summary-table th, .kv-table th {
  background: #f8fafc;
}

.plot-panel {
  display: grid;
  grid-template-columns: 220px minmax(0, 1fr);
  gap: 24px;
  align-items: start;
}

.metric-sidebar {
  border: 1px solid #e5e7eb;
  border-radius: 12px;
  background: #fbfcfe;
  overflow: hidden;
}

.metric-sidebar-title {
  padding: 14px 16px;
  font-weight: 700;
  border-bottom: 1px solid #e5e7eb;
  background: #f6f8fc;
}

.metric-list {
  max-height: 320px;
  overflow-y: auto;
}

.metric-item {
  display: block;
  width: 100%;
  text-align: left;
  padding: 12px 16px;
  border: 0;
  border-bottom: 1px solid #eef2f7;
  background: transparent;
  color: #334155;
  cursor: pointer;
  font: inherit;
}

.metric-item:hover {
  background: #eef5ff;
}

.metric-item.active {
  background: linear-gradient(90deg, #e8f1ff, #f8fbff);
  color: #1d4ed8;
  font-weight: 700;
}

.plot-area {
  min-width: 0;
}

.plot-caption {
  display: flex;
  justify-content: space-between;
  gap: 16px;
  align-items: baseline;
  margin-bottom: 12px;
  flex-wrap: wrap;
}

.plot-title {
  font-size: 1.05rem;
  font-weight: 700;
}

.plot-meta {
  color: #64748b;
}

.plot-scroller {
  overflow-x: auto;
  padding-bottom: 8px;
}

.plot-root {
  min-height: 400px;
  min-width: max-content;
}

.plot-svg {
  display: block;
}

.svg-axis-line {
  stroke: #cfd8e3;
  stroke-width: 1;
}

.svg-grid-line {
  stroke: #e5ebf2;
  stroke-width: 1;
  stroke-dasharray: 4 4;
}

.svg-axis-text {
  fill: #64748b;
  font-size: 13px;
}

.svg-value-text {
  fill: #475569;
  font-size: 13px;
  font-weight: 600;
}

.svg-label-text {
  fill: #1f2937;
  font-size: 13px;
}

.svg-bar {
  fill: url(#barGradient);
}

.figure-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 18px;
}

.compare-grid {
  margin-top: 18px;
}

.compare-grid-images {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 16px;
}

.compare-grid-files {
  display: grid;
  gap: 10px;
}

.compare-card {
  display: flex;
  flex-direction: column;
  min-width: 0;
}

.compare-card figcaption {
  display: flex;
  flex-direction: column;
  gap: 4px;
  margin: 0 0 10px;
}

.figure-sample-name {
  font-size: 1rem;
  font-weight: 700;
  margin-bottom: 12px;
}

.card {
  border: 1px solid #e5e7eb;
  border-radius: 12px;
  padding: 16px;
  background: #fcfdff;
}

.figure-thumb {
  display: flex;
  align-items: center;
  justify-content: center;
  height: 260px;
  padding: 12px;
  border-radius: 10px;
  background: #ffffff;
  border: 1px solid #e8edf3;
  overflow: hidden;
  text-decoration: none;
}

.figure-thumb img {
  max-width: 100%;
  max-height: 100%;
  width: auto;
  height: auto;
  display: block;
  object-fit: contain;
}

.figure-thumb-empty {
  color: #1d4ed8;
  font-weight: 600;
}

.compare-thumb {
  height: 220px;
}

.compare-file-row {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 16px;
  min-width: 0;
  padding: 12px 14px;
  border: 1px solid #e5e7eb;
  border-radius: 8px;
  background: #fcfdff;
}

.compare-file-row .figure-sample-name {
  margin-bottom: 2px;
}

.file-button {
  flex: 0 0 auto;
  padding: 7px 11px;
  border-radius: 8px;
  border: 1px solid #bfdbfe;
  background: #eff6ff;
  color: #1d4ed8;
  font-weight: 700;
  text-decoration: none;
}

.figure-card figcaption {
  margin-top: 10px;
  font-weight: 600;
}

.figure-card.compare-card figcaption {
  margin: 0 0 10px;
}

.compare-card .figure-sample-name {
  margin-bottom: 0;
}

.link-list {
  margin: 0;
  padding-left: 20px;
}

@media (max-width: 900px) {
  .plot-panel {
    grid-template-columns: 1fr;
  }
}
"""


def clean_legacy_report_layout(results_dir: Path) -> None:
    for legacy_dir in ("report_assets", "report_samples"):
        path = results_dir / legacy_dir
        if path.exists():
            shutil.rmtree(path)
    legacy_page = results_dir / "sample_objects.html"
    if legacy_page.exists():
        legacy_page.unlink()


def main() -> int:
    args = parse_args()
    looper_dir = Path.cwd()
    results_dir = resolve_results_dir(looper_dir, args)
    pipeline_name = infer_pipeline_name(results_dir, args.pipeline_name)
    statuses = discover_statuses(results_dir / "flags", pipeline_name)
    sample_records, project_records = collect_records(results_dir, pipeline_name)

    sample_scalars = {
        name: split_record_values(record)[0]
        for name, record in sorted(sample_records.items())
    }
    sample_objects = {
        name: split_record_values(record)[1]
        for name, record in sorted(sample_records.items())
    }
    project_records = dict(sorted(project_records.items()))
    report_dir = results_dir / "report"
    clean_legacy_report_layout(results_dir)
    sample_object_catalog, sample_object_order = collect_sample_object_catalog(
        sample_objects, report_dir, results_dir
    )

    report_dir.mkdir(exist_ok=True)
    (report_dir / "report.css").write_text(stylesheet(), encoding="utf-8")

    sample_names = list(sample_scalars.keys())
    for index, sample_name in enumerate(sample_names):
        previous_name = sample_names[index - 1] if index > 0 else None
        next_name = sample_names[index + 1] if index + 1 < len(sample_names) else None
        page = render_sample_page(
            results_dir,
            report_dir,
            sample_name,
            pipeline_name,
            sample_scalars[sample_name],
            sample_objects[sample_name],
            previous_name,
            next_name,
            bool(sample_object_order),
        )
        (report_dir / f"{sample_name}.html").write_text(page, encoding="utf-8")

    if sample_object_order:
        sample_objects_page = render_sample_object_summary_page(
            results_dir,
            report_dir,
            pipeline_name,
            sample_object_catalog,
            sample_object_order,
        )
        (report_dir / "sample_objects.html").write_text(sample_objects_page, encoding="utf-8")

    main_page = render_main_page(
        results_dir,
        report_dir,
        pipeline_name,
        sample_scalars,
        project_records,
        statuses,
        bool(sample_object_order),
    )
    (results_dir / "report.html").write_text(main_page, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
