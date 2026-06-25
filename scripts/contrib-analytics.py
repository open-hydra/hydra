#!/usr/bin/env python3

import csv
import json
import os
import subprocess
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd

COMMIT_WEIGHT = 0.0
LINE_WEIGHT = 1.0

IGNORED_AUTHORS = {
    "github-actions[bot] <41898282+github-actions[bot]@users.noreply.github.com>",
}

IGNORED_EXTENSIONS = {
    ".tec",
    ".vts",
    ".dat",
}

# Repo names (as produced by repo_name()) whose initial commit is excluded.
SKIP_INITIAL_COMMITS = {
    "lib_MOSKA",
}


def run(cmd, cwd=None):
    return subprocess.check_output(
        cmd,
        cwd=cwd,
        text=True,
        stderr=subprocess.STDOUT,
    )


def discover_submodules(root):
    """
    Returns:
        initialized_repos
        skipped_repos
    """

    repos = [root]
    skipped = []

    gitmodules = os.path.join(root, ".gitmodules")

    if not os.path.exists(gitmodules):
        return repos, skipped

    try:
        out = run(
            [
                "git",
                "config",
                "--file",
                ".gitmodules",
                "--get-regexp",
                "path",
            ],
            cwd=root,
        )
    except Exception:
        return repos, skipped

    paths = []

    for line in out.splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) == 2:
            paths.append(parts[1].strip())

    for relpath in paths:

        fullpath = os.path.join(root, relpath)

        if not os.path.exists(fullpath):
            skipped.append(
                (relpath, "path missing")
            )
            continue

        if not os.path.exists(
            os.path.join(fullpath, ".git")
        ):
            skipped.append(
                (relpath, "not initialized")
            )
            continue

        repos.append(fullpath)

    return repos, skipped


def parse_repository(repo_path, mailmap=None, skip_initial=False):

    stats = defaultdict(
        lambda: {
            "commits": 0,
            "additions": 0,
            "deletions": 0,
        }
    )

    cmd = [
        "git",
        "log",
        "--all",
        "--numstat",
        "--format=@@@%aN <%aE>",
    ]

    if mailmap and os.path.isfile(mailmap):
        cmd = ["git", "-c", f"mailmap.file={mailmap}"] + cmd[1:]

    try:
        run(["git", "fetch", "--all", "--quiet"], cwd=repo_path)
    except Exception:
        pass

    if skip_initial:
        try:
            roots = run(
                ["git", "rev-list", "--max-parents=0", "--all"],
                cwd=repo_path,
            ).split()
            cmd.extend(f"^{r}" for r in roots)
        except Exception:
            pass

    output = run(cmd, cwd=repo_path)

    current_author = None

    for line in output.splitlines():

        if line.startswith("@@@"):
            current_author = line[3:].strip()
            if current_author not in IGNORED_AUTHORS:
                stats[current_author]["commits"] += 1
            continue

        parts = line.split()

        if len(parts) != 3:
            continue

        if parts[0] == "-" or parts[1] == "-":
            continue

        try:
            add = int(parts[0])
            delete = int(parts[1])
        except ValueError:
            continue

        if current_author in IGNORED_AUTHORS:
            continue

        filepath = parts[2].split("=>")[-1].strip().rstrip("}")
        _, ext = os.path.splitext(filepath.split("/")[-1])
        if ext.lower() in IGNORED_EXTENSIONS:
            continue

        stats[current_author]["additions"] += add
        stats[current_author]["deletions"] += delete

    rows = []

    for author, s in stats.items():

        changed = s["additions"] + s["deletions"]

        rows.append(
            {
                "author": author,
                "commits": s["commits"],
                "additions": s["additions"],
                "deletions": s["deletions"],
                "changed_lines": changed,
            }
        )

    return pd.DataFrame(rows)


def score_dataframe(df):

    total_commits = df["commits"].sum()
    total_additions = df["additions"].sum()
    total_deletions = df["deletions"].sum()

    if total_commits == 0 or total_additions == 0 or total_deletions == 0:
        df["score"] = 0.0
        return df

    df["commit_share"] = (
        df["commits"] / total_commits
    )

    df["line_share"] = (
        0.9 * df["additions"] / total_additions + 0.1 * df["deletions"] / total_deletions
    )

    df["score"] = (
        COMMIT_WEIGHT * df["commit_share"]
        + LINE_WEIGHT * df["line_share"]
    )

    return df.sort_values(
        "score",
        ascending=False,
    )


def _build_palette():
    colors = []
    for name in ("tab20", "tab20b", "tab20c"):
        colors.extend(plt.cm.get_cmap(name).colors)
    return colors


_PALETTE = _build_palette()


def save_pie_chart(df, title, filename, color_map=None):

    if df.empty:
        return

    chart = df.copy()
    total = chart["score"].sum()
    chart["pct"] = chart["score"] / total * 100

    if color_map:
        colors = [color_map.get(a, "#aaaaaa") for a in chart["author"]]
    else:
        colors = _PALETTE[:len(chart)]

    height = max(6, len(chart) * 0.45)
    fig, ax = plt.subplots(figsize=(11, height))

    wedges, _, autotexts = ax.pie(
        chart["score"],
        autopct=lambda p: f"{p:.1f}%" if p >= 3 else "",
        colors=colors,
        startangle=90,
        wedgeprops={
            "linewidth": 1.0,
            "edgecolor": "white",
            "width": 0.6,
        },
        pctdistance=0.75,
    )

    for t in autotexts:
        t.set_fontsize(9)
        t.set_fontweight("bold")
        t.set_color("white")

    labels = [
        f"{row.author}  ({row.pct:.1f}%)"
        for row in chart.itertuples()
    ]

    ax.legend(
        wedges,
        labels,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=8,
        frameon=False,
    )

    ax.set_title(title, fontsize=13, fontweight="bold", pad=16)

    fig.tight_layout()
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close(fig)


def repo_name(root, repo):

    if repo == root:
        return "main"

    return os.path.relpath(repo, root).replace("/", "_")


def main():

    root = os.getcwd()

    mailmap_path = os.path.join(root, ".mailmap")

    repos, skipped = discover_submodules(root)

    analyzed = []
    aggregate = []
    repo_results = []

    for repo in repos:

        name = repo_name(root, repo)

        try:

            print(f"Processing {name}")

            df = parse_repository(
                repo,
                mailmap=mailmap_path,
                skip_initial=(name in SKIP_INITIAL_COMMITS),
            )

            if df.empty:

                skipped.append(
                    (name, "no commits")
                )

                continue

            df = score_dataframe(df)

            df.to_csv(
                f"{name}_contributors.csv",
                index=False,
            )

            analyzed.append(name)
            repo_results.append((name, df))

            aggregate.append(
                df[
                    [
                        "author",
                        "commits",
                        "additions",
                        "deletions",
                        "changed_lines",
                    ]
                ].copy()
            )

        except Exception as e:

            skipped.append(
                (name, str(e))
            )

    if not aggregate:
        print("No repositories analyzed")
        return

    overall = (
        pd.concat(aggregate)
        .groupby("author")
        .sum()
        .reset_index()
    )

    overall = score_dataframe(overall)

    overall.to_csv(
        "contributors.csv",
        index=False,
    )

    overall.to_json(
        "contributors.json",
        orient="records",
        indent=2,
    )

    all_authors = set(overall["author"])
    for _, df in repo_results:
        all_authors.update(df["author"])

    color_map = {
        author: _PALETTE[i % len(_PALETTE)]
        for i, author in enumerate(sorted(all_authors))
    }

    for name, df in repo_results:
        save_pie_chart(
            df,
            f"Contributions - {name}",
            f"{name}_pie.png",
            color_map=color_map,
        )

    save_pie_chart(
        overall,
        "Hydra Codes Contributions",
        "overall_pie.png",
        color_map=color_map,
    )

    with open(
        "not_analyzed_modules.csv",
        "w",
        newline="",
    ) as f:

        writer = csv.writer(f)

        writer.writerow(
            ["module", "reason"]
        )

        writer.writerows(skipped)

    with open(
        "analysis_summary.txt",
        "w",
    ) as f:

        f.write(
            "Analyzed modules:\n"
        )

        for x in analyzed:
            f.write(f"  {x}\n")

        f.write("\n")

        f.write(
            "Not analyzed:\n"
        )

        for name, reason in skipped:
            f.write(
                f"  {name}: {reason}\n"
            )

    print()
    print(
        f"Analyzed: {len(analyzed)}"
    )
    print(
        f"Skipped: {len(skipped)}"
    )
    print(
        "Generated contributors.csv, contributors.json, overall_pie.png"
    )


if __name__ == "__main__":
    main()