from flask import Flask, request, render_template, jsonify, send_from_directory
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
import json
import random
import os
import matplotlib.colors as mcolors  # Add this import
import sqlite3

app = Flask(__name__, static_url_path="/static")
__version__ = "0.3.0"

DATABASE = "database.db"


def init_db():
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS counter (
                id INTEGER PRIMARY KEY,
                count INTEGER NOT NULL
            )
        """
        )
        cursor.execute(
            """
            INSERT INTO counter (count) 
            SELECT 0 WHERE NOT EXISTS (SELECT 1 FROM counter)
        """
        )
        conn.commit()
        # Debugging log
        cursor.execute("SELECT count FROM counter")
        count = cursor.fetchone()[0]
        print(f"Database initialized with counter value: {count}")


def increment_counter():
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("UPDATE counter SET count = count + 1")
        conn.commit()


def get_counter():
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT count FROM counter")
        count = cursor.fetchone()[0]
        return count


@app.route("/api/lists_enriched", methods=["GET"])
def get_lists_enriched():
    lists_count = get_counter()
    print(lists_count)
    return jsonify({"lists_enriched": lists_count})


def get_version():
    with open("static/version.txt", "r") as f:
        version = f.read().strip()
    return version


@app.context_processor
def inject_version():
    return {"app_version": get_version()}


# Route to serve robots.txt
@app.route("/robots.txt")
def robots_txt():
    return send_from_directory(app.static_folder, "robots.txt")


# Load GMT file into a dictionary
def load_gmt(file_path):
    gene_sets = {}
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            go_term = parts[0]
            description = parts[1]
            wikipedia_url = parts[2]
            genes = parts[3:]
            gene_sets[go_term] = {
                "description": description,
                "wikipedia_url": wikipedia_url,
                "genes": set(genes),
            }
    return gene_sets


# Load genes.json
with open(os.path.join("static", "genes.json")) as f:
    genes_data = json.load(f)


# Load process information from JSON file
def load_processes(file_path):
    with open(file_path, "r") as f:
        return json.load(f)


# Function to select a random process and get up to 50 genes
def select_random_process(processes):
    process = random.choice(processes)
    genes = process["gene_symbol"][:50]
    return genes


def select_random_separator():
    separators = [", ", "\t", "\n", "; "]
    probs = [0.8, 0.1, 0.3, 0.1]
    separator = random.choices(separators, probs)[0]
    return separator


# Function to perform hypergeometric test
def hypergeometric_test(M, n, N, x):
    return hypergeom.sf(x - 1, M, n, N)


# Function to parse gene list
def parse_gene_list(gene_list):
    genes = re.split(r"[\s,;]+", gene_list.strip())
    return [gene for gene in genes if gene]


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/example_genes", methods=["POST"])
def example_genes():
    organism = request.form["organism"]
    analysis_type = request.form["analysis_type"]
    if analysis_type == "biological_processes":
        if organism == "human":
            processes = load_processes(
                "static/processes_human_biological_processes.json"
            )
        else:
            processes = load_processes(
                "static/processes_mouse_biological_processes.json"
            )
    elif analysis_type == "molecular_functions":
        if organism == "human":
            processes = load_processes(
                "static/processes_human_molecular_functions.json"
            )
        else:
            processes = load_processes(
                "static/processes_mouse_molecular_functions.json"
            )
    elif analysis_type == "cellular_components":
        if organism == "human":
            processes = load_processes(
                "static/processes_human_cellular_components.json"
            )
        else:
            processes = load_processes(
                "static/processes_mouse_cellular_components.json"
            )
    else:  # cell_type_markers is default
        if organism == "human":
            processes = load_processes("static/processes_human_cell_type.json")
        else:
            processes = load_processes("static/processes_mouse_cell_type.json")

    default_genes = select_random_process(processes)
    separator = select_random_separator()
    default_genes_str = separator.join(default_genes)
    return default_genes_str


@app.route("/api/enrich", methods=["GET", "POST"])
def api_enrich():
    organism = (
        request.args.get("organism")
        if request.method == "GET"
        else request.form.get("organism")
    )
    analysis_type = (
        request.args.get("analysis_type")
        if request.method == "GET"
        else request.form.get("analysis_type")
    )
    gene_list = (
        request.args.get("gene_list")
        if request.method == "GET"
        else request.form.get("gene_list")
    )
    genes = parse_gene_list(gene_list)
    results = []

    if analysis_type == "biological_processes":
        if organism == "human":
            gene_sets = load_gmt("static/gene_sets_human_biological_processes.gmt")
        else:
            gene_sets = load_gmt("static/gene_sets_mouse_biological_processes.gmt")
    elif analysis_type == "molecular_functions":
        if organism == "human":
            gene_sets = load_gmt("static/gene_sets_human_molecular_functions.gmt")
        else:
            gene_sets = load_gmt("static/gene_sets_mouse_molecular_functions.gmt")
    elif analysis_type == "cellular_components":
        if organism == "human":
            gene_sets = load_gmt("static/gene_sets_human_cellular_components.gmt")
        else:
            gene_sets = load_gmt("static/gene_sets_mouse_cellular_components.gmt")
    else:  # cell_type_markers is default
        if organism == "human":
            gene_sets = load_gmt("static/gene_sets_human_cell_type.gmt")
        else:
            gene_sets = load_gmt("static/gene_sets_mouse_cell_type.gmt")

    # Total number of genes in the background
    unique_genes = set()
    for details in gene_sets.values():
        unique_genes.update(details["genes"])
    M = len(unique_genes)
    N = len(genes)
    total_gene_sets = len(gene_sets)

    all_results = []
    for term, details in gene_sets.items():
        overlap = set(genes).intersection(details["genes"])
        x = len(overlap)
        n = len(details["genes"])
        p_value = hypergeometric_test(M, n, N, x)
        odds_ratio = (1.0 * x * (M - n - N + x)) / max(1.0 * (n - x) * (N - x), 1)
        combined_score = -np.log10(p_value) * odds_ratio
        overlap_info = []
        for gene in overlap:
            gene_info = genes_data.get(gene, {})
            gene_link = gene_info.get(
                "wikipediaLink", f"https://en.wikipedia.org/wiki/{gene}"
            )
            gene_status = gene_info.get("pageStatus", "red")
            overlap_info.append(
                {"gene": gene, "link": gene_link, "status": gene_status}
            )
        all_results.append(
            {
                "Term": term,
                "Description": details["description"],
                "Wikipedia URL": details["wikipedia_url"],
                "Overlap": overlap_info,
                "Count": x,
                "Gene Set Size": n,  # Include Gene Set Size in the results
                "p-value": p_value,
                "Odds Ratio": odds_ratio,
                "Combined Score": combined_score,
                "Gene Ratio": x / n,  # Add Gene Ratio metric
            }
        )

    if all_results:
        results_df = pd.DataFrame(all_results)
        # Apply Benjamini-Hochberg correction considering the total number of gene sets
        results_df["q-value"] = multipletests(results_df["p-value"], method="fdr_bh")[1]
        # Filter to include only those sets with overlap >= 1
        results_df = results_df[results_df["Count"] >= 1]
        results_df = results_df.sort_values(by="p-value").head(
            10
        )  # Show only the top 10 enriched sets with lowest p-values
        results = results_df.to_dict(orient="records")

    return jsonify(results)


@app.route("/enrich", methods=["GET", "POST"])
def enrich():
    organism = (
        request.args.get("organism")
        if request.method == "GET"
        else request.form.get("organism")
    )
    analysis_type = (
        request.args.get("analysis_type")
        if request.method == "GET"
        else request.form.get("analysis_type")
    )
    gene_list = (
        request.args.get("gene_list")
        if request.method == "GET"
        else request.form.get("gene_list")
    )
    response = api_enrich()
    results = response.get_json()
    increment_counter()
    return render_template("results.html", results=results)


@app.route("/download")
def download():
    return render_template("download.html")


@app.route("/about")
def about():
    return render_template("about.html")


# Not currently in use. May be used in the future as an API for download.
def plot_results(df):
    # Calculate -logP for the barplot
    df["-logP"] = -np.log10(df["q-value"])

    # Cap the Count values at 30 for dot sizes
    df["Count"] = np.minimum(df["Count"], 30)

    # Sort by Gene Ratio for dotplot
    df_sorted_by_gene_ratio = df.sort_values(by="Gene Ratio", ascending=False)
    plt.figure(figsize=(12, 12))

    cmap = plt.cm.Blues_r

    # Plot barplot
    plt.subplot(2, 1, 1)
    colors_barplot = cmap(df["q-value"])
    barplot = sns.barplot(x="-logP", y="Description", data=df, palette=colors_barplot)
    plt.axvline(-np.log10(0.05), color="red", linestyle="--", linewidth=1)
    plt.xlabel("-log(q-value)", fontsize=14)
    plt.ylabel("")
    plt.title("Top 10 Enriched Terms by q-value", fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Adding colorbar legend to barplot
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.0, vmax=1.0))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=barplot, orientation="vertical", pad=0.01)
    cbar.set_label("q-value", fontsize=12)

    # Plot dotplot
    plt.subplot(2, 1, 2)
    colors_dotplot = cmap(df_sorted_by_gene_ratio["q-value"])
    dotplot = plt.scatter(
        df_sorted_by_gene_ratio["Gene Ratio"],
        df_sorted_by_gene_ratio["Description"],
        s=df_sorted_by_gene_ratio["Count"] ** 0.7 * 23,
        c=colors_dotplot,
        edgecolors="black",  # Set the border color to black
        linewidth=0.8,
    )

    cbar = plt.colorbar(sm, ax=plt.gca(), orientation="vertical", pad=0.01)
    cbar.set_label("q-value", fontsize=12)

    plt.xlabel("Gene Ratio (overlap/set length)", fontsize=14)
    plt.ylabel("")
    plt.title("Top 10 Enriched Terms by Gene Ratio", fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.gca().invert_yaxis()  # Invert the y-axis

    # Adding dot size legend
    example_sizes = [2, 10, 30]
    for size in example_sizes:
        plt.scatter(
            [],
            [],
            c="k",
            alpha=0.5,
            s=size**0.7 * 23,
            label=f'{size}{"+" if size == 30 else ""}',
        )

    plt.legend(
        scatterpoints=1,
        frameon=True,
        labelspacing=1,
        title="Count",
        loc="lower right",
        fontsize=12,
    )

    plt.tight_layout()
    plt.savefig("static/enrichment_plot.png", dpi=300)


if __name__ == "__main__":
    init_db()
    debug_mode = os.environ.get("FLASK_DEBUG", "False").lower() in ["true", "1", "t"]
    app.run()
