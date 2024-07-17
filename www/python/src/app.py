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

app = Flask(__name__, static_url_path="/static")
__version__ = "0.1.2"


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
    default_genes_str = ", ".join(default_genes)
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
    M = sum(len(details["genes"]) for details in gene_sets.values())
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
                "p-value": p_value,
                "Odds Ratio": odds_ratio,
                "Combined Score": combined_score,
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
        plot_results(results_df)
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
    return render_template("results.html", results=results)


@app.route("/download")
def download():
    return render_template("download.html")


@app.route("/about")
def about():
    return render_template("about.html")


# Function to plot the results
def plot_results(df):
    plt.figure(figsize=(10, 8))
    df["-logP"] = -np.log10(df["q-value"])
    sns.barplot(x="-logP", y="Description", data=df, palette="viridis")
    plt.axvline(-np.log10(0.05), color="red", linestyle="--", linewidth=1)
    plt.xlabel("-log(corrected q-value)")
    plt.ylabel("Gene Set Description")
    plt.title("Top 10 Enriched Terms by q-value (FDR 0.05)")
    plt.tight_layout()
    plt.savefig("static/enrichment_plot.png")


if __name__ == "__main__":
    app.run(debug=True)
