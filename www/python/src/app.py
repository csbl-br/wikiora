from flask import Flask, request, render_template, send_from_directory
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
import json
import random

app = Flask(__name__, static_url_path="/static")


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


@app.route("/enrich", methods=["GET"])
def enrich():
    organism = request.args.get("organism")
    analysis_type = request.args.get("analysis_type")
    gene_list = request.args.get("gene_list")
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

    for go_term, details in gene_sets.items():
        overlap = set(genes).intersection(details["genes"])
        if overlap:
            x = len(overlap)
            n = len(details["genes"])
            p_value = hypergeometric_test(M, n, N, x)
            results.append(
                {
                    "GO Term": go_term,
                    "Description": details["description"],
                    "Wikipedia URL": details["wikipedia_url"],
                    "Overlap": ", ".join(overlap),
                    "Count": x,
                    "p-value": p_value,
                }
            )

    if results:
        results_df = pd.DataFrame(results)
        results_df["corrected p-value"] = multipletests(
            results_df["p-value"], method="bonferroni"
        )[1]
        results_df = results_df.sort_values(by="p-value").head(
            10
        )  # Show only the top 10 enriched sets
        plot_results(results_df)
        results = results_df.to_dict(orient="records")
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
    df["-logP"] = -np.log10(df["corrected p-value"])
    sns.barplot(x="-logP", y="Description", data=df, palette="viridis")
    plt.axvline(-np.log10(0.05), color="red", linestyle="--", linewidth=1)
    plt.xlabel("-log(corrected p-value)")
    plt.ylabel("GO Term Description")
    plt.title("Top 10 Enriched GO Terms by Corrected p-value")
    plt.tight_layout()
    plt.savefig("static/enrichment_plot.png")


if __name__ == "__main__":
    app.run(debug=True)
