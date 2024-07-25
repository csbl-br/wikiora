import os
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON
from jinja2 import Template
from pathlib import Path

HERE = Path(__file__).parent.resolve()

STATIC = HERE.joinpath("static").resolve()


# Template SPARQL query for GO terms using Jinja2
go_query_template = """
SELECT DISTINCT
  ?sitelink
  ?item ?itemLabel
  ?go
  ?gene_symbol 
WHERE 
{
  ?item wdt:P686 ?go . 
  ?sitelink schema:about ?item .
  ?sitelink schema:isPartOf <https://en.wikipedia.org/> .
  
  ?protein wdt:{{ property }} ?item .
  ?protein wdt:P703 wd:{{ organism }} . 
  
  ?protein wdt:P702 ?gene . 
  
  ?gene wdt:{{ gene_symbol }} ?gene_symbol . 
            
  ?item rdfs:label ?itemLabel .
  FILTER (LANG (?itemLabel) = "en")
}
"""

# SPARQL query for human cell type markers
human_cell_type_query = """
SELECT DISTINCT
  ?sitelink
  ?item ?itemLabel
  ?gene_symbol 
WHERE 
{
  ?item wdt:P8872 ?gene .  
  ?item wdt:P703 wd:Q15978631 . 
  {?item wdt:P279 ?superclass .}
  UNION
  {?item wdt:P10019 ?superclass .}
  ?sitelink schema:about ?superclass .
  ?sitelink schema:isPartOf <https://en.wikipedia.org/> .
  
  ?gene wdt:P353 ?gene_symbol . 
            
  ?item rdfs:label ?itemLabel .
  FILTER (LANG (?itemLabel) = "en")
}
"""

# SPARQL query for mouse cell type markers
mouse_cell_type_query = """
SELECT DISTINCT
  ?sitelink
  ?item ?itemLabel
  ?gene_symbol 
WHERE 
{
  ?item wdt:P8872 ?gene .  
  ?item wdt:P703 wd:Q83310 . 
  {?item wdt:P279 ?superclass .}
  UNION
  {?item wdt:P10019 ?superclass .}
  ?sitelink schema:about ?superclass .
  ?sitelink schema:isPartOf <https://en.wikipedia.org/> .
  
  ?gene wdt:P2394 ?gene_symbol . 
            
  ?item rdfs:label ?itemLabel .
  FILTER (LANG (?itemLabel) = "en")
}
"""

# SPARQL query for molecular functions (Wikipedia links)
molecular_function_wikipedia_query_template = """
SELECT DISTINCT
  ?sitelink
  ?item ?itemLabel
  ?go
WHERE 
{
   ?something p:P680 ?statement . 
    ?statement ps:P680 ?item ;
                pq:P4390 wd:Q39893449 .
  
  # Get Wikipedia sitelinks
  ?sitelink schema:about ?something ;
            schema:isPartOf <https://en.wikipedia.org/>  .
  
  ?item wdt:P686 ?go . 
  ?item rdfs:label ?itemLabel .
  FILTER (LANG (?itemLabel) = "en")
}
"""

# SPARQL query for molecular functions (gene mappings)
molecular_function_gene_query_template = """
SELECT DISTINCT
  ?item ?itemLabel
  ?go
  ?gene_symbol 
WHERE 
{
  ?item wdt:P686 ?go . 

  ?protein wdt:P680 ?item .
  ?protein wdt:P703 wd:{{ organism }}. 
  
  ?protein wdt:P702 ?gene . 
  
  ?gene wdt:{{ gene_symbol }} ?gene_symbol . 
 
  ?item rdfs:label ?itemLabel .
  FILTER (LANG (?itemLabel) = "en")
}
"""


def fetch_data(query):
    sparql = SPARQLWrapper("https://query.wikidata.org/sparql")
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results["results"]["bindings"]


def process_data(results, include_go=True):
    data = []
    for result in results:
        sitelink = result.get("sitelink", {}).get("value", "")
        item = result["item"]["value"]
        itemLabel = result["itemLabel"]["value"]
        go = result.get("go", {}).get("value", "") if include_go else ""
        gene_symbol = result.get("gene_symbol", {}).get("value", "")
        data.append(
            {
                "sitelink": sitelink,
                "item": item,
                "itemLabel": itemLabel,
                "go": go,
                "gene_symbol": gene_symbol,
            }
        )

    return pd.DataFrame(data)


def generate_gmt(df, output_file, use_item_label=False):
    with open(output_file, "w") as f:
        group_col = "itemLabel" if use_item_label else "go"
        for term, group in df.groupby(group_col):
            genes = group["gene_symbol"].tolist()
            line = (
                f"{term}\t{group['itemLabel'].iloc[0]}\t{group['sitelink'].iloc[0]}\t"
                + "\t".join(genes)
                + "\n"
            )
            f.write(line)


def save_processes(df, output_file):
    processes = df.groupby("itemLabel")["gene_symbol"].apply(list).reset_index()
    processes.to_json(output_file, orient="records")


def main():
    organisms = {
        "human": {"id": "Q15978631", "gene_symbol": "P353"},
        "mouse": {"id": "Q83310", "gene_symbol": "P2394"},
    }

    go_properties = {
        "biological_processes": "P682",
        "cellular_components": "P681",
    }

    template = Template(go_query_template)

    for organism, details in organisms.items():
        for category, prop in go_properties.items():
            query = template.render(
                property=prop,
                organism=details["id"],
                gene_symbol=details["gene_symbol"],
            )
            results = fetch_data(query)
            df = process_data(results)
            output_file = STATIC / f"gene_sets_{organism}_{category}.gmt"
            generate_gmt(df, output_file)
            save_processes(df, STATIC / f"processes_{organism}_{category}.json")

    # Human cell type markers
    human_cell_type_results = fetch_data(human_cell_type_query)
    human_cell_type_df = process_data(human_cell_type_results)
    generate_gmt(
        human_cell_type_df,
        STATIC / "/gene_sets_human_cell_type.gmt",
        use_item_label=True,
    )
    save_processes(human_cell_type_df, STATIC / "/processes_human_cell_type.json")

    # Mouse cell type markers
    mouse_cell_type_results = fetch_data(mouse_cell_type_query)
    mouse_cell_type_df = process_data(mouse_cell_type_results)
    generate_gmt(
        mouse_cell_type_df,
        STATIC / "gene_sets_mouse_cell_type.gmt",
        use_item_label=True,
    )
    save_processes(mouse_cell_type_df, STATIC / "processes_mouse_cell_type.json")

    # Custom handling for molecular functions
    for organism, details in organisms.items():
        wikipedia_query = Template(molecular_function_wikipedia_query_template).render(
            organism=details["id"]
        )
        gene_query = Template(molecular_function_gene_query_template).render(
            organism=details["id"], gene_symbol=details["gene_symbol"]
        )

        wikipedia_results = fetch_data(wikipedia_query)
        gene_results = fetch_data(gene_query)

        df_wikipedia = process_data(wikipedia_results, include_go=False)
        df_gene = process_data(gene_results)

        combined_df = pd.merge(df_gene, df_wikipedia, on="item", suffixes=("", "_wiki"))
        combined_df["sitelink"] = combined_df["sitelink_wiki"]
        combined_df = combined_df.drop(columns=["sitelink_wiki", "itemLabel_wiki"])

        output_file = STATIC / f"gene_sets_{organism}_molecular_functions.gmt"
        generate_gmt(combined_df, output_file)
        save_processes(
            combined_df, STATIC / f"processes_{organism}_molecular_functions.json"
        )


if __name__ == "__main__":
    main()
