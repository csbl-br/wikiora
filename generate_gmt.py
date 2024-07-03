import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON

# SPARQL query for human genes
human_query = """
SELECT DISTINCT
  ?sitelink
  ?item ?itemLabel
  ?go
  ?gene_symbol 
  ?entrez 
  ?ensembl_gene_id 
WHERE 
{
  ?item wdt:P686 ?go . 
  ?sitelink schema:about ?item .
  ?sitelink schema:isPartOf <https://en.wikipedia.org/> .
  
  ?protein wdt:P682 ?item .
  ?protein wdt:P703 wd:Q15978631 . 
  
  ?protein wdt:P702 ?gene . 
  
  ?gene wdt:P353 ?gene_symbol . 
  ?gene wdt:P351 ?entrez .
  ?gene wdt:P594 ?ensembl_gene_id . 
            
  ?item rdfs:label ?itemLabel .
  FILTER (LANG (?itemLabel) = "en")
}
"""

# SPARQL query for mouse genes
mouse_query = """
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
  
  ?protein wdt:P682 ?item .
  ?protein wdt:P703 wd:Q83310 .  # MGI for mouse
  
  ?protein wdt:P702 ?gene . 
  
  ?gene wdt:P2394 ?gene_symbol .  # MGI gene symbol
            
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


def process_data(results):
    data = []
    for result in results:
        sitelink = result["sitelink"]["value"]
        item = result["item"]["value"]
        itemLabel = result["itemLabel"]["value"]
        go = result["go"]["value"]
        gene_symbol = result["gene_symbol"]["value"]
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


def generate_gmt(df, output_file):
    with open(output_file, "w") as f:
        for go_term, group in df.groupby("go"):
            genes = group["gene_symbol"].tolist()
            line = (
                f"{go_term}\t{group['itemLabel'].iloc[0]}\t{group['sitelink'].iloc[0]}\t"
                + "\t".join(genes)
                + "\n"
            )
            f.write(line)


def save_processes(df, output_file):
    processes = df.groupby("itemLabel")["gene_symbol"].apply(list).reset_index()
    processes.to_json(output_file, orient="records")


def main():
    # Human genes
    human_results = fetch_data(human_query)
    human_df = process_data(human_results)
    generate_gmt(human_df, "gene_sets_human.gmt")
    save_processes(human_df, "processes_human.json")

    # Mouse genes
    mouse_results = fetch_data(mouse_query)
    mouse_df = process_data(mouse_results)
    generate_gmt(mouse_df, "gene_sets_mouse.gmt")
    save_processes(mouse_df, "processes_mouse.json")


if __name__ == "__main__":
    main()
