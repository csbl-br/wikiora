import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON
import random

# SPARQL query
query = """
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


def fetch_data():
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
        entrez = result["entrez"]["value"]
        ensembl_gene_id = result["ensembl_gene_id"]["value"]

        data.append(
            {
                "sitelink": sitelink,
                "item": item,
                "itemLabel": itemLabel,
                "go": go,
                "gene_symbol": gene_symbol,
                "entrez": entrez,
                "ensembl_gene_id": ensembl_gene_id,
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
    results = fetch_data()
    df = process_data(results)
    generate_gmt(df, "gene_sets.gmt")
    save_processes(df, "processes.json")


if __name__ == "__main__":
    main()
