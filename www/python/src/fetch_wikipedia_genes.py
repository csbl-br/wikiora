import requests
import json
from pathlib import Path

HERE = Path(__file__).parent.resolve()
STATIC = HERE.joinpath("static").resolve()

SPARQL_ENDPOINT = "https://query.wikidata.org/sparql"
SPARQL_QUERY = """
SELECT ?gene ?geneLabel (IF(BOUND(?wikipediaPage), "blue", "red") AS ?pageStatus) 
  (IF(!BOUND(?wikipediaPage), URI(CONCAT("https://en.wikipedia.org/wiki/", ?geneLabel)), ?wikipediaPage) AS ?wikipediaLink)
WHERE {
  {?gene wdt:P353 ?geneLabel.}
  UNION
  {?mouse_gene wdt:P2394 ?geneLabel.
   ?mouse_gene wdt:P684 ?gene . 
   ?gene wdt:P353 ?any .  
  }
  
  OPTIONAL {
    {
      ?wikipediaPage schema:about ?gene;
                     schema:isPartOf <https://en.wikipedia.org/>.
    }
    UNION
    {
      ?gene wdt:P688 ?protein . 
      ?wikipediaPage schema:about ?protein;
                     schema:isPartOf <https://en.wikipedia.org/>.
    }
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
"""


def fetch_sparql_results(endpoint, query):
    headers = {"Accept": "application/sparql-results+json"}
    response = requests.get(endpoint, headers=headers, params={"query": query})
    response.raise_for_status()
    return response.json()


def process_results(results):
    genes = {}
    for result in results["results"]["bindings"]:
        gene_label = result["geneLabel"]["value"]
        gene_info = {
            "gene": result["gene"]["value"],
            "pageStatus": result["pageStatus"]["value"],
            "wikipediaLink": result["wikipediaLink"]["value"],
        }
        genes[gene_label] = gene_info
    return genes


def main():
    results = fetch_sparql_results(SPARQL_ENDPOINT, SPARQL_QUERY)
    genes = process_results(results)
    with open(STATIC / "genes.json", "w") as f:
        json.dump(genes, f, indent=4)


if __name__ == "__main__":
    main()
