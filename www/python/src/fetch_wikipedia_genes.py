import requests
import json
import time
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


def fetch_sparql_results(endpoint, query, retries=3):
    headers = {
        "Accept": "application/sparql-results+json",
        "User-Agent": "Mozilla/5.0 (compatible; myscript/1.0; +https://example.com/bot)",
    }
    for _ in range(retries):
        try:
            response = requests.get(endpoint, headers=headers, params={"query": query})
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError as e:
            if response.status_code == 403:
                print(f"Access denied: {e}. Retrying...")
            else:
                print(f"HTTP error occurred: {e}")
            time.sleep(5)  # Wait before retrying
    raise requests.exceptions.HTTPError(
        f"Failed to fetch data after {retries} attempts"
    )


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
        json.dump(genes, f, indent=4, sort_keys=True


if __name__ == "__main__":
    main()
