# WikiORA - Gene Over-Representation Analysis

![WikiORA Example Workflow](https://wikiora.toolforge.org/static/example.png)

WikiORA is a tool designed to simplify the process of gene set over-representation analysis by integrating data from Wikidata and Wikipedia. Our aim is to provide an easy-to-use platform for researchers to identify significantly enriched gene sets in their data, using a combination of curated gene sets from various sources.

The tool provides a learning dimension during exploratory analysis, aiding bioinformaticians in making sense out of gene lists. 

## Features

- Integrates information from Wikidata, which has been enriched with Wikipedia, Gene Ontology, and PanglaoDB
- Supports human and mouse gene sets
- Provides over-representation analysis with hypergeometric test and Bonferroni correction
- Interactive results with links to Wikipedia pages for enriched terms
- Downloadable gene set files (GMT format)

## How It Works

WikiORA uses the following steps to perform over-representation analysis:

1. **Input:** A list of genes is provided by the user.
2. **Background Gene Sets:** The background gene sets are defined using data curated into Wikidata.
3. **Overlap Calculation:** For each gene set, the overlap between the user-provided gene list and the genes associated with the gene set (and its Wikipedia page) is calculated.
4. **p-value Calculation:** The p-value is calculated using the hypergeometric test, representing the probability of observing at least as many overlapping genes by chance.
5. **Correction:** The Bonferroni correction is applied to account for multiple testing and adjust the p-values.
6. **Results:** Results are sorted by p-value to highlight the most significantly over-represented terms.

## Data Sources

WikiORA uses Wikidata as a data source for the gene sets. It combines community curation with imports from sources such as:

- [Wikipedia](https://en.wikipedia.org)
- [Gene Ontology](http://geneontology.org)
- [PanglaoDB](https://panglaodb.se)

## Citing

(Manuscript in preparation)

While gene sets include more information than the original data sources, when using the cell type marker data, we recommend also citing [PanglaoDB](https://panglaodb.se). For gene ontology gene sets, we recommend citing also the [Gene Ontology Annotation (GOA) Database](https://www.ebi.ac.uk/GOA/) and the [Gene Ontology Resource](https://geneontology.org/docs/go-citation-policy/).

## Team

WikiORA is developed in Brazil by a team of bioinformaticians passionate about open knowledge. The project is led by [Tiago Lubiana](https://tiago.bio.br) at the [Computational Systems Biology Laboratory](https://www.csbiology.org/), headed by Prof. Helder Nakaya.

## Contact Us

If you have any questions, feedback, or suggestions, please feel free to contact us via [GitHub](https://github.com/lubianat/wikiora/issues).


## Usage and installation

WikiORA is available as a web-server at https://wikiora.sysbio.tools.

To run WikiORA locally, clone the repository and install the required dependencies:

```bash
git clone https://github.com/lubianat/wikiora.git
cd wikiora/www/python/src
pip install -r requirements.txt
```

Start the local server:

```bash
flask run
```

Open your web browser and go to `http://127.0.0.1:5000` to access WikiORA.

## Hosting

This project is hosted on Toolforge at [wikiora.toolforge.org](https://wikiora.toolforge.org).

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

Made in 🇧🇷 in 2024. Content from Wikidata is under CC0.

If you like our project, give us a ⭐ on [GitHub](https://github.com/lubianat/wikiora)!
