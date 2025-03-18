# Cell Knowledge Network Extraction, Translation, and Loading of Results

## Motivation

The Cell Knowledge Network (Cell KN) pilot aims to create a
comprehensive cell phenotype knowledge network that integrates
knowledge about diseases and drugs to facilitate discovery of new
biomarkers and therapeutic targets. The Cell KN captures single cell
genomics data from existing data repositories, such as CELLxGENE, and
uses NSForest to identify cell type-specific marker genes. The cell
types are manually mapped to the Cell Ontology (CL), and the marker
genes are linked to data from external sources, such as Open Targets,
to provide relationships to diseases and drugs. In addition, the Cell
KN uses natural language processing to extract information about cell
type-specific marker genes, and their association with disease state
from open access peer-reviewed publications.

To maximize interoperability of the Cell KN with with knowledge from
other NLM/NCBI resources, the knowledge will be derived in the form of
semantically-structured assertions of subject-predicate-object triple
statements which are compatible with storage using semantic web
technologies, and graph databases, such as the
[ArangoDB](https://arangodb.com/) database system.

The NCBI Information Resources Branch has extensive experience with
ArangoDB, including performance comparison testing with Neo4j,
interaction with the ArangoDB developers, and use in production.

## Purpose

The `cell-kn-etl-results` repository provides Python modules for
fetching data from external sources, and creating semantic triples
from NSForest results, their manual mapping to the CL, and external
data and NLP results. Note that these triples are created consistent
with the Cell KN schema, and a Python module is provided for creating
triples that represent the schema. Finally, a Java package is provided
for loading these triples into an ArangoDB instance. Note that this
package can accept quadruples which represent edge annotations.

## External Sources

Data can be fetched from the following external sources:

- [Open Targets](https://www.opentargets.org/): The
  [gget](https://pachterlab.github.io/gget/en/introduction.html)
  [opentargets](https://pachterlab.github.io/gget/en/opentargets.html)
  command can be used to obtain the `diseases`, `drugs`,
  `interactions`, `pharmacogenetics`, `tractability`, `expression`,
  and `depmap` Open Targets resources
- [Search API](https://www.ebi.ac.uk/ebisearch/documentation/rest-api): Hosted
  by the [European Bioinformatics Institute (EBI)](https://www.ebi.ac.uk/),
  part of the [European Molecular  Biology Laboratory (EMBL)](https://www.embl.org/),
  the search API can be used to obtain drug information
- [RxNav](https://lhncbc.nlm.nih.gov/RxNav/): The RxNav NLM resource,
  which supports RxNorm (the NLM standard terminology for drugs), can
  be used to obtain drug names, properties, and RXCUI identifiers

## Dependencies

### Submodule

The `cell-kn-etl-results` repository includes the
`cell-kn-etl-ontologies` repository as a submodule. After cloning
`cell-kn-etl-results`, initialize and update the submodule as follows:
```
git submodule init
git submodule update
```

### Java

Java SE 21 and Maven 3 or compatible are required to generate the
Javadocs, and package. Do each of these as follows:
```
$ mvn javadoc:javadoc
$ mvn clean package -DskipTests
```
Note that implementation of testing is planned for a later phase.

### Python

Python 3.12 and Poetry are required to generate the Sphinx
documentation, and run the modules. Install the dependencies as
follows:
```
$ python3.12 -m venv .poetry
$ source .poetry/bin/activate
$ python -m pip install -r .poetry.txt
$ deactivate
$ python3.12 -m venv .venv
$ source .venv/bin/activate
$ .poetry/bin/poetry install
```
Generate the Sphinx documentation as follows:
```
$ cd docs/python
$ make clean html
```
Again, implementation of testing is planned for a later phase.

### Data

The Python and Java classes require the ontology files to reside in
`data/obo`. Create and populate this directory as follows:
```
$ mkdir data/obo
$ cd src/main/python
$ python OntologyParserLoader.py --update
```

### Docker

Install [Docker Desktop](https://docs.docker.com/desktop/).

### ArangoDB

An ArangoDB docker image can be downloaded and a container started as
follows:
```
$ cd src/main/shell
$ export ARANGO_DB_HOME="<some-path>/arangodb"
$ export ARANGO_DB_PASSWORD="<some-password>"
$ ./start-arangodb.sh
```

## Usage

Run the Python modules for fetching data from external sources, and
creating semantic triples from NSForest results, their manual mapping
to the CL, external data and NLP results, and the Cell KN schema as follows
```
$ cd src/main/python
$ python ExternalApiResultsFetcher.py
$ python NSForestResultsTupleWriter.py
$ python AuthorToClResultsTupleWriter.py
$ python ExternalApiResultsTupleWriter.py
$ python CellKnSchemaTupleWriter.py
```
Note that quadruples are also created which represent edge
annotations.

Run the Java package for loading these tuples (triples and quadruples)
into an ArangoDB instance as follows:
```
$ export ARANGO_DB_HOST=127.0.0.1
$ export ARANGO_DB_PORT=8529
$ export ARANGO_DB_HOME="<some-path>/arangodb"
$ export ARANGO_DB_PASSWORD="<some-password>"
$ java -cp target/cell-kn-etl-ontologies-1.0.jar gov.nih.nlm.ResultsTupleLoader
```
