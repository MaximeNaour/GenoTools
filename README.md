# GenoCompare

```mermaid
%%{init: {'theme': 'base', 'themeVariables': { 'primaryColor': '#f8eeec', 'edgeLabelBackground':'#ffffff'}}}%%

flowchart TB;
linkStyle default interpolate basis
classDef exClass font-style:bold;
classDef exClass font-size:16px;

A("<font size=3>Identify genome/stain<br>script: strain2ref.py") --> B("<font size=3>Download genome from RefSeq/Assembly<br>script: dl_genomes.py");
B --> C("<font size=3>Gene identification (Prokka)");
C --> D("<font size=3>Gene annotation (UniProtKB/KEGG)");
D --> E("<font size=3>Pathway analysis (Pathway Tools)");
E --> F("<font size=3>Manual verification (Artemis)");
F --> G("<font size=3>Store annotation (GenBank)");

```

The "unique_seq.py" script allows comparing two or more genomes to identify unique portions in the genome of the organism of interest in comparison to the genomes of other organisms.

The "strain2ref.py" script uses NCBI's EDirect tool to retrieve the corresponding NCBI taxonomy ID, the latest GenBank and RefSeq assembly accession numbers, and the URL for a given organism and strain name and writes them to a file.

The "dl_genomes.py" script allows downloading and creating directories for each taxon from "strain2ref.py" output file.
