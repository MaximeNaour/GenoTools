# GenoCompare
Script allowing to compare multiple genomes

The strain2ref.py script uses NCBI's EDirect tool to retrieve the corresponding NCBI taxonomy ID, the latest GenBank and RefSeq assembly accession numbers, and the URL for a given organism and strain name and writes them to a file.

The dl_genomes.py script allows downloading and creating directories for each taxon from strain2ref.py output file.
