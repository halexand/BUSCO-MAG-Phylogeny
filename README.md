# Concatenated BUSCO Phylogeny for Eukaryotic MAG taxonomy

Workflow that based loosely off of the publication:  BUSCO applications from quality assessments to gene prediction and phylogenomics. Robert M. Waterhouse, Mathieu Seppey, Felipe A. Simão, Mose Manni, Panagiotis Ioannidis, Guennadi Klioutchnikov, Evgenia V. Kriventseva, and Evgeny M. Zdobnov. Mol Biol Evol, published online Dec 6, 2017. doi: 10.1093/molbev/msx319.

Includes:
- `BUSCO-Gene-Selection.ipynb`: Jupyter notebook to enable the selection of well covered BUSCO orthologs from a given set of MAGs. The notebook then pulls the orthologs from a given set of MAGs and reference genomes and create combined fasta files. 
- `Snakefile`: Snakefile that automates alignment of each individual fasta file with `MAAFT`, trimming with `trimal`, and phylogeny generation with `FastTree` or `RAxML`. 
- `submit_snakemake.sh`: Sumission file for slurm operations
- `envs/`: Required environment fiels
- `cluster.yaml`: Slurm details
- `config.yaml`: configuration file for the Snakefile that defines the extracted protein file locations
- `reference-genomes/`:
    - `Snakefile`: Run the BUSCO alignment of reference genomes. Currently designed to work with predicted sets of proteins. 
    - `genomes/`: All reference genomes listed for the project



Note: predicted proteins from reference genomes pulled from JGI with [jgi-query](https://github.com/glarue/jgi-query).  
 
