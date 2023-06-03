## Dovetail Clover Genome (DCG)

- TODO: Link manuscript/preprint when available
- TODO: Link to Haplotype Fasta and AGP files

### Manuscript: Chromosome-level, haplotype-resolved assembly of white clover (*Trifolium repens*)
#### Authors: James S. Santangelo, Paul Battlay, Brandon Hendrickson, Wen-Hsi Kuo, Kenneth M. Olsen, Nicholas J. Kooyers, Johnson, M.T. Johnson, Kathryn A. Hodgins, Rob W. Ness



### Description

This repository contains the code used for the assembly and annotation of a
chromosome-level, haplotype-resolved assembly of white clover (*Trifolium
repens*). The pipeline begins with input of the Dovetail haplotype assemblies,
associated AGP (i.e., “A Golden Path”) files (LINK), and linkage map data, 
and ends with the generation of the phased diploid assembly in FASTA format
(NCBI BioProjects PRJNA957817 and PRJNA957816), the annotated haploid mapping
assembly in FASTA, NCBI Sequin, and GFF3 formats (BioProject PRJNA951196), and
manuscript figures.

Graphical depiction of pipeline can be found [here](./workflow/rulegraph.pdf).
Docstrings for each `Snakemake` rule can be found [here](./docstrings.txt)

In addition to the raw dovetail haplotype FASTAs, the following
dependencies/resources need to be obtained prior to running the pipeline.

1. `GeneMark-ES/ET/EP+ ver 4.71_lic` needs to be obtained from
   [here](http://topaz.gatech.edu/Genemark/license_download.cgi). `GeneMark` is
   licensed software and could not be included as part of the pipeline. Once
   downloaded, install `GeneMark` as recommended in their README file. Once
   installed, the following needs to be done:  
    1. Include the `.gm_key` in the [workflow/](./workflow) directory containing the main Snakefile used for running the pipeline
    2. Point to the `GeneMark` and `ProtHint` (included with `GeneMark`) programs in the [config](./config/hpcnode.yaml) file
2. A copy of the RepBase repeat library (I used
   `RepBaseRepeatMaskerEdition-20181026.tar.gz`) from
   [here](https://www.girinst.org/repbase/). A license is required to use this
   database so it could not be included. Once obtained, it needs to be included
   in [resources/](./resources) directory
3. A copy of the `InterProScan` data needs to be obtained and included in the
   [resources/](./resources) directory. The following commands can be run to
   download and setup the database:

```
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.61-93.0/alt/interproscan-data-5.61-93.0.tar.gz
tar -xvzfp interproscan-data-5.61-93.0.tar.gz
chmod -R 777 interproscan-data-5.61-93.0
```

#### Using the pipeline

This pipeline requires `Conda` and `Singularity`:

A minimal installation of `Conda` (i.e., `Miniconda`) can be installed by
following the instructions for your platform
[here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Installation of `Singularity` requires Admin privileges, but using
`Singularity` to run pre-created containers does not. Installation instructions
can be found [here](https://docs.sylabs.io/guides/latest/admin-guide/). All
`Singularity` containers used in this pipeline are avalaible in [this public
repository](https://cloud.sylabs.io/library/james-s-santangelo), though they
will be automatically pulled and executed by the pipeline.
Assuming `Conda` is installed, the this repository's Conda environment can be
replicated by running the following command:

```
conda env create -f environment.yaml -n dcg 
```

This will create a `Conda` environment named dcgcontaining a minimal set of dependencies required to run the pipeline. 

After activating the environment (`conda activate dcg`), the pipeline can be
executed from the workflow directory by running a command that looks something
like:

```
snakemake --use-conda --use-singularity --singularity-args "--bind <path> --bind <path/to/interproscan/data>:/opt/interproscan-5.61-93.0/data" --configfile ../config/<configfile> --notemp -j <cores>
```

for local execution. Here, `<path>` is the path on the cluster from which files
will be read/written (e.g., `/scratch`), `<path/to/interproscan/data>` is
the full path to the interproscan data in [resources/](./resources),
`<configfile>` is one of the configfiles in the config directory that needs
to be modified to match the paths on your system, and `<cores>` is the
number of cores available for executing parallel processes.
