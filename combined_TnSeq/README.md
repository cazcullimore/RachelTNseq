# TNseq pipeline
this is basically Tyler's setup with a couple modifications
- I added some SLURM submission scripts so that it runs everything together
- It removes the 
- All my preprocessing is done in preprocess.sh (the other files are just helpers to run on slurm)


## Preprocessing
1. we remove the 3' adapter and any 3' C or G tailing with cut adapt
    - min 16 length of trimming (pipeline truncates to 16 later)
1. it has an option to do 5' trimming as well but that isn't used
1. then it finds the transposon, trims it and removes reads withuot the transposon.
1. We quality filter (90% of bases in a read have to have <1% chance of being incorrect)

## PTT file generation
I added code that generates a PTT file from a GFF since PTT is ancient. This file is used to know where the genes start and end.

## Existing pipeline (TnSeq-Pipeline_SKCmodifications)
1. aligns reads against the reference
1. maps read counts to the genes
1. normalizes counts to intergeneic inserts as well as between samples (CHECK)

## Dependencies and environment setup

`conda activate TnSeqPreProcess` selects an existing environment; it does not
create it or install packages. The preprocessing scripts expect this exact name.
Create it once before submitting jobs:

```bash
module load miniforge3  # Cluster-specific module name used by these scripts.
conda create -n TnSeqPreProcess -c conda-forge -c bioconda python=3 cutadapt fastx_toolkit
conda activate TnSeqPreProcess
```

| Dependency | Used for |
| --- | --- |
| Cutadapt (`cutadapt`) | Removing adapters, poly-C/poly-G tails, and the transposon. |
| FASTX-Toolkit (`fastx_toolkit` Conda package) | Provides `fastq_quality_filter` for read quality filtering. |
| Bowtie2 (`bowtie2` and `bowtie2-build`) | Building the reference index and aligning reads. |
| Python 2.7 | Running the legacy `TnSeq-Pipeline_SKCmodifications.py`; Python 3 cannot run it unchanged. |
| Python 3 | Running `gff_to_ptt.py` to convert annotations. |
| Bash, gzip, and standard Unix commands | Running the shell scripts, handling compressed files, and managing outputs. |
| Slurm (`sbatch`) | Submitting preprocessing jobs and the dependent combined analysis job. |
| Conda and the cluster's `module` command | Loading software and activating the environments in the job scripts. |

### Combined analysis environment

`analyze_samples.sbatch` activates a **different**, case-sensitive environment
named `Tnseq` and loads the `bowtie2` module. A minimal legacy environment is:

```bash
conda create -n Tnseq -c conda-forge python=2.7
```

The analysis job must have **both `python2.7` and `python3` on PATH after
activation**, as well as `bowtie2` and `bowtie2-build`. The Python pipeline and
GFF converter use only their respective Python standard libraries and the local
`objects.py`; they do not need extra Python packages. Conda environments normally
contain one Python version, so creating `Tnseq` alone does not provide Python 3.
Make Python 3 available through your cluster's software setup, or adjust the
`python3` invocation in `analyze_samples.sbatch` to an installed Python 3 path.

The module names `miniforge3` and `bowtie2` are cluster-specific. Update the
`module load` lines if your cluster uses different names. Slurm and the module
system are supplied by the cluster; installing Python packages does not provide
them. Python 2.7 is a legacy requirement; these setup commands are not a pinned
or tested environment specification.

### Optional downstream analysis and plotting

The scripts in the parent directory use additional packages. These are separate
from preprocessing and should run in a modern Python 3 environment:

| Scripts or notebooks | Additional Python packages |
| --- | --- |
| Plotting and replicate/insertion analysis scripts | `matplotlib`, `numpy`, `pandas`, `scipy` (requirements vary by script). |
| `plot_zero_count_genes.py` | `matplotlib`; optional `adjustText` for label collision adjustment. |
| `run_deseq2.py` | `pydeseq2`, `scikit-learn`, `matplotlib`, `numpy`, `pandas`. Uses PyDESeq2; no R installation is invoked. |
| `run_gseapy_kegg.py` | `gseapy`, `numpy`, `pandas`; internet access to the KEGG REST API. |
| `rep_sim.ipynb`, `pathways.ipynb` | A Jupyter environment; `pandas`, plus `numpy`, `matplotlib`, `seaborn`, `scipy`, and `scikit-learn` for the replicate notebook. |

For these downstream tools, create a separate environment:

```bash
conda create -n TnSeqAnalysis -c conda-forge python=3 pip
conda activate TnSeqAnalysis
python -m pip install matplotlib numpy pandas scipy scikit-learn pydeseq2 gseapy adjustText seaborn jupyterlab
```

`../requirements_plotting.txt` lists only `matplotlib` and `adjustText`; it is
not a dependency list for the entire pipeline or all downstream scripts.



# Usage
  `bash path/to/submit_samples.sh --output-dir OUTPUT_DIR REFERENCE.fasta[.gz] ANNOTATION.gff[.gz] READS.fastq[.gz]...`

Example:
  `bash submit_samples.sh --output-dir experiment_1 \
    ref.fasta.gz genes.gff.gz reads/*.fastq.gz`

## Plain English Usage
`bash` runs terminal commands stored in a file
`submit_samples.sh` this is the file we are running. You need to provide the path to the file or else the terminal assumes the file is in your current folder

`--output-dir` - this is a marker that the next word will be the output folder

`OUTPUT_DIR` - this is a placeholder which you should replace with the output folder you want

`REFERENCE.fasta` is the path to your reference fasta file

`ANNOTATION.gff` is the path to your reference annotation file (needs to match the fasta file)

`READS.fastq..` add all the files you want to anlzye separated by spaces or just use a `*` to have the terminal find the paths that match. `*` matches anything so `*.fastq.gz` would match all the files in your current folder that end in `.fastq.gz`

## usage note
Each sample is preprocessed in a separate Slurm job. One combined analysis job starts only after every preprocessing job succeeds, preserving multi-sample normalization in the legacy Python pipeline. OUTPUT_DIR must not already exist.

All inputs can either be gzipped or not compressed. Include all the input files you want to run. 

The pipeline will normalize things together so run all samples together



# What I ran
I ran two main versions:
- one with all samples and one dropping the low read count sample
   - the low read count sample messes up the normalization of the pipeline because it rounds to the nearest integer. I removed the rounding but haven't rerun things (it shouldn't matter since I run DEseq on the unnormalized counts)
- Each of these I ran aligned to M1 or M108
    - prbably need to run the combined run (M1 to M1 and M108 to M108)
- I ran DEseq on both
- I also checked swapping samples but it doesn't really improve the correlation

## DEseq analysis
- compare HIS vs NHS 
    - I did run a paired version but that is incorrect because the replicates aren't done in different batches. 


# Technical Notes
## modifications to TnSeq-Pipeline_SKCmodifications
By default, Bowtie2 now aligns the filtered FASTQs from `preprocess.sh` directly,
including `.fastq.gz` files. This preserves the remaining read lengths and
Phred+33 base-quality scores; the Python pipeline skips its internal transposon
removal, read truncation, reverse complementation, and FASTA conversion.
External preprocessing still runs when using `submit_samples.sh`.

The generated `tnSeq.config` includes these settings:

```ini
[parameters]
Transposon =
ReadLength = -1

[options]
DirectFastq = True
ReverseComplementReads = False
```

`DirectFastq` defaults to `True` when omitted. It requires the other settings
shown above and FASTQ input; conflicting settings cause an error. Set
`DirectFastq = False` to use the legacy internal preprocessing and FASTA
alignment path. When submitting with `analyze_samples.sbatch`, change its config
template rather than an old generated config.

Direct alignment does not change the MAPQ >=40 filter, gene-end trimming,
minimum insertion-count filter, or normalization. Gene counting still runs
after alignment. Existing input FASTQs are never deleted by analysis cleanup.
