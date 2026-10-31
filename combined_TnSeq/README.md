# TNseq pipeline
this is basically Tyler's setup with a couple modifications
- I added some SLURM submission scripts so that it runs everything together
- It removes the 
- All my preprocessing is done in preprocess.sh (the other files are just helpers to run on slurm)


## Preprocessing
    1. we remove the 3' adapter and any 3' C or G tailing with cut adapt
        - min 16 length of trimming
    1. it has an option to do 5' trimming as well but that isn't used
    1. then it finds the transposon, trims it and removes reads withuot the transposon.
    1. We quality filter (90% of bases in a read have to have <1% chance of being incorrect)

## ptt file generation
I added code that generates a PTT file from a GFF since PTT is ancient. This file is used to know where the genes start and end.

## Then everything is fed into TnSeq-Pipeline_SKCmodifications
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

1. aligns reads against the reference
1. maps read counts to the genes
1. normalizes counts to intergeneic inserts as well as between samples (CHECK)


it assumes a conda envirnoment called TnSeqPreProcess exists to have all the libraries (have chat add dpeendances)
conda activate TnSeqPreProcess



# Usage
  `bash path/to/submit_samples.sh --output-dir OUTPUT_DIR REFERENCE.fasta[.gz] ANNOTATION.gff[.gz] READS.fastq[.gz]...`

Example:
  bash submit_samples.sh --output-dir experiment_1 \
    ref.fasta.gz genes.gff.gz reads/*.fastq.gz

## plain english usage
`bash` runs terminal commands stored in a file
`submit_samples.sh` this is the file we are running. You need to provide the path to the file or else the terminal assumes the file is in your current folder

`--output-dir` - this is a marker that the next word will be the output folder
`OUTPUT_DIR` - this is a placeholder which you should replace with the output folder you want
`REFERENCE.fasta` - is the path to your reference fasta file
`ANNOTATION.gff` - is the path to your reference annotation file (needs to match the fasta file)
`READS.fastq..` - add all the files you want to anlzye separated by spaces or just use a `*` to have the terminal find the paths that match. `*` matches anything so `*.fasta.gz` would match all the files in your current folder that end in `.fasta.gz`

## usage note
Each sample is preprocessed in a separate Slurm job. One combined analysis job
starts only after every preprocessing job succeeds, preserving multi-sample normalization in the legacy Python pipeline. OUTPUT_DIR must not already exist.

All inputs can either be gzipped or not compressed. Include all the input files you want to run. 

The pipeline will normalize things together so run all the samples in your analysis



# What I ran
I ran two different versions:
- one with all samples and one dropping the low read count sample
   - the low read count sample messes up the normalization of the pipeline because it rounds to the nearest integer. I removed the rounding but haven't rerun things (it shouldn't matter since I run DEseq on the unnormalized counts)
- I ran a version aligned to M1 and M108
    - prbably need to run the combined run (M1 to M1 and M108 to M108)
- I ran DEseq on both
- I also checked swapping samples but it doesn't really improve the correlation

## DEseq analysis
- compare HIS vs NHS 
    - I did run a paired version but that is incorrect because the replicates aren't done in different batches and don't 

