# Genomic characterization of human respiratory syncytial virus circulating in Islamabad, Pakistan, during an outbreak in 2022 to 2023

Reproducible workflow that mirrors the analysis in the study, with an explicit CoVpipe2-inspired short read pipeline and clear manual steps for annotation and tree building.

DOI: https://doi.org/10.1007/s00705-024-06036-0

## Program summary
End to end analysis for RSV-A and RSV-B. Steps are automated where useful and documented for manual replication when a GUI is simpler or when external resources are required.

## Wet lab and sequencing context
- Input: throat or nasopharyngeal swabs collected in VTM from RSV positive patients.
- Library prep: Nextera DNA Flex or Nextera XT. Follow vendor instructions for dual indexed libraries.
- Platform: Illumina MiSeq, 2x150 bp. Use a run configuration that yields enough depth per sample for whole genome consensus. The original study used MiSeq and short paired end reads; 2x150 works well in this pipeline.

## CoVpipe2-inspired short read path
Minimal, explicit path that mirrors CoVpipe2 style. This is the reference based consensus path for each sample.

Trimmomatic -> Picard MarkDuplicates -> BWA-MEM alignment to reference -> samtools sort and index -> bcftools mpileup and call -> masked consensus

### Requirements
- FastQC v0.11.9
- Trimmomatic v0.39
- bwa, samtools, bcftools
- Picard
- Python 3.11+
- Optional: Snakemake, Conda, Docker

### References for mapping
- RSV-A: NC_038235.1
- RSV-B: NC_001781.1

### Quick start one command
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate rsv-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

### Manual pipeline for a single sample
Assume these variables:
```bash
SAMPLE="RSV_A_001"
SUBTYPE="A"     # A or B
R1="data-private/${SAMPLE}_R1.fastq.gz"
R2="data-private/${SAMPLE}_R2.fastq.gz"

# Adapter file for Nextera
ADAPTERS="env/NexteraPE-PE.fa"

# Reference by subtype
REF_A="refs/NC_038235.1.fasta"
REF_B="refs/NC_001781.1.fasta"
REF=$( [ "$SUBTYPE" = "A" ] && echo "$REF_A" || echo "$REF_B" )

MIN_DEPTH=10
MIN_QUAL=20
OUTDIR="results/${SAMPLE}"
mkdir -p "$OUTDIR"
```

1) Quality control and trimming
```bash
fastqc -o "$OUTDIR" "$R1" "$R2"

trimmomatic PE -threads 4   "$R1" "$R2"   "$OUTDIR/${SAMPLE}_R1.trim.fastq.gz" "$OUTDIR/${SAMPLE}_R1.unpaired.fastq.gz"   "$OUTDIR/${SAMPLE}_R2.trim.fastq.gz" "$OUTDIR/${SAMPLE}_R2.unpaired.fastq.gz"   ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```

2) Index reference once
```bash
mkdir -p refs/index
cp "$REF" refs/index/ref.fasta
bwa index refs/index/ref.fasta
samtools faidx refs/index/ref.fasta
picard CreateSequenceDictionary R=refs/index/ref.fasta O=refs/index/ref.dict
```

3) Alignment, sorting, duplicate marking
```bash
bwa mem -t 4 refs/index/ref.fasta   "$OUTDIR/${SAMPLE}_R1.trim.fastq.gz" "$OUTDIR/${SAMPLE}_R2.trim.fastq.gz"   | samtools sort -o "$OUTDIR/${SAMPLE}.sorted.bam"

picard MarkDuplicates I="$OUTDIR/${SAMPLE}.sorted.bam" O="$OUTDIR/${SAMPLE}.dedup.bam" M="$OUTDIR/${SAMPLE}.dupmetrics.txt" REMOVE_DUPLICATES=true
samtools index "$OUTDIR/${SAMPLE}.dedup.bam"
```

4) Depth based mask and consensus creation
```bash
# Depth
samtools depth -a "$OUTDIR/${SAMPLE}.dedup.bam" > "$OUTDIR/${SAMPLE}.depth.txt"

# Positions with coverage < MIN_DEPTH masked to N
awk -v min="$MIN_DEPTH" 'BEGIN{OFS="\t"} {if($3<min) print $1,$2-1,$2}' "$OUTDIR/${SAMPLE}.depth.txt" > "$OUTDIR/${SAMPLE}.mask.bed"

# Variant calling
bcftools mpileup -f refs/index/ref.fasta -Q "$MIN_QUAL" -Ou "$OUTDIR/${SAMPLE}.dedup.bam"   | bcftools call -mv -Oz -o "$OUTDIR/${SAMPLE}.vcf.gz"
bcftools index "$OUTDIR/${SAMPLE}.vcf.gz"

# Consensus with masked low depth
cat refs/index/ref.fasta   | bcftools consensus -m "$OUTDIR/${SAMPLE}.mask.bed" -f - "$OUTDIR/${SAMPLE}.vcf.gz"   > "$OUTDIR/${SAMPLE}.consensus.fasta"
```

Outputs
- Consensus FASTA: `results/<sample>/<sample>.consensus.fasta`
- BAM plus QC: `results/<sample>/*.bam, *.dupmetrics.txt, *.depth.txt`
- VCF: `results/<sample>/<sample>.vcf.gz`

## Annotation with Nextclade
You can do this manually in the web app or via CLI.

### Web workflow
1. Open https://clades.nextstrain.org.
2. Choose an RSV dataset that matches your subtype.
   - RSV A dataset for A samples.
   - RSV B dataset for B samples.
3. Upload your `*.consensus.fasta` files.
4. Review clade and genotype assignment, QC, and mutation calls.
5. Export CSV and aligned FASTA for your records.

### CLI workflow
```bash
# Get datasets
nextclade dataset get --name rsv-a --output-dir datasets/rsv-a
nextclade dataset get --name rsv-b --output-dir datasets/rsv-b

# Run on your consensus
nextclade run --input-dataset datasets/rsv-a --output-all results/nextclade/a results/consensus_rsv_a.fasta
nextclade run --input-dataset datasets/rsv-b --output-all results/nextclade/b results/consensus_rsv_b.fasta
```

## Context sequences and tree building
Follow the manuscript criteria to curate context sequences and then build a maximum likelihood tree with your samples plus context.

### Manual BLAST to collect context
1. Go to https://blast.ncbi.nlm.nih.gov and select **nucleotide BLAST**.
2. Upload your consensus FASTA.
3. Database: **nt**. Program selection: **Highly similar sequences (megablast)**.
4. Organism filter: for RSV-A use **Human orthopneumovirus A**. For RSV-B use **Human orthopneumovirus B**.
5. Increase **Max target sequences** to a generous number like 1000.
6. Run BLAST and sort by percent identity and date. Download FASTA for a reasonable set of close matches.
7. Apply manuscript-like filters before finalizing the context set:
   - Remove sequences with more than 2 percent ambiguous bases (Ns).
   - Remove nearly identical duplicates from the same place and time window.
   - Prefer complete or near complete genomes and contemporary sampling windows relevant to your data.
8. Save the curated context FASTA as `results/context/context_curated.fasta`.

### Align and build the tree
```bash
# Combine your sample consensuses and curated context
cat results/*/*.consensus.fasta > results/consensus/all_consensus.fasta
cat results/consensus/all_consensus.fasta results/context/context_curated.fasta > results/aln/combined.fasta

# MAFFT alignment
mkdir -p results/aln
mafft --auto results/aln/combined.fasta > results/aln/combined.aln.fasta

# IQ-TREE with GTR+G model and 1000 ultrafast bootstraps
mkdir -p results/iqtree
iqtree2 -s results/aln/combined.aln.fasta -m GTR+G -bb 1000 -nt AUTO -pre results/iqtree/rsv_tree

# Outputs
# - results/iqtree/rsv_tree.treefile
# - results/iqtree/rsv_tree.iqtree
# - results/iqtree/rsv_tree.log
```

### Optional visualization
- Load `results/iqtree/rsv_tree.treefile` in FigTree or iTOL.
- Color samples by subtype or location for clarity.

## Mutational analysis for G and F
Generate amino acid differences vs the appropriate reference to replicate the manuscript style tables.

```bash
# RSV-A
python analysis/scripts/mutations.py --subtype A --genomes results/consensus/all_consensus.fasta --out_tsv results/mutations_rsv_a.tsv
# RSV-B
python analysis/scripts/mutations.py --subtype B --genomes results/consensus/all_consensus.fasta --out_tsv results/mutations_rsv_b.tsv
```

The scripts translate G and F coding regions, compare to the reference for the subtype, and export a TSV with nonsynonymous changes.

## Outputs checklist
- results/<sample>/<sample>.consensus.fasta
- results/consensus/all_consensus.fasta
- results/context/context_curated.fasta
- results/aln/combined.aln.fasta
- results/iqtree/rsv_tree.treefile and associated files
- results/mutations_rsv_*.tsv

## Configuration
Edit `config/config.yaml` for batch runs. Example:
```yaml
pairs:
  - sample: RSV_A_001
    subtype: A
    r1: data-private/RSV_A_001_R1.fastq.gz
    r2: data-private/RSV_A_001_R2.fastq.gz
  - sample: RSV_B_001
    subtype: B
    r1: data-private/RSV_B_001_R1.fastq.gz
    r2: data-private/RSV_B_001_R2.fastq.gz

ref_accessions:
  A: NC_038235.1
  B: NC_001781.1

params:
  threads: 4
  min_depth_consensus: 10
  min_qual: 20
  iqtree_model: GTR+G
  bootstrap: 1000
```

## How to cite
- Cite the manuscript: DOI 10.1007/s00705-024-06036-0
- Cite the tools you used in your analysis pipeline.

## Automated pipeline with Snakemake

- Place your FASTQs in `fastq/` using the pattern `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`.
- Place your reference FASTA files in `refs/` as `RSV_A.fasta` and `RSV_B.fasta` or change the paths in `config/config.yaml`.
- Choose how to assign references to samples in `config/config.yaml`:
  - `single`: all samples use one subtype.
  - `auto`: infer subtype from tokens in the sample name.
  - `map`: explicit sample to subtype mapping in `config/samples.tsv`.

### Run
```bash
conda env create -f env/environment.yml
conda activate rsv-env
snakemake -s workflow/Snakefile -c 4 --configfile config/config.yaml --printshellcmds
```

Outputs will be in `results/`. If you set `tree.context_fasta`, MAFFT and IQ-TREE will also run automatically to build a tree.
