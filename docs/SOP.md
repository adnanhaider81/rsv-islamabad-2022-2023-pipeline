# Standard operating procedure

This repository provides a reproducible Snakemake pipeline for RSV genomes.
You provide the references in `refs/` and place your FASTQ files in `fastq/`.
The pipeline autodetects samples by basename and pairs reads using suffixes.

## Autodetection
- Files must be named like `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`.
- The list of samples is computed from the R1 files in `fastq/`.
- Reference assignment supports three modes:
  - single: all samples use the same subtype reference.
  - auto: subtype inferred from tokens in the sample name, set in `config.yaml`.
  - map: override with `config/samples.tsv` mapping sample to subtype.

## Outputs
- `results/consensus/<sample>.consensus.fasta`
- `results/consensus/all_consensus.fasta`
- `results/qc/` FastQC reports and `multiqc_report.html`
- `results/bam/` sorted and deduplicated BAMs
- `results/vcf/` variant calls
- If you provide `tree.context_fasta`, alignment and tree are created.
