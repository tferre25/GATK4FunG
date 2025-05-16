# Fungi Variant Calling

A reproducible Snakemake-based pipeline for calling variants in fungal sequencing data, from raw FASTQ to phylogenetic tree reconstruction.

## Description

This workflow automates standard steps in fungal variant calling:

1. **Quality control & trimming** of raw FASTQ files with Trimmomatic and FastQC.
2. **Alignment** of reads to a reference genome using BWA-MEM.
3. **Preprocessing**: marking adapters, adding read groups, sorting, and duplicate marking (Picard, Samtools).
4. **Variant discovery**: GATK HaplotypeCaller in GVCF mode, joint genotyping, and hard filtration.
5. **VCF post-processing**: genotype filters, conversion to FASTA (vcfSnpsToFasta).
6. **Phylogenetic inference** with IQ-TREE.

All tools are managed via Conda (`fungi` environment) to ensure reproducibility. Configuration is driven by a single `config.json` file in the project root.

## Usage

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/fungi-variant-calling.git
cd fungi-variant-calling
```

### 2. Prepare configuration

Edit `config.json` to set your input directory, output directory, reference genome, and tool options. Example:

```json
{
  "general_path": {
    "INPUT_PATH": "/path/to/raw_fastq",
    "OUTPUT_PATH": "/path/to/results",
    "filterGatkGenotypes_PATH": "/path/to/filterGatkGenotypes.py",
    "vcfToFasta_PATH": "/path/to/vcfSnpsToFasta.py"
  },
  "bwamem": {
    "DB": "/path/to/reference.fasta",
    "OPTIONS": "-M -t 8"
  },
  "trimmomatic": {
    "OPTIONS": "PE -threads 4 -phred33",
    "PARAMS": "ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
  },
  "iqtree": {
    "OPTIONS": "-m GTR+G -nt AUTO"
  }
}
```

### 3. Set up the Conda environment

```bash
conda env create --file environment.yml
conda activate fungi
```

### 4. Run the pipeline

From the project root:

```bash
snakemake --use-conda --cores 8
```

This will generate QC reports, BAMs, VCFs, and finally a phylogenetic tree (`<run_name>.treefile`).

### 5. Inspect results

* QC: `results/multiqc/multiqc_report.html`
* Variants: `results/<run_name>.hard_filtered.vcf.gz`
* Alignment stats & metrics: `results/logs/`
* Phylogeny: `results/<run_name>.treefile`

## Contributing

We welcome contributions! Please follow these steps:

1. Fork the repository and create your feature branch:

   ```bash
   git checkout -b feature/my-new-feature
   ```
2. Commit your changes with clear messages:

   ```bash
   git commit -m "Add support for multi-sample VCF merging"
   ```
3. Push to your fork and open a Pull Request on GitHub.
4. Ensure your code follows existing style (Snakefile formatting, JSON config structure) and include tests or example data if applicable.
5. We will review and provide feedback. Thank you for your contribution!

For bug reports and feature requests, please open an issue on the GitHub issue tracker.


## Contributing

We welcome contributions to our application! 

## Contact

If you have any questions or concerns about our application, please contact us at julien.robert@aphp.fr or theo.ferreira@aphp.fr

