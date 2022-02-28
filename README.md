# Epigenomics genotyping pipeline

Nextflow pipeline for genotyping from epigenomics data

## Requirements
- Nextflow (https://www.nextflow.io/)
- samtools (http://www.htslib.org/)
- bcftools (http://www.htslib.org/)
- pyfaidx (https://github.com/mdshw5/pyfaidx)
- biopython (`conda install -c conda-forge biopython`)
- bedops (`conda install -c bioconda bedops`)
- plink2 (`conda install -c bioconda plink2`)
- packaging (`conda install -c conda-forge packaging`)

## Pipeline overview

Samples BAM files are merged by corresponding individual and then used for a ``bcftools``-based genotyping pipeline.

## Usage
```
[jvierstra@dev0 ~]$ nextflow run main.nf -config nextflow.config -profile Altius
```

## Input

<details><summary>Primary Sample file [--samples_file]</summary>
<p></p>
<p>
A tab-delimited file containing information about each sample. The file must contain a header and the following columns (other columns are permitted and ignored):

- **indiv_id**: Individual identifier for each sample; many samples can refer to one individual (if running for TF ChIPseq use the TF name)
- **bam_file**: Absolute path the BAM-formated file
- **ln_number**: Individual identifier for each file (if running for TF ChIPseq use the ENCODE dataset ID with the replicate number)
</p>
</details>

<details><summary>Filter Variants Sample file [--samples_file]</summary>
<p></p>
<p>
A tab-delimited file containing information about each sample. The file must contain a header and the following columns (other columns are permitted and ignored):

- **indiv_id**: Individual identifier for each sample; many samples can refer to one individual (if running for TF ChIPseq use the TF name)
- **cell_type**: Cell type or ENCODE accession number
- **hotspots_file**: Peak call file in BED format
</p>
</details>


<details><summary>Genome reference [--genome]</summary>
<p></p>
Reference genome file types:

  - `.fa`
  - `.nuclear.txt` (example in `supplemental_files` directory)
  - `.chrom_sizes`
  - BWA indexed (use the same prefix as for the `.fa` file)
<p></p>
</details>

<details><summary>dbSNP reference [--dbsnp_file]</summary>
<p></p>
Find the latest release of the dbSNP reference [here](https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz)
<p></p>
</details>

<details><summary>Ancestral genome [--genome_ancestral_fasta_file]</summary>
<p></p>
Ancestral genome file. Latest for GRCh38 available [here](http://ftp.ensembl.org/pub/release-86/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38_e86.tar.gz)
<p></p>
</details>

### Additonal Parameters:
<details><summary>Chunk size [--chunksize 5000000]</summary>
<p></p>
<p>Specificies the size (in base-pairs) to use when dividing the genome into chunks for parallel processing.</p>
</details>

<details><summary>SNP quality [--min_SNPQ 10]</summary>
<p></p>
<p>Filter variants with poor quality</p>
</details>

<details><summary>Genotype quality [--min_GQ 50]</summary>
<p></p>
<p>Set genotype for an individual to ./. (missing) when genotyping score (FORMAT/GQ) is less than this value.</p>
</details>

<details><summary>Sequencing depth [--min_DP 12]</summary>
<p></p>
<p>Minimum sequencing depth per individual to call heterozygous sites.</p>
</details>

<details><summary>Hardy-Weinberg equilbrium [--hwe_cutoff 0.01]</summary>
<p></p>
<p>Filter variants that are out of Hardy-Weinberg equilibrium (p-value threshold)</p>
</details>

<details><summary>Output directory [--outdir .]</summary>
<p></p>
<p>Specify output directory</p>
</details>


## Output

The pipeline outputs a single VCF-formated file containing the called and filtered genotypes for each distinct individual in the samples file. Each variant is annotated with the following extra information:

- **ID field:** dbSNP rs number
- **INFO/CAF:** 1000 genomes project allele frequency (from dbSNP annotation file)
- **INFO/TOPMED:** TOPMED project allele frequency (from dbSNP annotation file)
- **INFO/AA:** Inferred ancenstral allele from EPO/PECAN alignments (see "Input" for information about how this is obtained)
