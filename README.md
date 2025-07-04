# WisecondorX

![Conda Version](https://img.shields.io/conda/v/bioconda/wisecondorx)
[![GitHub Actions CI Status](https://github.com/matthdsm/wisecondorx/actions/workflows/ci.yml/badge.svg)](https://github.com/matthdsm/wisecondorx/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/matthdsm/wisecondorx/actions/workflows/linting.yml/badge.svg)](https://github.com/matthdsm/wisecondorx/actions/workflows/linting.yml)

## 1. Background

After extensively comparing different (shallow) whole-genome sequencing-based copy number detection tools,
including [WISECONDOR](https://github.com/VUmcCGP/wisecondor), [QDNAseq](https://github.com/ccagc/QDNAseq),
[CNVkit](https://github.com/etal/cnvkit), [Control-FREEC](https://github.com/BoevaLab/FREEC),
[BIC-seq2](http://compbio.med.harvard.edu/BIC-seq/) and
[cn.MOPS](https://bioconductor.org/packages/release/bioc/html/cn.mops.html),
WISECONDOR appeared to normalize sequencing data in the most consistent way, as shown by
[our paper](https://www.ncbi.nlm.nih.gov/pubmed/30566647). Nevertheless, WISECONDOR has limitations:
Stouffer's Z-score approach is error-prone when dealing with large amounts of aberrations, the algorithm
is extremely slow (24h) when operating at small bin sizes (15 kb), and sex chromosomes are not part of the analysis.
Here, we present WisecondorX, an evolved WISECONDOR that aims at dealing with previous difficulties, resulting
in overall superior results and significantly lower computing times, allowing daily diagnostic use. WisecondorX is
meant to be applicable not only to NIPT, but also gDNA, PGT, FFPE, LQB, ... etc.

## 2. Manual

### 2.1 Data Alignment and Post-Processing

We found superior results through WisecondorX when using [bowtie2](https://github.com/BenLangmead/bowtie2) with the `--local` and `--fast-local` options as alignment tool.
Post alignment, you can use a duplicate marking tool of choice.
Note that it is important that **no** read quality filtering is executed prior to running WisecondorX: this software
requires low-quality reads to distinguish informative bins from non-informative ones.

### 2.2 WisecondorX

#### 2.2.1 Installation

WisecondorX is available through [Conda](https://conda.io/docs/), which install the binary and all necessary dependencies.

```bash
conda install -f -c conda-forge -c bioconda wisecondorx
```

#### 2.2.2 Building from source

WisecondorX is written in Go. To build from source, you need to have Go installed on your system.
```bash
git clone https://github.com/matthdsm/wisecondorx.git
cd wisecondorx
go build -o wisecondorx .
```

Note: When building from source, you need to manually install the `samtools` requirement.

#### 2.2.3 Running WisecondorX

There are three main stages (converting, reference creating and predicting) when using WisecondorX:

- Convert aligned reads to .npz files (for both reference and test samples)
- Create a reference (using reference .npz files)
  - **Important notes**
    - Automated sex prediction, required to consistently analyze sex chromosomes, is based on a Gaussian mixture
      model. If few samples (<20) are included during reference creation, or not both male and female samples (for
      NIPT, this means male and female feti) are represented, this process might not be accurate. Therefore,
      alternatively, one can manually tweak the `--yfrac` parameter.
    - It is of paramount importance that the reference set consists of exclusively negative control samples that
      originate from the same sequencer, mapper, reference genome, type of material, ... etc, as the test samples.
      As a rule of thumb, think of all laboratory and _in silico_ steps: the more sources of bias that can be omitted,
      the better.
    - Try to include at least 50 samples per reference. The more the better, yet, from 500 on it is unlikely to
      observe additional improvement concerning normalization.
- Predict copy number alterations (using the reference file and test .npz cases of interest)

An overview of the WisecondorX commands can be found in the [CLI documentation](./cli.md).
Default parameters are optimized for shallow whole-genome sequencing data (0.1x - 1x coverage) and reference bin sizes ranging from 50 to 500 kb.


## 3. Algorithm

To understand the underlying algorithm, we highly recommend reading
[Straver et al (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24170809); and its update shortly introduced in
[Huijsdens-van Amsterdam et al (2018)](https://www.nature.com/articles/gim201832.epdf). Numerous adaptations to this
algorithm have been made, yet the central normalization principles remain. Changes include e.g. the inclusion of a sex
prediction algorithm, sex handling prior to normalization (ultimately enabling X and Y predictions), between-sample
Z-scoring, inclusion of a weighted circular binary segmentation algorithm and improved codes for obtaining tables and
plots.


## 4. Outputs

**\<prefix\>_bins.bed**

This file contains all bin-wise information. When data is 'NaN', the corresponding bin is included in the blacklist.
The Z-scores are calculated as default using the within-sample reference bins as a null set.

**\<prefix\>_segments.bed**

This file contains all segment-wise information. A combined Z-score is calculated using a between-sample Z-scoring
technique (the test case vs the reference cases).

**\<prefix\>_aberrations.bed**

This file contains aberrant segments, defined by the `--beta` or `--zscore` parameters.

**\<prefix\>_statistics.bed**

This file contains some interesting statistics (per chromosome and overall). The definition of the Z-scores matches the one from
the `<prefix>_segments.bed`. Particularly interesting for NIPT.

## 5. Dependencies

WisecondorX depends on `samtools` to read input bam/cram files.
