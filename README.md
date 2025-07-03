# 1. Background

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

# 2. Manual

## 2.1 Data Alignment and Post-Processing

We found superior results through WisecondorX when using [bowtie2](https://github.com/BenLangmead/bowtie2) with the `--local` and `--fast-local` options as alignment tool.
Post alignment, you can use a duplicate marking tool of choice.
Note that it is important that **no** read quality filtering is executed prior to running WisecondorX: this software
requires low-quality reads to distinguish informative bins from non-informative ones.

## 2.2 WisecondorX

### 2.2.1 Installation

WisecondorX is available through [Conda](https://conda.io/docs/), which install the binary and all necessary dependencies.

```bash
conda install -f -c conda-forge -c bioconda wisecondorx
```

### 2.2.2 Building from source

WisecondorX is written in Go. To build from source, you need to have Go installed on your system.
```bash
git clone https://github.com/matthdsm/wisecondorx.git
cd wisecondorx
go build -o wisecondorx .
```

Note: When building from source, you need to manually install the `samtools` requirement.

### 2.2.3 Running WisecondorX

There are three main stages (converting, reference creating and predicting) when using WisecondorX:

- Convert aligned reads to .npz files (for both reference and test samples)
- Create a reference (using reference .npz files)
  - **Important notes**
    - Automated sex prediction, required to consistently analyze sex chromosomes, is based on a Gaussian mixture
      model. If few samples (<20) are included during reference creation, or not both male and female samples (for
      NIPT, this means male and female feti) are represented, this process might not be accurate. Therefore,
      alternatively, one can manually tweak the [`--yfrac`](#stage-2-create-reference) parameter.
    - It is of paramount importance that the reference set consists of exclusively negative control samples that
      originate from the same sequencer, mapper, reference genome, type of material, ... etc, as the test samples.
      As a rule of thumb, think of all laboratory and _in silico_ steps: the more sources of bias that can be omitted,
      the better.
    - Try to include at least 50 samples per reference. The more the better, yet, from 500 on it is unlikely to
      observe additional improvement concerning normalization.
- Predict copy number alterations (using the reference file and test .npz cases of interest)

### 2.2.4 Convert aligned reads (bam/cram) to .npz

```bash
wisecondorx convert input.bam/cram output.npz [--optional arguments]
```

| <br>Optional argument <br><br> | Function                                                                                                                                                                                                 |
| :----------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--reference x`                | Fasta reference to be used with cram inputs                                                                                                                                                              |
| `--binsize x`                  | Size per bin in bp; the reference bin size should be a multiple of this value. Note that this parameter does not impact the resolution, yet it can be used to optimize processing speed (default: x=5e3) |
| `--no-remove-duplicates`                    | Use this flag to avoid duplicate removal                                                                                                                                                                 |
| `--exclude-contigs` | Glob pattern to exclude certain contigs from conversion (default: "{*_alt,*_decoy,__random,chrUn_,HLA*,chrM,chrEBV}"|
| `--gonosomes x`                | Gonosome chromosomes to be used in the analysis; should generally not be tweaked (default: x=chrX, chrY)                                                                                                 |

### 2.2.5. Create reference

```bash
WisecondorX newref reference_input_dir/*.npz reference_output.npz [--optional arguments]
```

| <br>Optional argument <br><br> | Function                                                                                                                                    |
| :----------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------ |
| `--nipt`                       | **Always include this flag for the generation of a NIPT reference**                                                                         |
| `--binsize x`                  | Size per bin in bp, defines the resolution of the output (default: x=1e5)                                                                   |
| `--refsize x`                  | Amount of reference locations per target; should generally not be tweaked (default: x=300)                                                  |
| `--yfrac x`                    | Y read fraction cutoff, in order to manually define sex. Setting this to 1 will treat all samples as female                              |
| `--plotyfrac x`                | plots Y read fraction histogram and Gaussian mixture fit to file x, can help when setting `--yfrac` manually; software quits after plotting |
| `--cpus x`                     | Number of threads requested (default: x=1)                                                                                                  |

### 2.2.6. Predict copy number alterations

```bash
WisecondorX predict test_input.npz reference_input.npz output_id [--optional arguments]
```

| <br>Optional argument <br><br> | Function                                                                                                                                                                                                                          |
| :----------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--minrefbins x`               | Minimum amount of sensible reference bins per target bin; should generally not be tweaked (default: x=150)                                                                                                                        |
| `--maskrepeats x`              | Bins with distances > mean + sd \* 3 in the reference will be masked. This parameter represents the number of masking cycles and defines the stringency of the blacklist (default: x=5)                                           |
| `--zscore x`                   | Z-score cutoff to call segments as aberrations (default: x=5)                                                                                                                                                                     |
| `--alpha x`                    | P-value cutoff for calling circular binary segmentation breakpoints (default: x=1e-4)                                                                                                                                             |
| `--beta x`                     | When beta is given, `--zscore` is ignored. Beta sets a ratio cutoff for aberration calling. It's a number between 0 (liberal) and 1 (conservative) and, when used, is optimally close to the purity (e.g. fetal/tumor fraction)   |
| `--blacklist x`                | Blacklist for masking additional regions; requires headerless .bed file. This is particularly useful when the reference set is too small to recognize some obvious loci (such as centromeres; examples at `./example.blacklist/`) |
| `--sex x`                   | Force WisecondorX to analyze this case as male (M) or female (F). Useful when e.g. dealing with a loss of chromosome Y, which causes erroneous sex predictions (choices: x=F or x=M)                                           |
| `--seed`| Random seed for segmentation algorithm (default:None)                                                                                                                                                                                  |

# 2. Parameters

The default parameters are optimized for shallow whole-genome sequencing data (0.1x - 1x coverage) and reference bin
sizes ranging from 50 to 500 kb.

# 3. Algorithm

To understand the underlying algorithm, I highly recommend reading
[Straver et al (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24170809); and its update shortly introduced in
[Huijsdens-van Amsterdam et al (2018)](https://www.nature.com/articles/gim201832.epdf). Numerous adaptations to this
algorithm have been made, yet the central normalization principles remain. Changes include e.g. the inclusion of a sex
prediction algorithm, sex handling prior to normalization (ultimately enabling X and Y predictions), between-sample
Z-scoring, inclusion of a weighted circular binary segmentation algorithm and improved codes for obtaining tables and
plots.

# 4. Outputs

## \<prefix\>_bins.bed

This file contains all bin-wise information. When data is 'NaN', the corresponding bin is included in the blacklist.
The Z-scores are calculated as default using the within-sample reference bins as a null set.

## \<prefix\>_segments.bed

This file contains all segment-wise information. A combined Z-score is calculated using a between-sample Z-scoring
technique (the test case vs the reference cases).

## \<prefix\>_aberrations.bed

This file contains aberrant segments, defined by the `--beta` or
[`--zscore`] parameters.

## \<prefix\>_statistics.bed

This file contains some interesting statistics (per chromosome and overall). The definition of the Z-scores matches the one from
the `<prefix>_segments.bed`. Particularly interesting for NIPT.

# 5. Dependencies

WisecondorX depends on `samtools` to read input bam/cram files.
