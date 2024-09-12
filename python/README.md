# Background

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

# Manual

## Mapping

We found superior results through WisecondorX when using [bowtie2](https://github.com/BenLangmead/bowtie2) as a mapper.
Note that it is important that **no** read quality filtering is executed prior to running WisecondorX: this software
requires low-quality reads to distinguish informative bins from non-informative ones.

## WisecondorX

### Installation

Stable releases can be installed through pip install. This option ascertains the latest version is
downloaded, however, it does not install R [dependencies](#dependencies).

```bash

pip install -U git+https://github.com/CenterForMedicalGeneticsGhent/WisecondorX
```

Alternatively, [Conda](https://conda.io/docs/) additionally installs all necessary [depedencies](#dependencies),
however, the latest version might not be downloaded.

```bash

conda install -f -c conda-forge -c bioconda wisecondorx
```

### Running WisecondorX

There are three main stages (converting, reference creating and predicting) when using WisecondorX:

- Convert aligned reads to .npz files (for both reference and test samples)
- Create a reference (using reference .npz files)
  - **Important notes**
    - Automated gender prediction, required to consistently analyze sex chromosomes, is based on a Gaussian mixture
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

### Stage (1) Convert aligned reads (bam/cram) to .npz

```bash

WisecondorX convert input.bam/cram output.npz [--optional arguments]
```

| <br>Optional argument <br><br> | Function                                                                                                                                                                                                 |
| :----------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--reference x`                | Fasta reference to be used with cram inputs                                                                                                                                                              |
| `--binsize x`                  | Size per bin in bp; the reference bin size should be a multiple of this value. Note that this parameter does not impact the resolution, yet it can be used to optimize processing speed (default: x=5e3) |
| `--normdup`                    | Use this flag to avoid duplicate removal                                                                                                                                                                 |

&rarr; Bash recipe at `docs/include/pipeline/convert.sh`

### Stage (2) Create reference

```bash

WisecondorX newref reference_input_dir/*.npz reference_output.npz [--optional arguments]
```

| <br>Optional argument <br><br> | Function                                                                                                                                    |
| :----------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------ |
| `--nipt`                       | **Always include this flag for the generation of a NIPT reference**                                                                         |
| `--binsize x`                  | Size per bin in bp, defines the resolution of the output (default: x=1e5)                                                                   |
| `--refsize x`                  | Amount of reference locations per target; should generally not be tweaked (default: x=300)                                                  |
| `--yfrac x`                    | Y read fraction cutoff, in order to manually define gender. Setting this to 1 will treat all samples as female                              |
| `--plotyfrac x`                | plots Y read fraction histogram and Gaussian mixture fit to file x, can help when setting `--yfrac` manually; software quits after plotting |
| `--cpus x`                     | Number of threads requested (default: x=1)                                                                                                  |

&rarr; Bash recipe at `docs/include/pipeline/newref.sh`

### Stage (3) Predict copy number alterations

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
| `--gender x`                   | Force WisecondorX to analyze this case as male (M) or female (F). Useful when e.g. dealing with a loss of chromosome Y, which causes erroneous gender predictions (choices: x=F or x=M)                                           |
| `--bed`                        | Outputs tab-delimited .bed files (trisomy 21 NIPT example at `./example.bed`), containing all necessary information **(\*)**                                                                                                      |
| `--plot`                       | Outputs custom .png plots (trisomy 21 NIPT example at `./example.plot`), directly interpretable **(\*)**                                                                                                                          |
| `--ylim [a,b]`                 | Force WisecondorX to use y-axis interval [a,b] during plotting, e.g. [-2,2]                                                                                                                                                       |
| `--cairo`                      | Some operating systems require the cairo bitmap type to write plots                                                                                                                                                               |
| `--seed`| Random seed for segmentation algorithm (default:None)                                                                                                                                                                                  |

<sup>**(\*)** At least one of these output formats should be selected</sup>

&rarr; Bash recipe at `docs/include/pipeline/predict.sh`

### Additional functionality

```bash

WisecondorX gender test_input.npz reference_input.npz
```

Returns gender according to the reference.

# Parameters

The default parameters are optimized for shallow whole-genome sequencing data (0.1x - 1x coverage) and reference bin
sizes ranging from 50 to 500 kb.

# Underlying algorithm

To understand the underlying algorithm, I highly recommend reading
[Straver et al (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24170809); and its update shortly introduced in
[Huijsdens-van Amsterdam et al (2018)](https://www.nature.com/articles/gim201832.epdf). Numerous adaptations to this
algorithm have been made, yet the central normalization principles remain. Changes include e.g. the inclusion of a gender
prediction algorithm, gender handling prior to normalization (ultimately enabling X and Y predictions), between-sample
Z-scoring, inclusion of a weighted circular binary segmentation algorithm and improved codes for obtaining tables and
plots.

# Interpretation results

## Plots

Every dot represents a bin. The dots range across the X-axis from chromosome 1 to X (or Y, in case of a male). The
vertical position of a dot represents the ratio between the observed number of reads and the expected number of reads,
the latter being the 'normal' state. As these values are log2-transformed, copy neutral dots should be close-to 0.
Importantly, notice that the dots are always subject to Gaussian noise. Therefore, segments, indicated by horizontal
white lines, cover bins of predicted equal copy number. The size of the dots represents the variability at the reference
set. Thus, the size increases with the certainty of an observation. The same goes for the line width of the segments.
Vertical grey bars represent the blacklist, which matches mostly hypervariable loci and repeats. Finally, the horizontal
colored dotted lines show where the constitutional 1n and 3n states are expected (when constitutional DNA is at 100%
purity). Often, an aberration does not reach these limits, which has several potential causes: depending on your type
of analysis, the sample could be subject to tumor fraction, fetal fraction, a mosaicism, ... etc. Sometimes, the
segments do surpass these limits: here it's likely you are dealing with 0n, 4n, 5n, 6n, ...

## Tables

### ID_bins.bed

This file contains all bin-wise information. When data is 'NaN', the corresponding bin is included in the blacklist.
The Z-scores are calculated as default using the within-sample reference bins as a null set.

### ID_segments.bed

This file contains all segment-wise information. A combined Z-score is calculated using a between-sample Z-scoring
technique (the test case vs the reference cases).

### ID_aberrations.bed

This file contains aberrant segments, defined by the [`--beta`](#stage-3-predict-copy-number-alterations) or
[`--zscore`](#stage-3-predict-copy-number-alterations) parameters.

### ID_statistics.bed

This file contains some interesting statistics (per chromosome and overall). The definition of the Z-scores matches the one from
the 'ID_segments.bed'. Particularly interesting for NIPT.

# Dependencies

- R (v4.3.3) packages
  - jsonlite (v1.8.8)
- R Bioconductor packages
  - DNAcopy (v1.76.0)
- Python (> v3.6) libraries
  - scipy (v1.13.0)
  - scikit-learn (v1.4.2)
  - pysam (v0.22.0)
  - numpy (v1.26.4)
  - matplotlib (v3.8.4)

And of course, other versions are very likely to work as well.
