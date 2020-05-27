# checker_complete.2.pl <!-- omit in toc -->

Finding outlier sequences in amino acid multiple sequence alignments.

## Table of Contents <!-- omit in toc -->

- [General description](#general-description)
- [Usage](#usage)
  - [Input](#input)
    - [Default mode](#default-mode)
      - [Checking multiple files](#checking-multiple-files)
    - [1KITE compatibility mode](#1kite-compatibility-mode)
  - [Content of the subjects file](#content-of-the-subjects-file)
    - [Default mode](#default-mode-1)
    - [1KITE compatibility mode](#1kite-compatibility-mode-1)
  - [Substitution matrix](#substitution-matrix)
  - [Output](#output)
- [Test data](#test-data)
- [List of changes](#list-of-changes)
- [Authors/Citation/Copyright](#authorscitationcopyright)
- [History](#history)

## General description

`checker_complete.2.pl` is a complete reimplementation of the outlier check routine provided by the `checker_complete.1.3.x.2.pl` scripts. It includes a number of changes, bug fixes, and improvements, such as the possibility to define arbitrary reference sets for the outlier distance calculation (for a detailed list of fixes and improvements, [see below](#list-of-changes)).

The script offers two different execution modes. The **default mode** allows the definition of arbitrary reference sequence sets and subject sequence sets, respectively. The **1KITE compatibility mode** defines reference and subject species based on the FASTA headers in the alignment (two pipe symbols `|` for reference species, three for subject species) and integrates flawlessly in the current 1KITE pipeline.

## Usage

How to use `checker_complete.2.pl` is probably best explained with a few examples.

Calling `checker_complete.2.pl` without any parameters (or `-h`), will show all available options:

```Text
$ ./checker_complete.2.pl

Usage: ./checker_complete.2.pl -a alignment(s) [further options]

Options:
      -a    Path to alignment file(s)
            [Required]

      -s    Path to file with taxa of interest to test

      -r    Path to reference taxa file

      -k    Ensure 1KITE compatible input and output

      -m    Substitution matrix [possible values: b62|b80|b62n|b80n]
            [Default: b62 (= original BLOSUM62 matrix)]

            WARNING: b62n and b80n are not tested.

      -l    Name of log file
            [Default: log.txt]

      -o1   Name of first output file
            [Default: outlier_1.txt]

            Contains for each alignment a list of outlier taxa and
            their distance score.

    -o2     Name of second output file
            [Default: outlier_2.txt]

            Contains for each alignment a list of outlier FASTA
            headers.

    -h      Print this help message

    -v      Print version
```

The only required input is one or more alignment file(s) that should be checked for outlier sequences (option: `-a`). Additionally, depending on the [execution mode](#input), you can provide a file to the script that lists the headers of all sequences you would like to use as reference taxa (option: `-r`) and/or a file that contains the headers of all sequences you would like to filter (*i.e.*, the subject file; option: `-s`). The content of the subject file is explained in detail [below](#content-of-the-subjects-file).

### Input

#### Default mode

The script only requires the name of one or several alignment files that should be checked for outliers. Reference sequences and subject sequences are then deduced automatically based on the name of the alignment file. For example, calling the script with an alignment file name `EOG090X002.aa.fa`

```bash
checker_complete.2.pl -a EOG090X002.aa.fa
```

assumes the presence of a file called `EOG090X002.aa.fa.refs` that contains the sequence headers (without the `>` symbol) of all taxa that should be used as reference species. All remaining sequences in the alignment will be treated as subject sequences and searched for outlier sequences.

You can, however, provide a **reference sequence file** explicitly (option: `-r`):

```bash
checker_complete.2.pl -a EOG090X002.aa.fa -r my.refs.txt
```

In this case, the sequences in `EOG090X002.aa.fa` that belong to the sequence headers provided in the file `my.refs.txt` (without the `>` symbol) will be used as reference species. All remaining sequences in the alignment will be treated as subject sequences and searched for outlier sequences.

You can modify this behavior by providing a list of FASTA headers of those **subject sequences** (option: `-s`) you would like to test:

```bash
checker_complete.2.pl -a EOG090X002.aa.fa -r my.refs.txt -s subjects.txt
```

Here, the sequences in `EOG090X002.aa.fa` that belong to the sequence headers provided in the file `my.refs.txt` will be used as reference species and all sequences with their FASTA headers listed in `subjects.txt` that are present in the alignment file (`EOG090X002.aa.fa`) will be checked for outliers.

##### Checking multiple files

`checker_complete.2.pl` supports the use of **wildcards** to check multiple alignment files in the same run:

```bash
checker_complete.2.pl -a EOG*.aa.fa
```

This will check all alignments found in files that match the pattern `EOG*.aa.fa` for outlier sequences. For each alignment file, the script requires the presence of a file called `EOG*.aa.fa.refs` that contains the sequence headers of all taxa which should be used as reference species for this specific alignment. For example, given two alignment files called `EOGxxx.aa.fa` and `EOGyyy.aa.fa`, the script requires the presence of two additional files `EOGxxx.aa.fa.refs` containing the reference sequence headers for `EOGxxx.aa.fa` and `EOGyyy.aa.fa.refs` containing the reference sequence headers for `EOGyyy.aa`.fa, respectively. All other sequences in the respective alignment will be treated as subject sequences and searched for outlier sequences.

If you want to check multiple alignment files and the **FASTA headers of all reference taxa are identical between alignment files**, you can also provide a single reference file that contains these headers (option: `-r`):

```bash
checker_complete.2.pl -a EOG*.aa.fa -r my.universal.refs.txt
```

All other sequences in the respective alignment will be treated as subject sequences and searched for outlier sequences.

In **default mode**, the following algorithmic decisions are enforced:

1. The **reference species** file should contain the full FASTA headers (without the `>` symbol) of the sequences in the alignment that should be treated as reference species.

2. If a **subjects file** is provided, it is assumed that it contains the full FASTA header (without the `>` symbol) of the sequences that should be checked for outliers.
If **no subject file** is provided, all sequences in the alignment that are not present in the reference species list are used as subject sequences.

3. The BLOSUM distances ([see below](#substitution-matrix)) of the amino-acid sequences of **all reference species** are calculated. Based on these distances, a cut-off value is derived. After that, for each transcript (subject) sequence the **minimal BLOSUM distance to any of the reference taxa** is calculated. Sequences with a minimal BLOSUM distance to their closest reference species larger than the cut-off value are classified as outliers.

4. The **output** format of the first output file (option: `-o1`) will be different to the format originally provided by `checker_complete.1.3.1.2.pl`. The content, however, is similar.

#### 1KITE compatibility mode

This mode is enforced by using the `-k` flag:

```bash
checker_complete.2.pl -a EOG090X002.aa.fa -s subjects.txt -k
```

This will change the algorithm in several ways:

1. A **subject file** that now only contains the name of all transcriptome taxa as written in the third field ([see below](#1kite-compatibility-mode-1)) is now mandatory.

2. **Reference species** in the alignment are automatically identified based on the number of pipe symbols (`|`) in their sequence header (*i.e.*, two pipe symbols, see next point).

3. The **FASTA headers** of the sequences in the alignment are required to be in `HaMStrad` format, *i.e.* reference species headers have two pipe symbols (*e.g.*, `>EOGxyz|reference_taxon_name|GenID`) and all subject (transcriptome) sequence headers have three pipe symbols (*e.g.*, `>EOGxyz|reference_taxon_name|subject_taxon_name|contig-header`).

4. As in default mode, the BLOSUM distances of the amino-acid sequences of all reference species are calculated. Based on these distances, the cut-off value is determined. After that, **the BLOSUM distance of each transcript** (subject) **sequence to the corresponding sequence of the reference taxon found as best reciprocal hit during orthology assignment** (`reference_taxon_name`) is calculated. Sequences with a minimal BLOSUM distance to their reference species larger than the cut-off value are classified as outliers.

5. The **output** format of the first output file (option: `-o1`) will be compatible to the format provided by the original `checker_complete.1.3.1.2.pl` script. This will ensure that existing scripts in the 1KITE community that make use of this file for further downstream analyses don't need to be modified.

Requesting the 1KITE compatibility mode and providing a reference species list simultaneously, will result in an error:

```bash
$ checker_complete.2.pl -a EOG090X002.aa.fa -s subjects.txt -k -r my.refs
[ERROR] Contradictory settings. 1KITE-compatibility requested and reference file given.
```

### Content of the subjects file

#### Default mode

In default mode, the subjects file is optional. If it is provided, however, it must contain the full FASTA header (without the `>` symbol) of the sequences that should be checked for outliers. For example, given the following FASTA headers in the alignment file:

```Text
>EOG508HX6|BMORI_3.2|Tricholepidion_gertschi|C20406_l5235
>EOG508HX6|BMORI_3.2|Drosophila_melanogaster|C20406_l6733
...
```

The content of the subjects file must be:

```Text
EOG508HX6|BMORI_3.2|Tricholepidion_gertschi|C20406_l5235
EOG508HX6|BMORI_3.2|Drosophila_melanogaster|C20406_l6733
```

If no subjects file is provided, all sequences in the alignment that are not present in the reference species list are used as subject sequences.

#### 1KITE compatibility mode

In this mode, the subjects file must only contain the name of all transcriptome taxa as written in the third field of their sequence headers, one taxon name per line. For example, given the following sequence headers in the alignment file:

```Text
>EOG508HX6|BMORI_3.2|Tricholepidion_gertschi|C20406_l5235
>EOG508HX6|BMORI_3.2|Drosophila_melanogaster|C20406_l6733
...
```

The content of the subjects file must be:

```Text
Tricholepidion_gertschi
Drosophila_melanogaster
```

### Substitution matrix

As in previous instances of this script, `checker_complete.2.pl` uses the BLOSUM62 log-odds matrix for the sequence distance calculation. However, the BLOSUM62 matrix might not be appropriate for very similar sequences. Thus, the BLOSUM80 matrix is also supported now and can be requested with the `-m` option:

```bash
checker_complete.2.pl -a EOG090X002.aa.fa -s subjects.txt -m b80
```

Note that the script also supports the [updated versions](https://www.nature.com/articles/nbt0308-274) of the original BLOSUM62 and BLOSUM80 matrices, respectively (`b62n` and `b80n`, respectively). However, their use is currently experimental and may not be include in future versions of this script.

### Output

The original `checker_complete.1.3.1.2.pl` and `checker_complete.1.3.2.2.pl` scripts differed only in their output. Both created identical (and therefore redundant) log files (`log_1.txt` and `log_2.txt`, respectively) that contained various details about the analyzed alignment:

- The name of the alignment file,
- the alignment length,
- the missing subject sequences,
- the total number of taxa in the alignment,
- the number of transcriptome (subject) sequences,
- the sequence overlap length,
- the minimum and maximum distance values,
- the median, the lower (Q1) and upper (Q3) quartile,
- the outlier cutoff,
- the BLOSUM62 distance and relative BLOSUM62 distance to the closest reference taxon (during orthology assignment step),
- whether a sequence was classified as outlier,
- the number of outliers in the alignment,
- the total number of checked alignments, and
- the total number of files with outliers.

Additionally, both scripts created an outlier file. The outlier file created by `checker_complete.1.3.1.2.pl` (called `outlier_1.txt`) contained a reduced set of information:

- The name of the alignment file,
- the minimum and maximum distance values,
- the median, the lower (Q1) and upper (Q3) quartile,
- the BLOSUM62 distance and relative BLOSUM62 distance to the closest reference taxon (during orthology assignment step), and
- the name of the sequences that were classified as outlier.

The outlier file created by `checker_complete.1.3.2.2.pl` (called `outlier_2.txt`) contained only:

- The name of the alignment file and
- the complete sequence headers of those sequences that were classified as outlier.

`checker_complete.2.pl` now creates all three files (a single log file and both outlier files) simultaneously. To ensure compatibility to existing 1KITE workflows, the default names of these files are `log.txt`, `outlier_1.txt`, and `outlier_2.txt`. All names, however, can be modified by the user:

```bash
checker_complete.2.pl -a file.fas -l my.log -o1 outlier1.log -o2 outlier2.log
```

## Test data

To test whether checker_complete.2.pl works as intended, some examples files are
provided in the `example` directory. To analyse the test data, simply change into the directory and execute the `runMe.sh` script:

```bash
cd example
chmod +x runMe.sh
./runMe.sh
```

## List of changes

In the following, a list of major changes and resolved bugs compared to  `checker_complete.1.3.x.2.pl` is given.

**Changes/Improvements:**

- Arbitrary sets for reference taxa and taxa to test can now be provided by the user
- BLOSUM80 can now be used for distance calculation
- Outlier calculation is now based on the interquartile range (IQR)
- File names for log files and outlier files can now be selected by the user
- Output files of `checker_complete.1.3.1.2.pl` and `checker_complete.1.3.2.2.pl` are now created simultaneously
- Output files are not overwritten

**Bugs fixed:**

- Now works with newer PERL versions
- Quartile calculation was inconsistent
- Not all reference taxa were included in the output
- Lowest distance between reference taxa was not calculated correct

## Authors/Citation/Copyright

`checker_complete.2.pl` is a reimplementation of the two scripts `checker_complete.1.3.1.2.pl` and `checker_complete.1.3.2.2.pl` originally written by Bernhard Misof (zmb, ZFMK) and published under the GNU General Public License version 2 (or any later version).

Both scripts have been published in:

> Misof B, Liu S, Meusemann K, Peters RS, Donath A, Mayer C, et al. Phylogenomics resolves the timing and pattern of insect evolution. Science. 2014;346:763–7.

`checker_complete.2.pl` is published under the GPLv3.

## History

v 2.0.0b - First release.
