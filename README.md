CRISPR Screening (computational part)
================================

This is a tutorial adapted from [Nature Protocol paper](https://www.nature.com/articles/nprot.2017.016).

## Install Packages

We will use [pip](https://pypi.org/) to install three packages, **twobitreader** and **biopython**.

If `pip` is not installed, we will update packages and install `pip` first.

```
sudo apt update
sudo apt install python-pip
pip --version
```

If `pip` is installed, upgrade `pip` to the latest version.

```
sudo pip install --upgrade pip
```

Then, we will use `pip` to install the following packages.

```
sudo pip install twobitreader
sudo pip install biopython
```

We also need to install [dryscrape](https://dryscrape.readthedocs.io/en/latest/) package.

```
sudo apt-get install qt5-default libqt5webkit5-dev build-essential python-lxml xvfb

sudo pip install dryscrape
```

For [seqmap](http://www-personal.umich.edu/~jianghui/seqmap/), install the version **1.0.13** source code for all platforms and compile with the command.

```
wget http://www-personal.umich.edu/~jianghui/seqmap/download/seqmap-1.0.13-src.zip

unzip seqmap-1.0.13-src.zip
```

## Generate Library

### a library targeting a custom set of genomic coordinates

* Install **seqmap**

Place seqmap in the same folder as the Python script `design_library.py`.

```
cd seqmap-1.0.13-src/

g++ -O3 -m64 -o seqmap match.cpp
```

We will add file `seqmap` to executable path in `.bashrc` file.

```
export PATH=$PATH:/home/lintian0616/seqmap-1.0.13-src
```

* Prepare target gene list

Prepare Gene List below in a *csv* format. The **target_genes.csv** file should contain the headers **name**, **chrom**, **start**, and **end**.

  name   |  chrom  |   start  |    end   |
:-------:|:-------:|:--------:|:--------:|
  *EGFR* |   chr7  | 55086525 | 55086725 |
 *LPAR5* |  chr12  |  6745297 |  6745497 |
 *GPR35* |   chr2  | 241544625| 241544825|
 
* run **design_library.py**

Download genome 2bit file (e.g. [hg38.2bit](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit)) to construct a database of off-target scores based on the position and distribution of mismatches between each spacer sequence and similar sequences in the genome.

For each region in the target genes *csv* file, the **design_library.py** Python script will identify potential sgRNAs and select a specified number of sgRNAs with fewer potential off-target sites using this database for the custom library.
 
Here is the example command:

```
python design_library.py -o final_guides.csv -i hg38 -g target_genes.csv -gc 25 -s 20 -n 3 -gecko
```

**1.** `-o`: Output **csv** file with names for target genes, spacer sequences, spacer orientations, chromosome locations, cleavage site locations, off-target scores, and oligo library sequences in columns from left to right (default: **final_guides.csv**)

**2.** `-i`: Prefix of input genome 2-bit file (default: **hg19**)

**3.** `-g`: Target-gene **csv** file

**4.** `-gc`: Minimum GC content required for an sgRNA spacer sequence (default: **25**)

**5.** `-s`: Minimum spacing required between cleavage sites of sgRNAs targeting the same genomic region (default: **20**)

**6.** `-n`: Maximum number of guides selected targeting each gene in the target-gene **csv** file (default: **3**)

**7.** `-gecko`: add flanking sequences to the spacers for the oligo library synthesis

### a library an existing library

* Prepare a `csv` file containing the names of the target genes, with each line corresponding to one gene.

* Prepare another `csv` file for the annotated genome-scale library, with the names of each gene in the first column and the respective spacer sequences in the second column. Each line contains a different spacer sequence. The gene names in the target genes file should be in the same format as the names of the annotated library file.

* run **design_targeted_library.py**

Here is the example command:

```
python design_targeted_library.py -o selected_sgRNAs_from_library.csv -l sabatini_library_targets.csv -g target_genes_2.csv -gecko
```

**1.** `-o`: Output **csv** file with names for target genes, corresponding spacer sequences, and oligo library sequences in columns from left to right (default: **oligos.csv**)

**2.** `-l`: Annotated library **csv** file with names in the first column and corresponding spacer sequences in the second column (default: **annotated_library.csv**)

**3.** `-g`: Target-gene **csv** file

**4.** `-gecko`: add flanking sequences to the spacers for the oligo library synthesis

## amplified oligo structure

![oligo structure](./Examples/oligo_structure.jpg)
