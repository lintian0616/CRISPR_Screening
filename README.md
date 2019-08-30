CRISPR Screening (computational part)
================================

## Install Packages

We will use [pip](https://pypi.org/) to install three packages, **twobitreader** and **biopython**.

First, upgrade `pip` to the latest version.

```
sudo pip install --upgrade pip
```

Then, we will use `pip` to install the following packages.

```
sudo pip install twobitreader
sudo pip install biopython
```

For [seqmap](http://www-personal.umich.edu/~jianghui/seqmap/), install the version **1.0.13** source code for all platforms and compile with the command.

```
wget http://www-personal.umich.edu/~jianghui/seqmap/download/seqmap-1.0.13-src.zip

unzip seqmap-1.0.13-src.zip

cd seqmap-1.0.13-src/

g++ -O3 -m64 -o seqmap match.cpp
```

You should see the binary file `seqmap` in the folder. Place seqmap in the same folder as the Python script `design_library.py`.

## Prepare Gene List

Here is an example.


  Name   |  Chrom  |   Start  |    End   |
:-------:|:-------:|:--------:|:--------:|
  *EGFR* |   chr7  | 55086525 | 55086725 |
 *LPAR5* |  chr12  |  6745297 |
 *GPR35* |   chr2  | 241544625|
 
 

