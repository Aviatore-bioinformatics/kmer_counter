# Environment preparation
## Install miniconda
```
conda create -n <environment-name> python=3.7
```
```
conda activate -n <environment-name>
```
## Install required software using conda
### Install jellyfish
Jellyfish is a k-mer counter based on a multi-threaded hash table implementation.
```
conda install -c bioconda jellyfish
```
### Install MEME Suit
[MEME Suit](http://web.mit.edu/meme_v4.11.4/share/doc/overview.html) - Motif-based sequence analysis tools. Kmer-counter uses two tools from this tools set:
* [iupac2meme](http://web.mit.edu/meme_v4.11.4/share/doc/iupac2meme.html) - converts an IUPAC string (in this case DNA sequences) to MEME format.
* [tomtom](http://web.mit.edu/meme_v4.11.4/share/doc/tomtom.html) - compares one or more motifs (in MEME format) against a database of known motifs (e.g., JASPAR) (also in MEME format). Below you can find links to some motifs databases:
  * [meme](https://meme-suite.org/meme/db/motifs)
  * [JASPAR](http://jaspar2018.genereg.net/downloads/)
> The motif database should be represented as a **single** batch file which localisation must be provided in the ***config.txt*** file.

To install the MEME Suit write the following command:

```
conda install -c bioconda meme
```
 ## Install Python modules
 To install all required Python modules just run below command:
 ```
 pip install -r requirements.txt
 ```
 The *requirements.txt* file is provided together with the Kmer-counter. Below are listed all required Python modules:
 * certifi
 * intervaltree
 * numpy
 * pandas
 * patsy
 * python-dateutil
 * pytz
 * scipy
 * six
 * sortedcontainers
 * statsmodels