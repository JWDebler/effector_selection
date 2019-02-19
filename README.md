# candidates
A Python 3 script that scrapes the output files of [deepsig](https://deepsig.biocomp.unibo.it/welcome/default/index), [effectorP](http://effectorp.csiro.au/), [Interproscan](https://www.ebi.ac.uk/interpro/search/sequence-search) (including [SignalP](http://www.cbs.dtu.dk/services/SignalP/)) and [dbCAN](http://cys.bios.niu.edu/dbCAN2/blast.php) (cazymes) and selects effector candidates on the basis of user selectable criteria, i.e. effectorP score, number of cysteins and molecular weight of the mature protein.

# Prerequisits

This script requires an input folder per sample/isolate. The foldername will be used as the sampleID. This folder needs to contain following files (which need to at least contain the following as suffix i.e. `yourname.proteins.fasta` will work, as long as it contains the string `proteins.fasta`):
- proteins.fasta (proteome file used as input for all the other programs)
- deepsig.out (deepsig output file)
- effectorP.tsv (effectorP output file)
- interproscan.tsv (interproscan output file)
- cazymes.txt (OPTIONAL - output of dbCAN)

# Usage

usage: 
```
candidates.py [-h] [-i INPUT] [-o OUTPUT] [-e EFFECTORPSCORE] [-m MWCUTOFF] [-c CYSTEINS]
```
example:
```
candidates.py -i /path/inputfolder -o /path/outputfolder -e 0.8 -m 25 -c 2
```
optional arguments:
```
-h, --help 
show this help message and exit

-i INPUT, --input INPUT
path to a folder containing another folder with input
files (default is a folder called "input" inside the
directory where the script is)

-o OUTPUT, --output OUTPUT
folder for output files (default "output" inside the directory where the script is)

-e EFFECTORPSCORE, --effectorPscore EFFECTORPSCORE
set minimum effectoP score (default = 0.8)

-m MWCUTOFF, --MWcutoff MWCUTOFF
set maximum molecular weight of mature candidate (default = 25)

-c CYSTEINS, --cysteins CYSTEINS
set minimum number of cysteins in mature candidate (default = 2)
```
