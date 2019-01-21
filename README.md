# candidates.py
A Python script that scrapes the output files of deepsig, effectorP, Interproscan (including SignalP) and dbCAN (cazymes) and selects effector candidates on the basis of user selectable criteria, i.e. effectorP score, number of cysteins, molecular weight of the mature protein.

# Usage

usage: 
```
candidates.py [-h] [-i INPUT] [-o OUTPUT] [-e EFFECTORPSCORE] [-m MWCUTOFF] [-c CYSTEINS]
```
```
optional arguments:

-h, --help 
show this help message and exit

-i INPUT, --input INPUT
path to a folder containing another folder with input
files (default is a folder called "input" inside the
directory where the script is)

-o OUTPUT, --output OUTPUT
folder for output files (default "output" inside the directory where the script is)

- -e EFFECTORPSCORE, --effectorPscore EFFECTORPSCORE
set minimum effectoP score (default = 0.8)

-m MWCUTOFF, --MWcutoff MWCUTOFF
set maximum molecular weight of mature candidate (default = 25)

-c CYSTEINS, --cysteins CYSTEINS
set minimum number of cysteins in mature candidate (default = 2)
```
