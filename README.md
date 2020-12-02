# chimeramutagenesis

## SETUP
### Prepare python environment
Requires python3 with numpy installed. In my case, installed miniconda to create an environment that included `python 3.8`, `numpy` and `pandas` following https://hpc.nih.gov/apps/python.html#envs

### Organize codon table

Codon table retrieved from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N

Body of the above table's data saved as text file `codonTable.txt` 

The contents of `codonTable.txt` is parsed and processed by [`formatCodons.py`](formatCodons.py) into `codonTableCleaned.txt`:

```
( if [ ! -f "codons.tab" ]; then python3 formatCodons.py > "codonTableCleaned.txt"; fi )
```

## RUNNING CODE
