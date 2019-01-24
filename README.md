# PhageContentCalculator
A Python script designed to calculate the prophage content within bacterial genomes. 

PhageContentCalculator script uses a database of previously known full-length and satellite prophage genomes in order to estimate the total percentage of phage-related ORFs in the bacterial genome of interest.


**Author**: Reza Rezaei Javan

**E-mail**: reza.rezaeijavan@ndm.ox.ac.uk

Copyright (C) 2019 University of Oxford

## Usage
Use the following command to start PhageMiner:
```
python PhageContentCalculator.py <input file>
```
## Input Files
### The bacterial genome of interest
The bacterial genome needs to be provided in a fasta format. A *Streptococcus pneumoniae* draft genome (multiple contigs) is provided as an example input file.

### Database of prophage genomes 
Prophage genomes for the database need to be annotated and provided in .GFF format. They must contain the nucleotide sequence at the end of the file. Files annotated using the [PROKKA](https://github.com/tseemann/prokka) (Rapid Prokaryotic Genome Annotation) pipeline are valid with PhageContentCalculator and this is the recommended way of generating these files. 

## Installation

### Required dependencies
PhageContentCalculator has the following dependencies, which need to be installed in advance:
* [PROKKA](https://github.com/tseemann/prokka)
* [Roary](https://sanger-pathogens.github.io/Roary/)

Before running the script, make sure that the script knows where these programmes are installed on your computer. The correct directories can be set at the top of the script.  

## Output files
PhageContentCalculator output the total percentage of full-length and satellite prophage related ORFs in the given bacterial genome. Furthermore, a tab-separated values (TSV) file named "results.txt" is produced, in which the following details are recorded:

| 1st column | 2nd column | 3rd column | 4th column | 5th column | 6th column | 7th column
| --- | --- | --- | --- | --- | --- | --- |
| Genome name | Percentage of full-length prophage related genes | Percentage of satellite prophage related genes | Total number of genes in the bacterial genome | Total number of genes in the bacterial genome | Number of full-length prophage related genes| Number of satellite prophage related genes| Unique bactrial genes|

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/RezaRezaeiJavan/PhageMiner/issues) or email reza.rezaeijavan@ndm.ox.ac.uk.
