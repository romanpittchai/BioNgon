# BioNgon 

The project BioNgon in the future, an expandable project, the task of which is to create simple, free tools for tasks related to biology and activities related to the field of biology. The project is developing and I hope it will benefit all those who are interested in this benefit.


## Information about the project

## _GC Content Calculator_

For testing, data obtained in the public domain was used (source: Institute of Genetics, University of California, Santa Cruz).

- Linux

#### GUI

A small program implemented in the C programming language, GTK+ graphical shell, cairo library and sqlite3 source code. In addition to the source code of the program itself, the repository also has a simple makefile for the convenience of building the program, as well as a service file (`gc_content_calculator.cbp`) for convenient operation in the CodeBlocks IDE. There is a file with the interface in a separate folder (created in Glade).

The result of the calculation:

```
Filename: Homo_sapiens.GRCh38.dna.chromosome.1.fa

File header: 1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF

Date of processing:
10:13 - 11.01.2024

A - 933 - 20,30%
U - 0 - 0,00%
G - 1223 - 26,62%
C - 1452 - 31,60%
T - 953 - 20,74%
Total number of characters - 4595
Number of ATGCU - 4561
Other characters - 34
Number of GC content in characters - 2675
Number of GC content in percent - 58,22%
```
Diagram:
![Example of drawing a diagram](examples/dna)



#### CLI

A small program written in Rust. At the moment, it just counts and outputs data to the console. In the future, it is planned to work with the database, text data, and graph output.

The result of the calculation:

```
*************************
Filename: Homo_sapiens.GRCh38.dna.chromosome.1.fa
-------------------------
File header: 1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
-------------------------
Date of processing: 21:37 - 12.02.2024
-------------------------
Total number of characters 248956422 - 100%
Number of ATGCU 230481012 - 92.58%
Number of GC content in characters 96166571 - 38.63%
A - 67070277 - 26.94%
C - 48055043 - 19.3%
G - 48111528 - 19.33%
T - 67244164 - 27.01%
U - 0 - 0%
Other characters - 18475410 - 7.42%
*************************
```




## Contact

Open to cooperation and assistance in development. 
Email: dnlitdil@gmail.com


> ***_Note:_*** _The project uses the MIT License._
