# vfkmrz
fast tool for kmerizing genomes or metagenomes
* optimized for speed
* tunable RAM use 
* support two popular formats: fasta and fastq

# vfkmer_fasta
the program vfkmrz_fasta.cpp scans its input (fasta format) for forward k mers,

usage:
```shell
    $ g++ -O3 --std=c++11 -o vfkmrz_fasta vfkmrz_fasta.cpp
    $ cat /path/exp.fasta | ./vfkmrz_fasta
```

# vfkmer_fastq
the program vfkmrz_fastq.cpp scans its input (fastq text stream) for forward k mers,

usage:

```shell
    $ g++ -O3 --std=c++11 -o vfkmrz_fastq vfkmrz_fastq.cpp
    $ gzip -dc /path/exp.fastq.gz | ./vfkmrz_fastq
```

-or-

```shell
    $ zcat /path/exp.fastq.gz | ./vfkmrz_fastq
```

-or- 

```shell
    $ cat /path/exp.fastq | ./vfkmrz_fastq
```

# notes
standard fasta and fastq format only for inputs, otherwise failure is almost guaranteed. 

