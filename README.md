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

# fastq2fasta
the program fastq2fasta.cpp scans fastq format inputs and convert them into fasta format outputs.

usage:

```shell
    $ g++ -O3 --std=c++11 -o fastq2fasta fastq2fasta.cpp
    $ gzip -dc /path/exp.fastq.gz | ./fastq2fasta
```

-or-

```shell
    $ zcat /path/exp.fastq.gz | ./fastq2fasta
```

-or- 

```shell
    $ cat /path/exp.fastq | ./fastq2fasta
```

# htfs_stats 
the program htsf_stats scans inputs (fasta or fastq format) and output basic stats information about the inputs. (e.g. number of lines, number of reads and number of bases)


```shell
    $ g++ -O3 --std=c++11 -o htsf_stats htsf_stats.cpp
    $ gzip -dc /path/exp.fastq.gz | ./htsf_stats -fq
```

-or-

```shell
    $ cat /path/exp.fasta | ./htsf_stats -fa
```

# vfkmrz_match
the prgram vfkmrz_match takes two list of kmers, db and query list, that have the same length, then search query list against db list and return the occurrences of db kmers in query kmers.

```shell
    $ g++ -O --std=c++14 -o vfkmrz_match vfkmrz_match.cpp 
    $ cat </path/to/kmer/list> | ./vfkmrz_match
```

# vfkmrz_bunion
the prgram vfkmrz_bunion takes two list of kmers, specified by option of -k1 and -k2, that have the same length, then return the unique kmers for the union of k1 and k2 kmers 

```shell
    $ g++ -O --std=c++11 -o vfkmrz_bunion vfkmrz_bunion.cpp
    $ ./vfkmrz_bunion -k1 </path/to/kmer/list1> -k2 </path/to/kmer/list2>
```


# notes
standard fasta and fastq format only for inputs, otherwise failure is almost guaranteed. 
the programs require the compilers that are compatible with C++ 14 standards.
all the tests have been done on clang-900.0.38

vfkmrz programs are now compatible for GNU C Compiler. All the programs were tested and compiled successfully under GCC 7.2.0; All the programs except for vfkmrz_match.cpp (requires C++ standard newer than 14) were compiled successfully under GCC 4.8.5  
