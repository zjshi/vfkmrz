# vfkmrz
fast tool for kmerizing genomes or metagenomes in FASTA or FASTQ format

the program vfkmrz_fasta.cpp scans its input (fasta format) for forward k mers,

usage:
   g++ -O3 --std=c++11 -o vfkmrz_fasta vfkmrz_fasta.cpp
   cat /path/exp.fasta | ./vfkmrz_fasta

standard fasta format only for input, otherwise failure is almost guaranteed. 

the program vfkmrz_fastq.cpp scans its input (fastq text stream) for forward k mers,

usage:
   g++ -O3 --std=c++11 -o vfkmrz_fastq vfkmrz_fastq.cpp
   gzip -dc /path/exp.fastq.gz | ./vfkmrz_fastq
   or 
   zcat /path/exp.fastq.gz | ./vfkmrz_fastq
   or 
   cat /path/exp.fastq | ./vfkmrz_fastq
   the output path can be specified in a variable below

standard fastq format only for input, otherwise failure is almost guaranteed. 

