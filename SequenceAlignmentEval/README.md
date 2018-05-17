# A evaluation tool for evaluating sequence alignment outcome
### Introduction

This project implement a tool to evaluate DNA sequence alignment result given a reference genbank file.

### Build Everything

Build everything with the following command
```Shell
make -j2
```

### Evaluating Sequence Alignment Result

Evaluate DNA sequence alignment result with the following command
```Shell
./evaluation -r genbank/NC_000909.gbk -d <path/to/DNA/sequence/file/in/FASTA/format> -o <path/to/evaluation/report>
```

The reference file is a genbank file specifying all open reading frames on DNA. The DNA sequence alignment result given to this program must be in FASTA format.

