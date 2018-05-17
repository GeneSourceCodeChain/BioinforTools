# BioinforTools
Several Tools set for Bioinformatics 

### Sequence Evaluation Tool

Tool for evaluating DNA alignment result. The evaluation result is output in a text report.

### Exon 

exon_gatk.sh

The script uses bwa, samtools, gatk and other tools to analyze exon mutation analysis.

### WGS

pipline.sh

The script is for full genome variation analysis. The program used is speedseq, which is fast and accurate.

### Python scripts

format.py

The scripting effect is to convert the dictionary format into rows of data.

### Shell scripts

unzip.sh
There are multiple compressed files in the directory, and the file names contained in each compressed file are the same and the contents are different.
At the same time decompress multiple compression files and rename the extracted files.

cat.sh
Merge multiple files into a file in order.

paste.sh
Merge multiple files with paste.
Append the next file to the back of the previous file rather than below.

grep.sh
Use grep to extract the same line from the beginning string in a file to the new file.

grep_diff.sh
Find a file that is different from other files, and write the results to the new file.
If you compare only one file with the other file, remove the comment from the line and use the script.
exec_grep_diff.sh
If you want to compare arbitrary files to other files, put two scripts in the same directory and execute the script.

compare.sh
The script is to compare the Y chromosome genotypes of multiple individuals, and to identify individuals with identical genotypes.

