# MIMEqualityCheck
Plots for the sanity check of MIME data, including coverage, mutation rate, etc.

Running MIMEqualityCheck
-------------------

Run the R script 
```
Rscript dataSanityCheck.R <count_directory> <referenceFile> <result_directory> [<sample_sheet_file>]
```

withe following parameters

| parameter       | type          | description  |
| :---  |:---:| :----------------|
| count_directory         | (string)      |   directory where the subdirectories of the MIMEAnTo/MIMEAn2 results are lying |
| referenceFile         | (string)      |   reference sequence in fasta format |
| result_directory         | (string)      |   directory to save the generated plots |
| sample_sheet_file          | (string)      |   (optional) table with information about each sample |

The count_directory has to contain the countfiles in the subdirectories 1d and 2d. The countfiles contain the nucleotide (co-)occurrences of mapped reads, which can be inferred with the tool [sam2counts](https://github.com/maureensmith/sam2counts).
The 1d and 2d count file for the respective sample is named with the respective id/barcode:

```
/path/to/counts
+-- 1d
|   +-- 1.txt
|   +-- 2.txt
|   +-- 3.txt
|   +-- 4.txt
+-- 2d
|   +-- 1.txt
|   +-- 2.txt
|   +-- 3.txt
|   +-- 4.txt
```

The optional parameter sample_sheet_file is a semilcolon (";") separated file, with additional information about each sample. It has to contain the id or barcode in a column "Encdoding" and the actual sample name in the column "Sample". 
This is to show the names of the samples in the plots instead of the id.

In the result_directory several plots are saved, such as 

* the coverage per position per sample
* boxplots of the mutations frequency per sample
* the mutation frequency per mutations type per position per sample
* the Shannon entropy per position per sample
* ...