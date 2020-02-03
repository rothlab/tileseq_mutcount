## TileSeq mutation count package

This package is made to parse input sequecning files (fastq) with user provided parameter.csv file.
Output of this pipeline is mutation counts for each pair of fastq files. 

## Dependencies
---

`python 3.6`
`R`

### Execution
---

To run this pipeline on BC2

`python main.py -p param_json -o output_folder -f fastq_file_folder`

```
Example: 

python main.py -p $HOME/dev/tilseq_mutcount/190506_param_MTHFR.json -o $HOME/dev/tilseq_mutcount/output/ -f $HOME/tileseq_data/WT/

Required parameters: #required from users

-p PARAM 
-o OUTPUT Output directory
-f FASTQ Path to fastq files you want to analyze

Optional parameters:
-h, --help Show help message

-n N_READS Number of reads to down-sample to (default = 30000)
-env ENVIRONMENT Name of cluster you want to run this script (default = BC2)
-log LOG_LEVEL Log level of the logging object: debug, info, warning, error, critical (default = info)

--skip_alignment Skip alignment for this run. Alignment output already exsist and the output path should be the output generated by a previous run

Example of skipping alignment: 

python main.py -p $HOME/dev/tilseq_mutcount/190506_param_MTHFR.json -o /home/rothlab1/rli/dev/tilseq_mutcount/output/190506_MTHFR_WT_2020-01-29-17-07-04/ -f $HOME/tileseq_data/WT/ --skip_alignment
```


### Input files
---

`/path/to/fastq/` - Full path to input fastq files 

`param.csv` - CSV file contains information for this run (please see example
[here](https://docs.google.com/spreadsheets/d/1tIblmIFgOApPNzWN2KUwj8BKzBiJ1pOL7R4AOUGrqvE/edit?usp=sharing)
).
This file is required to be comma-seperated and saved in csv format. 


### Output files
---

One ouptut folder is created for each run. The output folder are named with `project_name-time-stamp`

Within each output folder, the following files and folders:

`main.log` - main logging output 

`ref` - Reference fasta file and bowtie2 index

`ds_fastq` - Raw fastq files are downsampled to `n` reads and the downsampled fastq files are saved
in this folder

`sam_files` - Alignemnt output and log files for the raw fastq files

`ds_sam_files` - Alignment output and log files for the down-sampled fastq files


### Alignment
---

The pipeline takes the sequence in the parameter file as reference and try to align the fastq files
to the whole reference sequence. 
