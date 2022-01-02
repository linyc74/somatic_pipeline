# GATK Pipeline

**Custom-built GATK pipeline**

```commandline
python gatk_pipeline \
  -r path/to/reference_genome.fa \
  -1 path/to/tumor.1.fq.gz \
  -2 path/to/tumor.2.fq.gz \
  -3 path/to/normal.1.fq.gz \
  -4 path/to/normal.1.fq.gz \
  -o path/to/outdir \
  -t cpu_threads

required arguments:
  -r REF_FA, --ref-fa REF_FA
                        path to the reference genome fasta file
  -1 TUMOR_FQ1, --tumor-fq1 TUMOR_FQ1
                        path to the tumor read 1 fastq file
  -2 TUMOR_FQ2, --tumor-fq2 TUMOR_FQ2
                        path to the tumor read 2 fastq file
  -3 NORMAL_FQ1, --normal-fq1 NORMAL_FQ1
                        path to the normal read 1 fastq file
  -4 NORMAL_FQ2, --normal-fq2 NORMAL_FQ2
                        path to the normal read 2 fastq file

optional arguments:
  -o OUTDIR, --outdir OUTDIR
                        path to the output directory (default: gatk_pipeline_outdir)
  -t THREADS, --threads THREADS
                        number of CPU threads (default: 4)
  -d, --debug           debug mode
  -h, --help            show this help message
  -v, --version         show version
```
