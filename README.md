# iDARTS
### individualized Deep-learning Analysis of RNA Transcript Splicing
### Date: "07.01.2021"

### Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contact](#contact)
- [Copyright and License Information](#copyright-and-license-information)

### Installation (Under development)
The installation is made easy through [Anaconda](https://anaconda.org/idarts).

```bash
conda create -n idarts python=2.7  # optional
conda activate idarts              # optional
conda install -c idarts
```

### Usage
Detailed arguments:
```
usage: iDARTS [-h] [--version]
              {get_resources,build_feature,predict,parse_vcf}
positional arguments:

get_resources
    usage:iDARTS get_resources [-h] [-o OUT_DIR]

    optional arguments:
      -h, --help            show this help message and exit
      -o OUT_DIR, --out-dir OUT_DIR
                            Optional, default user home directory: Output folder
                            for downloaded data

parse_vcf
  usage: iDARTS parse_vcf [-h] [-t {SE,A5SS,A3SS,RI}] -i INPUT -v VCF_PATH -o
                          OUTPUT  

  optional arguments:
    -h, --help            show this help message and exit
    -t {SE,A5SS,A3SS,RI}, --type {SE,A5SS,A3SS,RI}
                          Optional, default SE: specify the alternative splicing
                          event type. SE: skipped exons, A3SS: alternative 3
                          splice sites, A5SS: alternative 5 splice sites, RI:
                          retained introns. A5SS, A3SS, and RI are under development.
    -i INPUT, --input INPUT
                          A list of alternative splicing events; iDARTS parse
                          SNVs from vcf for alternative splicing events
    -v VCF_PATH, --vcf_path VCF_PATH
                          vcf path
    -o OUTPUT, --out-file-name OUTPUT
                          parsed vcf output file name

build_feature
  usage: iDARTS build_feature [-h] [-t {SE,A5SS,A3SS,RI}] -i INPUT
                              [-m {True,False}] -o OUTPUT 

  optional arguments:
    -h, --help            show this help message and exit
    -t {SE,A5SS,A3SS,RI}, --type {SE,A5SS,A3SS,RI}
                          Optional, default SE: specify the alternative splicing
                          event type. SE: skipped exons, A3SS: alternative 3
                          splice sites, A5SS: alternative 5 splice sites, RI:
                          retained introns
    -i INPUT, --input INPUT
                          A list of alternative splicing events; iDARTS build
                          feature
    -m {True,False}, --mutate {True,False}
                          Annotate the sequence features with SNV (300nt within
                          exon-intron boundary or on exons)
    -o OUTPUT, --out-file-name OUTPUT
                          feature annotation output file name   

predict
  usage: iDARTS predict [-h] [-t {SE,A5SS,A3SS,RI}] -i INPUT [-e EXPR] -o OUTPUT  

  optional arguments:
    -h, --help            show this help message and exit
    -t {SE,A5SS,A3SS,RI}, --type {SE,A5SS,A3SS,RI}
                          Optional, default SE: specify the alternative splicing
                          event type. SE: skipped exons, A3SS: alternative 3
                          splice sites, A5SS: alternative 5 splice sites, RI:
                          retained introns
    -i INPUT, --input INPUT
                          A list of annotated alternative splicing features
    -e EXPR, --expression EXPR
                          Expressing file (TPM value from Kallisto);header
                          format 'Gene_ID\tExp1,Exp2,Exp3...(different
                          expression profiles separated by comma)'
    -o OUTPUT, --out-file-name OUTPUT
                          iDARTS prediction output file name                          
```

### Contact
[Mailing List](mailto:idarts-user-group@googlegroups.com) / [Group](https://groups.google.com/d/forum/idarts-user-group)

[GitHub Issues](https://github.com/Xinglab/iDARTS/issues)

Zhicheng Pan <zc.pan@ucla.edu>

Yi Xing <xingyi@chop.edu>

### Copyright and License Information
Copyright (C) 2021 The Children’s Hospital of Philadelphia

Authors: Zhicheng Pan and Yi Xing

This program is licensed with commercial restriction use license. Please see the attached LICENSE file for details.
