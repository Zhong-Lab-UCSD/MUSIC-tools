# MUSIC-tools (data analysis for MUSIC data)
<img src="https://github.com/Irenexzwen/MUSIC-tools/assets/24513822/fad5c93a-945f-461d-8021-d31bef367113" width="500">

# 1.Introduction
This repository includes major code we used to analyze MUSIC data associated with the manuscript: "Single-cell multiplex chromatin and RNA interactions in aging human brain" (Wen et al,. 2024)[].

# 2.MUSIC data library configuration
http://sysbiocomp.ucsd.edu/public/wenxingzhao/MUSIC_docker/intro.html#library-configuration.

# 3.Raw data processing
MUSIC derives a paired-end sequencing library which includes cell barcodes, complexes barcodes and DNA/RNA insert sequences. The first step to analyze MUSIC data requires extract all essential information from the paired-end fastq files. Step by step data processing documentation is here: http://sysbiocomp.ucsd.edu/public/wenxingzhao/MUSIC_docker/step.html. 

## 3.1 MUSIC-docker
We developed a docker container called MUSIC-docker which is a wrapper of a Snakemake pipeline that enables raw data processing. A full documentation of MUSIC-docker can be found here: [http://sysbiocomp.ucsd.edu/public/wenxingzhao/MUSIC_docker/index.html](http://sysbiocomp.ucsd.edu/public/wenxingzhao/MUSIC_docker/rundocker.html). You can also make your own updates and re-compile to build your own version of MUSIC-docker using source code here: https://github.com/Irenexzwen/MUSIC-docker

## 3.2 MUSIC Snakemake pipeline
If you are familiar with Snakemake pipeline, we also provided a Snakemake pipeline where you can update the config files according to your local working enviroment. 


