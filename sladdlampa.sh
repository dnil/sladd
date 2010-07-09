#!/bin/bash

export PERL5LIB=/home/ubelix/izb/nilsson/lib/lib64/perl5/site_perl/5.8.8:/home/ubelix/izb/nilsson/src/sladd

export BINDIR=~/src/sladd
export MAQBIN=~/bin/maq
export single=yes
export uorfutrcutofflen=2000
export runuorf=yes

../../src/sladd/run_maq.sh GDG-4.size24andup.fastq Tbrucei_TriTrypDB-1.1.fasta mRNA
