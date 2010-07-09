#!/bin/bash
export PERL5LIB=/home/ubelix/izb/nilsson/lib/lib64/perl5/site_perl/5.8.8:/home/ubelix/izb/nilsson/src/sladd

export BINDIR=~/src/sladd
export MAQBIN=~/bin/maq
export single=yes
export uorfutrcutofflen=2000

lib=GDG-9.Size23andup

for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do ../../src/sladd/extract_splicesites.pl -g ${lib}*maqmap.${seq}.gene.gff -s Tbrucei_TriTrypDB-1.1.${seq}.fasta -f maqmatch${lib%%.*} -1 -m -c mRNA -o $lib.maqmap.$seq.major.once.splicesites.out; done

echo -n > $lib.maqmap.major.once.splicesites.out; for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do cat $lib.maqmap.$seq.major.once.splicesites.out >> $lib.maqmap.major.once.splicesites.out; done

for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do ../../src/sladd/extract_splicesites.pl -g ${lib}*maqmap.${seq}.gene.gff -s Tbrucei_TriTrypDB-1.1.${seq}.fasta -f maqmatch${lib%%.*} -1  -c mRNA -o $lib.maqmap.$seq.once.splicesites.out; done

echo -n > $lib.maqmap.once.splicesites.out; for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do cat $lib.maqmap.$seq.once.splicesites.out >> $lib.maqmap.once.splicesites.out; done

for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do ../../src/sladd/extract_splicesites.pl -g ${lib}*maqmap.${seq}.gene.gff -s Tbrucei_TriTrypDB-1.1.${seq}.fasta -f maqmatch${lib%%.*} -1 -E -c mRNA -o $lib.maqmap.$seq.internal.once.splicesites.out; done

echo -n > $lib.maqmap.internal.once.splicesites.out; for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do cat $lib.maqmap.$seq.internal.once.splicesites.out >> $lib.maqmap.internal.once.splicesites.out; done

for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do ../../src/sladd/extract_splicesites.pl -g ${lib}*maqmap.${seq}.gene.gff -s Tbrucei_TriTrypDB-1.1.${seq}.fasta -f maqmatch${lib%%.*} -1 -I -c mRNA -o $lib.maqmap.$seq.external.once.splicesites.out; done

echo -n > $lib.maqmap.external.once.splicesites.out; for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do cat $lib.maqmap.$seq.external.once.splicesites.out >> $lib.maqmap.external.once.splicesites.out; done

for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do ../../src/sladd/extract_splicesites.pl -g ${lib}*maqmap.${seq}.gene.gff -s Tbrucei_TriTrypDB-1.1.${seq}.fasta -f maqmatch${lib%%.*} -1 -M -c mRNA -o $lib.maqmap.$seq.minor.once.splicesites.out; done

echo -n > $lib.maqmap.minor.once.splicesites.out; for seq in `cat Tbrucei_TriTrypDB-1.1.fastalist`; do cat $lib.maqmap.$seq.minor.once.splicesites.out >> $lib.maqmap.minor.once.splicesites.out; done

for file in *splicesites.out ; do echo $file ; $BINDIR/splicesite_stats.pl < $file ; done > $lib.splicesite_stats.summary

# ../../install/weblogo/seqlogo -f $lib.maqmap.major.once.splicesites.out -h 8 -w 26 -n -Y -c -k1 -s -90 -d 1 -t "Consensus splice site major once ${lib%%.*}" -F PDF -o $lib.major.once.splicesites.out
