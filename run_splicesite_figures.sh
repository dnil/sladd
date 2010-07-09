#!/bin/bash

SEQLOGOBIN=~/install/weblog/seqlogo

# extract splicesites

for lib in 1 2 3 4 ; do for seq in `cat ../fastalist` ; do ~/install/bin/extract_splicesites.pl -g $seq.tb_ttdb_1.GDG${lib}.maq.map.gff -s ../$seq.fasta -f maqmatchGDG${lib} -1 -m -o $seq.tb_ttdb_1.GDG-${lib}_24up.maq.map.once.major.splicesites ; done ; done

for lib in 1 2 3 4 ; do echo > tb_ttdb_1.GDG-${lib}_24up.maq.map.once.major.splicesites ; for seq in `cat ../fastalist` ; do cat $seq.tb_ttdb_1.GDG-${lib}_24up.maq.map.once.major.splicesites >> tb_ttdb_1.GDG-${lib}_24up.maq.map.once.major.splicesites; done ; done

../../install/weblogo/seqlogo -f Tb927_01_v4.tb_ttdb_1.GDG-1_24up.maq.map.once.splicesites -h 8 -w 26 -n -Y -c -k1 -s -150 -d 1 -t "Consensus splice site" -F PDF -o Tb927_01_v4.tb_ttdb_1.GDG-1_24up.maq.map.once.splicesites

for lib in 1 2 3 4 ; do ../../install/weblogo/seqlogo -f tb_ttdb_1.GDG-${lib}_24up.maq.map.once.major.splicesites -h 8 -w 26 -n -Y -c -k1 -s -150 -d 1 -t "Consensus splice site once major GDG-${lib}" -F PDF -o tb_ttdb_1.GDG-${lib}_24up.maq.map.once.major.splicesites.out; done

for lib in 1 2 3 4 ; do for seq in `cat ../fastalist` ; do ~/install/bin/extract_splicesites.pl -g $seq.tb_ttdb_1.GDG${lib}.maq.map.gff -s ../$seq.fasta -f maqmatchGDG${lib} -m -o $seq.tb_ttdb_1.GDG-${lib}_24up.maq.map.major.splicesites ; done ; done

for lib in 1 2 3 4 ; do ../../install/weblogo/seqlogo -f tb_ttdb_1.GDG-${lib}_24up.maq.map.major.splicesites -h 8 -w 26 -n -Y -c -k1 -s -150 -d 1 -t "Consensus splice site major GDG-${lib}" -F PDF -o tb_ttdb_1.GDG-${lib}_24up.maq.map.major.splicesites.out; done

for lib in 1 2 3 4 ; do echo > tb_ttdb_1.GDG-${lib}_24up.maq.map.internal.splicesites ; for seq in `cat ../fastalist` ; do cat $seq.tb_ttdb_1.GDG-${lib}_24up.maq.map.major.splicesites >> tb_ttdb_1.GDG-${lib}_24up.maq.map.major.splicesites; done ; done

for lib in 1 2 3 4 ; do for seq in `cat ../fastalist` ; do ~/install/bin/extract_splicesites.pl -g $seq.tb_ttdb_1.GDG${lib}.maq.map.gff -s ../$seq.fasta -f maqmatchGDG${lib} -E -o $seq.tb_ttdb_1.GDG-${lib}_24up.maq.map.internal.splicesites ; done ; done

for lib in 1 2 3 4; do echo GDG-$lib; sort tb_ttdb_1.GDG-${lib}_24up.maq.map.internal.splicesites |uniq |wc -l; done

for lib in 1 2 3 4 ; do for seq in `cat ../fastalist` ; do ~/install/bin/extract_splicesites.pl -g $seq.tb_ttdb_1.GDG${lib}.maq.map.gff -s ../$seq.fasta -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}_24up.maq.map.splicesites ; done ; done

for lib in 1 2 3 4 ; do ../../install/weblogo/seqlogo -f tb_ttdb_1.GDG-${lib}_24up.maq.map.once.splicesites -h 8 -w 26 -n -Y -c -k1 -s -150 -d 1 -t "Consensus splice site once GDG-${lib}" -F PDF -o tb_ttdb_1.GDG-${lib}_24up.maq.map.once.splicesites.out; done
