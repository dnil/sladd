#!/bin/bash

lib=1 ; 
for count in 1000 5000 25000 50000 100000 150000 200000 250000 300000 350000 400000 450000 500000 550000 600000 650000 700000 750000;
do
    echo GDG$lib - $count;
    ~/install/bin/randomly_pick_fastq.pl -n $count -N `grep -c ^@ ~/fasteris/GDG-${lib}/Size24andup.fastq` -o GDG-${lib}.Size24andup.rnd${count}.fastq /Users/daniel/fasteris/GDG-$lib/Size24andup.fastq
    maq fastq2bfq GDG-${lib}.Size24andup.rnd${count}.fastq GDG-${lib}.Size24andup.rnd${count}.bfq; maq map -n 2 GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map ~/tb/TbruceiGenomic_TriTrypDB-1.0.bfa GDG-${lib}.Size24andup.rnd${count}.bfq > GDG-${lib}.rnd${count}.maqmap.log; 
    maq mapview GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map > GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.txt; 
    ~/install/bin/maqmaptxt2gff.pl -s ~/tb/TbruceiGenomic_TriTrypDB-1.0.sizes.sillynames -l GDG${lib} -f maqmatchGDG${lib} < GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.txt > GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff ; 

    echo > tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab; 
    echo > GDG-${lib}.rnd${count}.maq.map.crunch.gff ;

    for seq in `cat ~/tb/fastalist` ; 
    do
	grep $seq GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff > ${seq}.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff ; 
	~/install/bin/per_gene_splicesite_info.pl -g ~/tb/$seq.TriTrypDB-1.0.gff -s ~/tb/$seq.fasta -t $seq.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.gene.gff;
	~/install/bin/tagged_gff_to_counts_tab.pl -g $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.gene.gff -s ~/tb/$seq.fasta -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab ; 
	echo $seq >> tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab; sort -k2n,2n $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab ;  
	~/install/bin/crunch_gff_with_counts.pl -g $seq.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff -s ~/tb/$seq.fasta -f maqmatchGDG${lib} -o $seq.GDG-${lib}.rnd${count}.maq.map.crunch.gff ;
	grep maqmatch $seq.GDG-${lib}.rnd${count}.maq.map.crunch.gff >> GDG-${lib}.rnd${count}.maq.map.crunch.gff; 

    done

    echo -n "GDG-${lib} ${count} " >> tb_ttdb_1.rnd.crunched_splicesites; 
    grep -c maqmatch GDG-${lib}.rnd${count}.maq.map.crunch.gff >> tb_ttdb_1.rnd.crunched_splicesites; 

    echo -n "GDG-${lib} ${count}" >> tb_ttdb_1.rnd.genes_with_counts; awk 'BEGIN {tick = 0} ($4 > 0) {tick=tick+1} END {print " ",tick}' tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.rnd.genes_with_counts;
    echo -n "GDG-${lib} ${count}" >> tb_ttdb_1.rnd.genes_with_counts_and_shortest_lt2k; awk 'BEGIN {tick = 0} ($4>0 && $11 < 2000) {tick=tick+1} END {print " ",tick}' tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.rnd.genes_with_counts_and_shortest_lt2k ;
done


lib=2 ; 
for count in 1000 5000 25000 50000 100000 150000 200000 250000 300000 350000 400000 450000 500000 550000 600000;
do
    echo GDG$lib - $count;
    ~/install/bin/randomly_pick_fastq.pl -n $count -N `grep -c ^@ ~/fasteris/GDG-${lib}/Size24andup.fastq` -o GDG-${lib}.Size24andup.rnd${count}.fastq /Users/daniel/fasteris/GDG-$lib/Size24andup.fastq
    maq fastq2bfq GDG-${lib}.Size24andup.rnd${count}.fastq GDG-${lib}.Size24andup.rnd${count}.bfq; maq map -n 2 GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map ~/tb/TbruceiGenomic_TriTrypDB-1.0.bfa GDG-${lib}.Size24andup.rnd${count}.bfq > GDG-${lib}.rnd${count}.maqmap.log; 
    maq mapview GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map > GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.txt; 
    ~/install/bin/maqmaptxt2gff.pl -s ~/tb/TbruceiGenomic_TriTrypDB-1.0.sizes.sillynames -l GDG${lib} -f maqmatchGDG${lib} < GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.txt > GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff ; 

    echo > tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab; 
    echo > GDG-${lib}.rnd${count}.maq.map.crunch.gff ;

    for seq in `cat ~/tb/fastalist` ; 
    do
	grep $seq GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff > ${seq}.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff ; 
	~/install/bin/per_gene_splicesite_info.pl -g ~/tb/$seq.TriTrypDB-1.0.gff -s ~/tb/$seq.fasta -t $seq.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.gene.gff;
	~/install/bin/tagged_gff_to_counts_tab.pl -g $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.gene.gff -s ~/tb/$seq.fasta -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab ; 
	echo $seq >> tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab; sort -k2n,2n $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab ;  
	~/install/bin/crunch_gff_with_counts.pl -g $seq.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff -s ~/tb/$seq.fasta -f maqmatchGDG${lib} -o $seq.GDG-${lib}.rnd${count}.maq.map.crunch.gff ;
	grep maqmatch $seq.GDG-${lib}.rnd${count}.maq.map.crunch.gff >> GDG-${lib}.rnd${count}.maq.map.crunch.gff; 

    done

    echo -n "GDG-${lib} ${count} " >> tb_ttdb_1.rnd.crunched_splicesites; 
    grep -c maqmatch GDG-${lib}.rnd${count}.maq.map.crunch.gff >> tb_ttdb_1.rnd.crunched_splicesites; 

    echo -n "GDG-${lib} ${count}" >> tb_ttdb_1.rnd.genes_with_counts; awk 'BEGIN {tick = 0} ($4 > 0) {tick=tick+1} END {print " ",tick}' tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.rnd.genes_with_counts;
    echo -n "GDG-${lib} ${count}" >> tb_ttdb_1.rnd.genes_with_counts_and_shortest_lt2k; awk 'BEGIN {tick = 0} ($4>0 && $11 < 2000) {tick=tick+1} END {print " ",tick}' tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.rnd.genes_with_counts_and_shortest_lt2k ;
done

lib=3 ;
for count in 1000 5000 25000 50000 100000 150000 200000 250000 300000 350000 400000 450000 500000 550000 600000 650000 700000 750000 800000 850000 900000 950000 1000000 1050000 1100000 1150000 ;
do
    echo GDG$lib - $count;
    ~/install/bin/randomly_pick_fastq.pl -n $count -N `grep -c ^@ ~/fasteris/GDG-${lib}/Size24andup.fastq` -o GDG-${lib}.Size24andup.rnd${count}.fastq /Users/daniel/fasteris/GDG-$lib/Size24andup.fastq
    maq fastq2bfq GDG-${lib}.Size24andup.rnd${count}.fastq GDG-${lib}.Size24andup.rnd${count}.bfq; maq map -n 2 GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map ~/tb/TbruceiGenomic_TriTrypDB-1.0.bfa GDG-${lib}.Size24andup.rnd${count}.bfq > GDG-${lib}.rnd${count}.maqmap.log; 
    maq mapview GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map > GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.txt; 
    ~/install/bin/maqmaptxt2gff.pl -s ~/tb/TbruceiGenomic_TriTrypDB-1.0.sizes.sillynames -l GDG${lib} -f maqmatchGDG${lib} < GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.txt > GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff ; 

    echo > tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab; 
    echo > GDG-${lib}.rnd${count}.maq.map.crunch.gff ;

    for seq in `cat ~/tb/fastalist` ; 
    do
	grep $seq GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff > ${seq}.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff ; 
	~/install/bin/per_gene_splicesite_info.pl -g ~/tb/$seq.TriTrypDB-1.0.gff -s ~/tb/$seq.fasta -t $seq.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.gene.gff;
	~/install/bin/tagged_gff_to_counts_tab.pl -g $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.gene.gff -s ~/tb/$seq.fasta -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab ; 
	echo $seq >> tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab; sort -k2n,2n $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab ;  
	~/install/bin/crunch_gff_with_counts.pl -g $seq.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff -s ~/tb/$seq.fasta -f maqmatchGDG${lib} -o $seq.GDG-${lib}.rnd${count}.maq.map.crunch.gff ;
	grep maqmatch $seq.GDG-${lib}.rnd${count}.maq.map.crunch.gff >> GDG-${lib}.rnd${count}.maq.map.crunch.gff; 

    done

    echo -n "GDG-${lib} ${count} " >> tb_ttdb_1.rnd.crunched_splicesites; 
    grep -c maqmatch GDG-${lib}.rnd${count}.maq.map.crunch.gff >> tb_ttdb_1.rnd.crunched_splicesites; 

    echo -n "GDG-${lib} ${count}" >> tb_ttdb_1.rnd.genes_with_counts; awk 'BEGIN {tick = 0} ($4 > 0) {tick=tick+1} END {print " ",tick}' tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.rnd.genes_with_counts;
    echo -n "GDG-${lib} ${count}" >> tb_ttdb_1.rnd.genes_with_counts_and_shortest_lt2k; awk 'BEGIN {tick = 0} ($4>0 && $11 < 2000) {tick=tick+1} END {print " ",tick}' tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.rnd.genes_with_counts_and_shortest_lt2k ;
done

lib=4 ;
for count in 1000 5000 25000 50000 100000 150000 200000 250000 300000 350000 400000 450000 500000 550000 600000 650000 700000 750000 800000 850000 900000;
do
    echo GDG$lib - $count;
    ~/install/bin/randomly_pick_fastq.pl -n $count -N `grep -c ^@ ~/fasteris/GDG-${lib}/Size24andup.fastq` -o GDG-${lib}.Size24andup.rnd${count}.fastq /Users/daniel/fasteris/GDG-$lib/Size24andup.fastq
    maq fastq2bfq GDG-${lib}.Size24andup.rnd${count}.fastq GDG-${lib}.Size24andup.rnd${count}.bfq; maq map -n 2 GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map ~/tb/TbruceiGenomic_TriTrypDB-1.0.bfa GDG-${lib}.Size24andup.rnd${count}.bfq > GDG-${lib}.rnd${count}.maqmap.log; 
    maq mapview GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map > GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.txt; 
    ~/install/bin/maqmaptxt2gff.pl -s ~/tb/TbruceiGenomic_TriTrypDB-1.0.sizes.sillynames -l GDG${lib} -f maqmatchGDG${lib} < GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.txt > GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff ; 

    echo > tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab; 
    echo > GDG-${lib}.rnd${count}.maq.map.crunch.gff ;

    for seq in `cat ~/tb/fastalist` ; 
    do
	grep $seq GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff > ${seq}.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff ; 
	~/install/bin/per_gene_splicesite_info.pl -g ~/tb/$seq.TriTrypDB-1.0.gff -s ~/tb/$seq.fasta -t $seq.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.gene.gff;
	~/install/bin/tagged_gff_to_counts_tab.pl -g $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.gene.gff -s ~/tb/$seq.fasta -f maqmatchGDG${lib} -o $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab ; 
	echo $seq >> tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab; sort -k2n,2n $seq.tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab ;  
	~/install/bin/crunch_gff_with_counts.pl -g $seq.GDG-${lib}.rnd${count}.tb_ttdb_1.maq.map.gff -s ~/tb/$seq.fasta -f maqmatchGDG${lib} -o $seq.GDG-${lib}.rnd${count}.maq.map.crunch.gff ;
	grep maqmatch $seq.GDG-${lib}.rnd${count}.maq.map.crunch.gff >> GDG-${lib}.rnd${count}.maq.map.crunch.gff; 

    done

    echo -n "GDG-${lib} ${count} " >> tb_ttdb_1.rnd.crunched_splicesites; 
    grep -c maqmatch GDG-${lib}.rnd${count}.maq.map.crunch.gff >> tb_ttdb_1.rnd.crunched_splicesites; 

    echo -n "GDG-${lib} ${count}" >> tb_ttdb_1.rnd.genes_with_counts; awk 'BEGIN {tick = 0} ($4 > 0) {tick=tick+1} END {print " ",tick}' tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.rnd.genes_with_counts;
    echo -n "GDG-${lib} ${count}" >> tb_ttdb_1.rnd.genes_with_counts_and_shortest_lt2k; awk 'BEGIN {tick = 0} ($4>0 && $11 < 2000) {tick=tick+1} END {print " ",tick}' tb_ttdb_1.GDG-${lib}.rnd${count}.maq.map.counts.tab >> tb_ttdb_1.rnd.genes_with_counts_and_shortest_lt2k ;
done
