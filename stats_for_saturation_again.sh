#!/bin/bash

# for ubelix
export PERL5LIB=/home/ubelix/izb/nilsson/lib/lib64/perl5/site_perl/5.8.8:/home/ubelix/izb/nilsson/src/sladd
export BINDIR=~/src/sladd
export MAQBIN=~/bin/maq

# general
export single=yes
export uorfutrcutofflen=2000
export runuorf=no

export forceupdate=no
 
ref=Tbrucei_TriTrypDB-1.1.fasta
lib=GDG-9.Size23andup.fastq

function needsUpdate()
{
    # USAGE: needsUpdate(target, prereq [, prereq]*)
    # return true (needsupdate=yes) if target does not yet exist, is older than its prereqs or forceupdate=yes is in effect.

    needsupdate="no"
    
    if [ "$forceupdate" = "yes" ] 
    then
	needsupdate="yes"
    fi

    target=$1;
    
    for prereq in ${@:2}
    do
	if [ $target -ot $prereq ]
	then
	    needsupdate="yes"
	fi
   done
    
    [ "$needsupdate" = "yes" ]
}

orig_lib_size=`grep -c ^@ $lib`

# for count in 1000 5000 25000 50000 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000 1100000 1200000 1300000 1400000 1500000 1600000 1800000 2000000 2200000;

declare -a rndlibs

# pick a COUNT random reads from lib fastq, perform ordinary mapping using these.
for count in 1000 5000 10000 25000 50000 100000 200000 300000 400000 500000 600000 800000 1000000 1200000 1400000 1600000 1800000;
do
    echo ${lib%%.fastq} - $count;
    
    # TODO: change to lib%.fastq, but only if there is to be a new major rerun..
    rndlib=${lib}.rnd${count}.fastq

    rndlibs[${#rndlibs[@]}]=$rndlib
    if needsUpdate $rndlib $lib $BINDIR/randomly_pick_fastq.pl
    then
	$BINDIR/randomly_pick_fastq.pl -n $count -N $orig_lib_size -o $rndlib $lib
    fi

    runmaq=no

    crunchfile=${rndlib%.fastq}_vs_${ref%.fasta}.maqmap.crunch.gff
    if needsUpdate $crunchfile $rndlib $ref
    then
	runmaq=yes
    fi

    genetab=${rndlib%.fastq}_vs_${ref%%.fasta}.maqmap.counts.genes_with_major_UTR_lt2k.tab
    if needsUpdate $genetab $rndlib $ref
    then
	runmaq=yes
    fi

    if [ "$runmaq" = "yes" ]
    then
	$BINDIR/run_maq.sh ${rndlib} ${ref} mRNA
    fi

done


# Update summary
rndsummary=${ref}_vs_${lib}.rnd.summary

if needsUpdate $rndsummary $lib $ref $rndlibs[@]
then
    echo "rnd_count=Splicesites=Genes" |sed -e 's/=/\t/g' > $rndsummary
    for rndlib in ${rndlibs[@]}
    do 
	counttmp=${rndlib##*.rnd}
	count=${counttmp%%.fastq}

	crunchfile=${rndlib%.fastq}_vs_${ref%.fasta}.maqmap.crunch.gff
	splicesites=`grep -c maqmatch ${crunchfile}`

	genetab=${rndlib%.fastq}_vs_${ref%%.fasta}.maqmap.counts.genes_with_major_UTR_lt2k.tab
	genes=`grep -v ^Name ${genetab}|cut -f1 |sort|uniq |wc -l`

	echo $count"="$splicesites"="$genes
    done |sed -e 's/=/\t/g' >> $rndsummary
fi
