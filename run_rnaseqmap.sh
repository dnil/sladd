#!/bin/bash
#
# Daniel Nilsson, 2009
#
# daniel.nilsson@izb.unibe.ch, daniel.k.nilsson@gmail.com
#
# USAGE: run_rnaseqmap.sh tagfile.fastq referenecesequence.fasta [gene_feature_type(gene)]
#
# Expects: * referenceseqeunce.gff in same dir as referencesequence.fasta
#          * maq binary installed in $MAQBIN
#          * several of the other "sladdlampa"-scripts in $BINDIR
#          * indirectly, bioPerl in the perl library path (and a perl 5.8.x interpreter at /usr/bin)
#
# Environment: Several bash environment variables influence the pipeline behaviour
#            BINDIR           (~/install/bin)
#            MAQBIN           (~/install/maq-0.7.1/maq)
#            alqt             (30)
#            single           (yes)
#            forceupdate      (no)
#

# check input environment variables and set unset ones to default values

# BINDIR holds the scripts/binaries for this pipeline
if [ -z "$BINDIR" ]
then   	
	#BINDIR=~/src/sladd
	BINDIR=~/install/sladd
fi

export PERL5LIB=$PERL5LIB:$BINDIR

# MAQBIN is the MAQ installation directory (path to the maq binary)
if [ -z "$MAQBIN" ]
then
	#MAQBIN=~/bin/maq
	MAQBIN=~/install/maq-0.7.1/maq
fi

# SEQLOGO is the weblogo path to the Weblogo seqlogo app
if [ -z "$SEQLOGOBIN" ]
then
	SEQLOGOBIN=~/install/weblogo/seqlogo
fi

# MAQ alignment quality cutoff
# if 0, multimapping is no longer a special case..
if [ -z "$alqt" ]
then
    alqt=30
fi

# single mapping filtered files to be used for tab generation
if [ -z "$single" ]
then
    single=yes
fi

# set forceupdate=yes to run all available analyses, even if the file modification times advise against it 
if [ -z "$forceupdate" ]
then
    forceupdate=no
fi

# set desired number of concurrent processes
# prefer an already set value of $NPROC, or use nproc to return it if available
if [ -z "$NPROC" ]
then
    NPROCBIN=`which nproc`
    if [ -x $NPROCBIN ] 
    then
	NPROC=`$NPROCBIN`
    fi

    if [ -z "$NPROC" ] 
    then 
	NPROC=1
    fi
fi

# CALLED will contain a complete pathname for this script
CALLED=$0

# command line arguments

if [ $# -lt 2 ]
then
        echo "USAGE: ${0##*/} tagfile.fastq refsequence.fasta [gene_gff_type]"
        exit 1
fi

# first argument : tagfilename.fastq : name convention is libname.sizefilter.fastq
tagfile=$1 

# second argument : refsequencename.fasta, also expects to find refsequencename.gff in the same directory
reffile=$2 

if [ ! -e "$tagfile" ]
then
    echo "Tagfile $tagfile not found."
    exit 1
fi

if [ ! -e "$reffile" ]
then
    echo "Reference file $reffile not found."
    exit 1
fi

refgff=${reffile%fasta}gff
if [ ! -e "$refgff" ]
then
    echo "Reference GFF $refgff not found."
    exit 1
fi

# GFF feature type to associate tags with (gene, CDS, mRNA or such)
genetype=gene

if [ $# -eq 3 ]
then
    genetype=$3
fi

updates=no

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

# check if reference bfa does not exist, or is older than fasta

refbfa=${reffile%fasta}bfa

if needsUpdate $refbfa $reffile
then
    $MAQBIN fasta2bfa $reffile $refbfa
    updates=yes
fi

# count the number of reads in tagfile
# [TODO] save result for repetitive script use?
libsize=`grep -c ^@ $tagfile`

# check if bfq does not exist, or is older than fastq

# maq performs optimally around 2M reads, according to its manual pages
# if the number of reads is considerably larger than that (3M), split before bfq generation, 
# map separately and mapmerge the alignments.

maqmap=${tagfile%.fastq}_vs_${reffile%.fasta}.maqmap

if [ $libsize -gt 3000000 ] 
then
    tagfileprefix=${tagfile}.split
    maqmapupdate="no"
    
    # split into 8M/4=2M line files
    if needsUpdate $tagfile.done.split $tagfile
    then
	split -l 8000000 $tagfile ${tagfileprefix}
	echo ${tagfileprefix}[a-z][a-z] > $tagfile.done.split
	updates=yes
    fi

    splitnr=0

    bfqupdate="no"
    arglist=""
    for splittagfile in ${tagfileprefix}[a-z][a-z]
    do
	bfqfile=${splittagfile}.bfq
	if needsUpdate $bfqfile $splittagfile
	then
	    arglist=$arglist"$splittagfile $bfqfile "
	    bfqupdate="yes"
	fi
    done

    if [ "$bfqupdate" = "yes" ]
    then
	echo "DEBUG: run xargs -P $NPROC -n 2 $MAQBIN fastq2bfq ON $arglist."
	echo $arglist | xargs -P $NPROC -n 2 $MAQBIN fastq2bfq
	updates="yes"
    fi

    declare -a logs
    lognr=0

    splitmaqmapupdate="no"
    arglist=""

    for splittagfile in ${tagfileprefix}[a-z][a-z]
    do
	bfqfile=${splittagfile}.bfq

	splitmaqmap=${splittagfile}_vs_${reffile%.fasta}.maqmap
	splitmaqmaplog=$splitmaqmap.log

	# keep separate logs, then concatenate in the next concurrent section
	if needsUpdate $splitmaqmap $refbfa $bfqfile
	then
	    log=$splitmaqmaplog
	    logs[$lognr]=$log
	    lognr=$(( $lognr + 1 ))

	    rundate=`date`
	    echo "$MAQBIN map -n 3 $splitmaqmap $refbfa $bfqfile : (${rundate})" >> ${log}

	    arglist=$arglist"$splitmaqmap $refbfa $bfqfile ${log} "
	    splitmaqmapupdate="yes"
	fi

	splittagfiles[$splitnr]=$splittagfile
	splitmaqmap[$splitnr]=$splitmaqmap

	splitnr=$(( $splitnr + 1 ))
    done
    
    if [ "$splitmaqmapupdate" = "yes" ]
    then
	echo $arglist | xargs -P $NPROC -n 4 sh -c "$MAQBIN map -n 3 \$1 \$2 \$3 2>> \$4" -
	
	for mylog in ${logs[@]}
	do
	    cat $mylog >> ${maqmap}.log
	    rm $mylog
	done

	updates="yes"
	maqmapupdate="yes"
    fi

    if [ "$maqmapupdate" = "yes" ] 
    then
	$MAQBIN mapmerge $maqmap ${splitmaqmap[@]}
    fi
else
    bfqfile=${tagfile%fastq}bfq

    if needsUpdate $bfqfile $tagfile
    then
	$MAQBIN fastq2bfq $tagfile $bfqfile
	updates="yes"
    fi

    # map
    # size-sep to avoid length cutoff?
    if needsUpdate $maqmap $refbfa $bfqfile
    then
	log=${maqmap}.log
	rundate=`date`
	echo "$MAQBIN map -n 3 $maqmap $refbfa $bfqfile : (${rundate})" >>${log}
	$MAQBIN map -n 3 $maqmap $refbfa $bfqfile 2>>${log}
	updates=yes
    fi
fi

# convert binary maq map file to text format

if needsUpdate ${maqmap}.txt $maqmap
then
    $MAQBIN mapview ${maqmap} > ${maqmap}.txt
    updates=yes
fi

maqmapsingle=${maqmap}.q${alqt}

if needsUpdate ${maqmapsingle}.txt ${maqmap}.txt
then
    awk '($7>='$alqt') {print;}' < ${maqmap}.txt > ${maqmapsingle}.txt
    awk '($7==0) {print;}' < ${maqmap}.txt > ${maqmap}.multi.txt 
fi

# prepare for gff conversions (sizes needed for GBrowse GFF3 import)

refsizes=${reffile%fasta}sizes
if needsUpdate $refsizes $reffile $BINDIR/count_fasta_seq_len.pl
then
    $BINDIR/count_fasta_seq_len.pl < $reffile > $refsizes
    updates=yes
fi

lib=${tagfile%%.*} # first part of tagfile name is library name
maqmapgff=${maqmap}.gff
if needsUpdate $maqmapgff $maqmap ${maqmap}.txt $BINDIR/maqmaptxt2gff.pl
then
    $BINDIR/maqmaptxt2gff.pl -s $refsizes -l ${lib} -f maqmatch${lib} < ${maqmap}.txt > ${maqmapgff} ; 
    updates=yes
fi

maqmapsinglegff=${maqmapsingle}.gff
if [ "$single" = "yes" ] 
then
    if needsUpdate $maqmapsinglegff $maqmap ${maqmapsingle}.txt $BINDIR/maqmaptxt2gff.pl
    then

	$BINDIR/maqmaptxt2gff.pl -s $refsizes -l ${lib} -f maqq${alqt}${lib} < ${maqmapsingle}.txt > ${maqmapsinglegff};

	updates=yes
    fi
fi

# prepare for assigning more genome info to mapping
# most is executed on a per chromosome basis, and final results merged into a complete file

fastalist=${reffile%fasta}fastalist

if needsUpdate $fastalist $reffile 
then
    grep \> $reffile |sed -e 's/>//;'|cut -f 1 -d' ' > $fastalist
    updates=yes
fi

# prepare per-sequence part (chromosome or such) prerequisite files

for seq in `cat $fastalist` ;
do

    # update refgff split if there is a new refgff available

    refseqgff=${reffile%.fasta}.${seq}.gff

    grepsafeseq=`echo ${seq}|sed -e 's/|/\\\|/g'`

    if needsUpdate ${refseqgff} ${refgff}
    then
	grep "${grepsafeseq}" ${refgff} > ${refseqgff}
	updates=yes
    fi

    # update ref fasta split if there is a new reffile available

    refseq=${reffile%.fasta}.${seq}.fasta 

    if needsUpdate ${refseq} ${reffile} $BINDIR/fasta_header_grep.pl
    then
	$BINDIR/fasta_header_grep.pl "${grepsafeseq}" ${reffile} > ${refseq}
	updates=yes
    fi

done

# produce main gff files

maqmaptab=${maqmap}.counts.tab
maqmapsingletab=${maqmapsingle}.counts.tab

# update flags for agglomerate multi sequence unit files 
maqmap_update=no
maqmapsingle_update=no
crunch_update=no
crunchsingle_update=no

for seq in `cat $fastalist` ;
do
    # split the maqmap data (gff format) into per sequence unit files, full set of mappers

    seqmaqmapgff=${maqmapgff%.gff}.${seq}.gff
    
#    grepsafeseq=`echo ${seq}|sed -e 's/|/\\\|/g'`

    if needsUpdate $seqmaqmapgff $maqmapgff
    then
	grep ${seq} ${maqmapgff} > ${seqmaqmapgff}
	updates=yes
    fi

    # split the maqmap data (gff format) into per sequence unit files, single mappers

    seqmaqmapsinglegff=${maqmapsinglegff%.gff}.${seq}.gff

    if [ "$single" = "yes" ]
    then
	if needsUpdate $seqmaqmapsinglegff $maqmapsinglegff
	then
	    grep ${seq} ${maqmapsinglegff} > ${seqmaqmapsinglegff}
	    updates=yes
	fi
    fi

    refseq=${reffile%.fasta}.${seq}.fasta

    # convert initial map gff's into naïve 1-space wiggle files for gbrowse display
    
    seqmaqmapwig=${seqmaqmapgff%.gff}.wig
    if needsUpdate $seqmaqmapwig $seqmaqmapgff $refseq $BINDIR/gff2wig.pl 
    then
	log=${maqmap}.log
	rundate=`date`
	echo "$seq : gff2wig ($rundate)" >> $log
	$BINDIR/gff2wig.pl -s $refseq -t $seqmaqmapgff -f maqmatch${lib} -o $seqmaqmapwig >> $log
	updates=yes
    fi

    if [ "$single" = "yes" ]
    then
	seqmaqmapsinglewig=${seqmaqmapsinglegff%.gff}.wig
	if needsUpdate $seqmaqmapsinglewig $seqmaqmapsinglegff $refseq $BINDIR/gff2wig.pl 
	then
	    log=${maqmap}.log
	    rundate=`date`
	    echo "$seq : gff2wig ($rundate)" >> $log
	    $BINDIR/gff2wig.pl -s $refseq -t $seqmaqmapsinglegff -f maqq${alqt}${lib} -o $seqmaqmapsinglewig >> $log
	    updates=yes
	fi
    fi
    
    # the initial mapped read to gene connection in the form of gene.gffs, all mappers

    refseqgff=${reffile%.fasta}.${seq}.gff
    genegff=${maqmapgff%.gff}.${seq}.gene.gff

    if needsUpdate $genegff $refseqgff $refseq $seqmaqmapgff "$BINDIR/per_gene_coverage_info.pl"
    then
	log=${maqmap}.log
	rundate=`date`
	echo "$seq : per_gene_coverage_info ($rundate)" >> $log
	$BINDIR/per_gene_coverage_info.pl -g ${refseqgff} -s ${refseq} -t ${seqmaqmapgff} -f maqmatch${lib} -o ${genegff} -c ${genetype} -N ${libsize} >> $log
	genegff_update=yes
	updates=yes
    fi

    # initial gene-gffs, single mappers

    genesinglegff=${maqmapsinglegff%.gff}.${seq}.gene.gff
    if [ "$single" = "yes" ]
    then
	if needsUpdate $genesinglegff $refseqgff $refseq $seqmaqmapsinglegff "$BINDIR/per_gene_coverage_info.pl"
	then
	    log=${maqmap}.log
	    rundate=`date`
	    echo "$seq : per_gene_coverage_info ($rundate) : single mappers : " >> $log
	    $BINDIR/per_gene_coverage_info.pl -g ${refseqgff} -s ${refseq} -t ${seqmaqmapsinglegff} -f maqq${alqt}${lib} -o ${genesinglegff} -c ${genetype} -N ${libsize} >> $log
	    genesinglegff_update=yes
	    updates=yes
	fi
    fi

    # antisense tags

    astaggff=${maqmapgff%.gff}.${seq}.astags.gff
    if needsUpdate $astaggff $refseqgff $refseq $seqmaqmapgff "$BINDIR/list_antisense_splicesites.pl"
    then
	aslog=${maqmap}.antisense.log
	rundate=`date`
	echo "$seq : list_antisense_splicesites ($rundate)" >> $aslog
	$BINDIR/list_antisense_splicesites.pl -g ${refseqgff} -s ${refseq} -t ${seqmaqmapgff} -f maqmatch${lib} -o ${astaggff} -c ${genetype} >> $aslog
	updates=yes
    fi

    # antisense tags, single mappers

    if [ "$single" = "yes" ]
    then
	astagsinglegff=${maqmapsinglegff%.gff}.${seq}.astags.gff
	if needsUpdate $astagsinglegff $refseqgff $refseq $seqmaqmapsinglegff "$BINDIR/list_antisense_splicesites.pl"
	then
	    aslogsingle=${maqmapsingle}.antisense.log
	    rundate=`date`
	    echo "$seq : list_antisense_splicesites ($rundate)" >> $aslogsingle
	    $BINDIR/list_antisense_splicesites.pl -g ${refseqgff} -s ${refseq} -t ${seqmaqmapsinglegff} -f maqq${alqt}${lib} -o ${astagsinglegff} -c ${genetype} >> $aslogsingle
	    updates=yes
	fi
    fi
done

# join split GFFs

# additional update criteria
allgenegff=${maqmapgff%.gff}.gene.gff
if needsUpdate $allgenegff $refgff $refseq $maqmapgff "$BINDIR/per_gene_coverage_info.pl"
then
    gff_update=yes
fi 
    
if [ "${gff_update}" = "yes" ] 
then
    cat /dev/null > $allgenegff.tmp

    for seq in `cat $fastalist` ;
    do
	genegff=${maqmapgff%.gff}.${seq}.gene.gff
	cat $genegff >> $allgenegff.tmp
    done

# will multiple ##gff-version 3 lines cause a problem? maybe, maybe not.
    grep gff-version $allgenegff.tmp |sort |uniq > $allgenegff
    grep -v gff-version $allgenegff.tmp >> $allgenegff
    rm $allgenegff.tmp

    updates=yes
fi

allgenesinglegff=${maqmapsinglegff%.gff}.gene.gff
if needsUpdate $allgenesinglegff $refgff $refseq $maqmapsinglegff "$BINDIR/per_gene_coverage_info.pl"
then
    gffsingle_update=yes
fi

if [ "${gffsingle_update}" = "yes" ] 
then
    cat /dev/null > $allgenesinglegff.tmp
    for seq in `cat $fastalist` ;
    do
	genesinglegff=${maqmapsinglegff%.gff}.${seq}.gene.gff
	cat $genesinglegff >> $allgenesinglegff.tmp
    done

# will multiple ##gff-version 3 lines cause a problem?
    grep gff-version $allgenesinglegff.tmp |sort |uniq > $allgenesinglegff
    grep -v gff-version $allgenesinglegff.tmp >> $allgenesinglegff
    rm $allgenesinglegff.tmp

    updates=yes
fi

# run rescoring routine based on complete genome (gene lengths and total number of mapped tags could be made available, but for per gene nuc fractions, this is probably easiest?
# or just work on the tabs? would be ok for now, but eventually you may want those gffs.

allgenemapinfo=${allgenegff%.gff}.info
if needsUpdate $allgenemapinfo $allgenegff "$BINDIR/get_globals_from_tagged_rnaseq_gff.pl"
then
    log=${maqmap}.log
    rundate=`date`
    
    $BINDIR/get_globals_from_tagged_rnaseq_gff.pl -g ${allgenegff} -o ${allgenemapinfo} -c ${genetype} -N ${libsize} >> $log
    updates=yes
fi

allgenesinglemapinfo=${allgenesinglegff%.gff}.info
if needsUpdate $allgenesinglemapinfo $allgenesinglegff "$BINDIR/get_globals_from_tagged_rnaseq_gff.pl"
then
    log=${maqmap}.log
    rundate=`date`
    
    $BINDIR/get_globals_from_tagged_rnaseq_gff.pl -g ${allgenesinglegff} -o ${allgenesinglemapinfo} -c ${genetype} -N ${libsize} >> $log

    genegff_update=yes
fi

totalmappedtagcount=`grep mapped ${allgenemapinfo} |cut -f2 -d:`
totalnuovergenelength=`grep nu ${allgenemapinfo} |cut -f2 -d:`

totalmappedsingletagcount=`grep mapped ${allgenesinglemapinfo} |cut -f2 -d:`
totalsinglenuovergenelength=`grep nu ${allgenesinglemapinfo} |cut -f2 -d:`

# correct split gffs
# ah, or make this in the tab file generation directly, and have a separate tagged_rnaseq_gff_to_counts.pl that makes gffs? same (bad) as before..

# anyway, make new counts tabs in a per fasta fashion

for seq in `cat $fastalist` ;
do
    refseqgff=${reffile%.fasta}.${seq}.gff
    refseq=${reffile%.fasta}.${seq}.fasta
    
    genegff=${maqmapgff%.gff}.${seq}.gene.gff
    
    # counts tabs, reporting most of the per gene mapping data in a tabular format, for all mapping reads
    seqmaqmaptab=${maqmaptab%.tab}.${seq}.tab

    if needsUpdate $seqmaqmaptab $genegff $refseq $allgenemapinfo $BINDIR/tagged_rnaseq_gff_to_counts_tab.pl
    then
	$BINDIR/tagged_rnaseq_gff_to_counts_tab.pl -g ${genegff} -s ${refseq} -f maqmatch${lib} -N $libsize -c ${genetype} -o ${seqmaqmaptab} -T $totalmappedtagcount -n $totalnuovergenelength
	maqmap_update=yes
	updates=yes
    fi
   
    # counts tabs, single mappers

    genesinglegff=${maqmapsinglegff%.gff}.${seq}.gene.gff
    seqmaqmapsingletab=${maqmapsingletab%.tab}.${seq}.tab
    if [ "$single" = "yes" ] 
    then
	if needsUpdate $seqmaqmapsingletab $genesinglegff $refseq $allgenesinglemapinfo $BINDIR/tagged_rnaseq_gff_to_counts_tab.pl
	then
	    $BINDIR/tagged_rnaseq_gff_to_counts_tab.pl -g ${genesinglegff} -s ${refseq} -f maqq${alqt}${lib} -c ${genetype} -o ${seqmaqmapsingletab} -T $totalmappedsingletagcount -n $totalsinglenuovergenelength
	    maqmapsingle_update=yes
	    updates=yes
	fi
    fi
done

# check update flags for any need to agglomerate multi sequence unit files

# join the tabs

function UpdateTab()
{
    local this_tab=$1
    local suffix=$2

    cat /dev/null > ${this_tab}

    for seq in `cat $fastalist` ;
    do
	local this_seqtab=${this_tab%%.${suffix}}.${seq}.${suffix}
    
	echo $seq >> $this_tab
	sort -k2n,2n ${this_seqtab} >> ${this_tab}
    done

    updates=yes
}

# update the joint tabs also if any of the subtabs has mysteriously updated.. (shortcut for development)
for seq in `cat $fastalist` ;
do
    # check counts tab seqtabs
    seqmaqmaptab=${maqmaptab%.tab}.${seq}.tab

    if needsUpdate $maqmaptab $seqmaqmaptab
    then
	maqmap_update=yes
    fi

    # check single mapper counts tab seqtabs    
    if [ "$single" = "yes" ]
    then
	seqmaqmapsingletab=${maqmapsingletab%.tab}.${seq}.tab

	if needsUpdate $maqmapsingletab $seqmaqmapsingletab
	then
	    maqmapsingle_update=yes
	fi
    fi    
done

if [ "${maqmap_update}" = "yes" ] 
then
    UpdateTab ${maqmaptab} tab
fi

if [ "${maqmapsingle_update}" = "yes" ] 
then
    UpdateTab ${maqmapsingletab} tab
fi

# function to make filtered sub-tables from the big ones

summary=${maqmaptab%.tab}.summary

if needsUpdate $summary $maqmap.txt $maqmaptab ${maqmapsingletab} ${maqmapsingle}.txt $BINDIR/median.pl
then
    echo -n > $summary

    # produce a nice summary stats file..
    echo "${lib} Number of reads: $libsize" >> $summary

    echo -n "${lib} Number of reads aligned ok: " >> $summary
    wc -l $maqmap.txt |awk '{ print "\t",$1 }' >> $summary

    # percent aligned ok.. 

    echo -n "${lib} number of maq singlemappers, q >= ${alqt}: " >> $summary
    wc -l ${maqmapsingle}.txt |awk '{ print "\t",$1 }' >> $summary

    echo -n "${lib} Number of maq multimappers, q == 0: " >> $summary
    wc -l ${maqmap}.multi.txt |awk '{ print "\t",$1 }' >> $summary

    # use a perl script for the median calculations 

    echo -n "${lib} Maq single mappers alignment quality: " >> $summary
    awk '($7>0) {print $7;}' ${maqmap}.txt | $BINDIR/median.pl >> $summary

    echo -n "${lib} Maq single mappers mismatches: " >> $summary
    awk '($7>0) {print $10;}' ${maqmap}.txt | $BINDIR/median.pl >> $summary
    
    echo -n "${lib} Maq all mapped reads mismatches: " >> $summary
    awk '{print $10;}' ${maqmap}.txt | $BINDIR/median.pl >> $summary
    
    echo -n "${lib} ">> $summary
    echo -n "Genes with counts">> $summary
    awk 'BEGIN {tick = 0} ($4 > 0) {tick=tick+1} END {print " ",tick}' $maqmaptab >> $summary;
 
    topcounts=${maqmaptab%tab}topcounts

    awk 'BEGIN { print "Name","Count","Desc"; } ($5 >= 0) { print $1,$5,$13; }' $maqmaptab > $topcounts

    cat > ${maqmap}.tags.R <<-RTAGS
	lib<-read.table("${topcounts}", header = TRUE, sep = " ");
        
	meancount <- mean(lib\$Count)
	meanstr <- sprintf("Mean gene tag count, $lib: %.2f",meancount)
	cat(meanstr,file="${summary}", sep="\\n", append=TRUE)

#	mediancount <- median(lib\$Count)
#	medianstr <- sprintf("Median gene tag count, $lib $liblabel: %.1f",mediancount)
#	cat(medianstr,file="${summary}", sep="\\n", append=TRUE)

	f<-fivenum(lib\$Count)
	fstr <- sprintf("$lib tag count for genes min %d, 1st quartile %d, median %d, 3rd quartile %d and max %d.",f[1],f[2],f[3],f[4],f[5])
	cat(fstr,file="${summary}", sep="\\n", append=TRUE)

	RTAGS
    
    R CMD BATCH --no-save --slave ${maqmap}.tags.R /dev/null

    showtoptagged=40
    echo >> $summary
    echo "Genes with the $showtoptagged highest tag counts" >> $summary    
    head -1 $topcounts >> $summary
    grep -v ^Name $topcounts |sort -k2,2nr | head -$showtoptagged >> $summary

    updates=yes
fi

if [ "$updates" = "no" ]
then
    log=${maqmap}.log
    rundate=`date`
    echo "Project was already up to date. Nothing to do. Please update date on appropriate input/program/intermediate results files to force results updates. ($rundate)" | tee -a $log
else
    log=${maqmap}.log
    rundate=`date`
    echo "Project has been brought up to date. ($rundate)" | tee -a $log
fi
exit
