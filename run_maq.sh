#!/bin/bash
#
# Daniel Nilsson, 2009
#
# daniel.nilsson@izb.unibe.ch, daniel.k.nilsson@gmail.com
#
# POD documentation to follow throughout the file - use e.g. perldoc to read
#

: <<'POD_INIT'

=head1 NAME

run_maq.sh - core of the SLadd pipeline for turning SLT data into results

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@izb.unibe.ch, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2009, 2010 held by Daniel Nilsson. The package is realesed for use under the Perl Artistic License.

=head1 SYNOPSIS

USAGE: C<run_maq.sh tagfile.fastq referenecesequence.fasta [gene_feature_type(gene)]>

=head1 DESCRIPTION

Takes SLT tags and an annotated reference genome and provides
annotation of spliced leader addition sites and expression levels.

This is a make like pipeline implemented in bash, calling mostly perl
and R programs. The need for (re)evaluation of a particular result is
based on the existance and modification time of the result file in
relation to the files needed to produce it. If you want to repeat an
analysis, ensure that your new/modified input files are dated more
recently than the particular analysis result you are interested in. Or
invoke the pipeline with forceupdate=yes, but keep in mind that that
will rerun everything.

Please see the F<sladd_howto.pdf>/F<sladd_howto.odt> for more information.

=head1 DEPENDENCIES

Expects:

=over 4

=item * 

F<referenceseqeunce.gff> in same dir as F<referencesequence.fasta>

=item *
          
maq binary installed in C<$MAQBIN>

=item *

Several of the other "sladdlampa"-scripts in C<$BINDIR>

=item *

Indirectly, bioPerl in the perl library path (and a perl 5.8.x interpreter at /usr/bin)

=back

=head1 SHELL ENVIRONMENT VARIABLES

Several environment variables influence the pipeline behaviour.

=over 4

=item BINDIR [path (~/install/sladd)]          

Directory where the rest of the SLADD pipeline lives.

=item ALIGN [<maq|bowtie> (maq)]

The short read aligner to use.  Currently limited to maq or bowtie.

=item MAQBIN [path (~/install/maq-0.7.1/maq)] 

Your maq binary.

=item BOWTIEBINDIR [path (~/install/bowtie-0.12.5/)]

The path to your bowtie installation.

=item SEQLOGOBIN [path (~/install/weblogo/seqlogo)]
                 
The path to the Weblogo seqlogo binary

=item single [<yes|no> (yes)]
             
Toggle runnig separate analysis with alignment quality cutoff for trusted (single) mappings.                    

=item alqt [integer (30)] 
           
Alignment quality cutoff for trusted (single) mappings.
                     
=item forceupdate [<yes|no> (no)]
                     
Force rerun of all analysis regardless of 
modification times.

=item runuorf [<yes|no> (no)]

Run the fairly slow variant of uORF analysis.

=item utrcutoff [integer (0)] 

Cutoff UTR length for splicesite analysis.
Set to the longest UTR you consider meaningful.
Set to 0 for no limit.

=item uorfutrcutofflen [integer (2000)]
          
Cutoff UTR length for uORF analysis.

=item plotsplicesites [<yes|no> (yes)]

Set plotsplicesites=no to avoid extracting and plotting splicesites.

=back

=head1 OPTIONS AND ARGUMENTS

=over 4

=cut
 
POD_INIT

#
# check input environment variables and set unset ones to default values
# 

# BINDIR is where the rest of the scripts and binaries in the SLadd pipeline reside
if [ -z "$BINDIR" ]
then
	#BINDIR=~/src/sladd
	BINDIR=~/install/sladd
fi

export PERL5LIB=$PERL5LIB:$BINDIR

# Set prefered short read aligner
if [ -z "$ALIGN" ]
then
    ALIGN=maq
    #ALIGN=bowtie
fi

# MAQBIN is the MAQ installation directory (path to the maq binary)
if [ -z "$MAQBIN" ]
then
	MAQBIN=/usr/bin/maq
	#MAQBIN=~/install/maq-0.7.1/maq
fi

# BOWTIEBINDIR is the bowtie installation directory
if [ -z "$BOWTIEBINDIR" ]
then
	BOWTIEBINDIR=~/install/bowtie-0.12.5
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

# set plotsplicesites=no to avoid extracting and plotting splicesites
if [ -z "$plotsplicesites" ]
then
    plotsplicesites=yes
fi

# splicesite plot window, bases before and after the 3'sa dinucleotide (typically AG)
beforeag=88
afterag=5

# utrlen to consider for splicesite analysis
if [ -z "$utrcutoff" ]
then
    utrcutoff=0
fi

# utrlen to consider for uorf analysis
if [ -z "$uorfutrcutofflen" ]
then
    uorfutrcutofflen=2000
fi

# actually run uorf check (takes quite some processing power on large chrs with many splicesites..)
if [ -z "$runuorf" ]
then
    runuorf=no
fi

# uses pipelinefunk.sh for needsUpdate, registerFile etc.

. $BINDIR/pipelinefunk.sh

# CALLED will contain a complete pathname for this script
CALLED=$0

# command line arguments
if [ $# -lt 2 ]
then
	perldoc $CALLED
        echo "USAGE: ${CALLED##*/} tagfile.fastq refsequence.fasta [gene_gff_type]"
        exit 1
fi

: <<'POD_ARG'

=item first argument: tagfilename.fastq

first argument : F<tagfilename.fastq> : name convention is F<libname.sizefilter.fastq>

=cut

POD_ARG

tagfile=$1

: <<'POD_ARG'

=item second argument: referencesequence.fasta

second argument : F<refsequencename.fasta>, also expects to find F<refsequencename.gff> in the same directory.

=cut

POD_ARG

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

: <<'POD_ARG'

=item third argument, optional: genetype

third argument: genetype : GFF feature type to associate tags with (gene, CDS, mRNA or such).

=back

=cut

POD_ARG

genetype=gene

if [ $# -eq 3 ]
then
    genetype=$3
fi

updates=no
: <<'POD_FILE'

=head1 FILES AND FUNCTIONS

=over 4

=cut 

POD_FILE

# check if reference bfa does not exist, or is older than fasta

if [ "$ALIGN" == "maq" ] 
then
    refbfa=${reffile%fasta}bfa

    if needsUpdate $refbfa $reffile
    then
	$MAQBIN fasta2bfa $reffile $refbfa
	registerFile $refbfa temp
	updates=yes
    fi
elif [ "$ALIGN" == "bowtie" ]
then
    refbuild=${reffile%fasta}bowtiebuild
    
    if needsUpdate $refbuild $reffile
    then
	$BOWTIEBINDIR/bowtie-build $reffile $refbuild
	registerFile $refbuild temp
	updates=yes
    fi
else
    echo "ERROR: sorry, run_maq.sh is not configured to run aligner ${ALIGN}. Consider using another value for \$ALIGN."
    exit 1
fi

# count the number of reads in tagfile
# [TODO] save result for repetitive script use?

libsize=`grep -c ^@ $tagfile`

map=${tagfile%.fastq}_vs_${reffile%.fasta}.${ALIGN}map

if [ "$ALIGN" == "maq" ]
then 

    # maq performs optimally around 2M reads, according to its manual pages
    # if the number of reads is considerably larger than that (3M), split before bfq generation, 
    # map separately and mapmerge the alignments.

    if [ $libsize -gt 3000000 ] 
    then
	tagfileprefix=${tagfile}.split
	mapupdate="no"
	
        # split into 8M/4=2M line files

	if needsUpdate $tagfile.done.split $tagfile
	then
	    split -l 8000000 $tagfile ${tagfileprefix}
	    echo ${tagfileprefix}[a-z][a-z] > $tagfile.done.split 
	    registerFile $tagfile.done.split temp
	fi

	splitnr=0

    # check if bfq does not exist, or is older than fastq
	bfqupdate="no"
	arglist=""
	for splittagfile in ${tagfileprefix}[a-z][a-z]
	do
	    #now register the split tagfile here, as we didn't have the name direcly available in the last loop
	    registerFile $splittagfile temp

	    bfqfile=${splittagfile}.bfq
	    if needsUpdate $bfqfile $splittagfile
	    then
		arglist=$arglist"$splittagfile $bfqfile "
		registerFile $bfqfile temp
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
	    
	splitmapupdate="no"
	arglist=""

	for splittagfile in ${tagfileprefix}[a-z][a-z]
	do
	    bfqfile=${splittagfile}.bfq
	    splitmap=${splittagfile}_vs_${reffile%.fasta}.${ALIGN}map
	    splitmaplog=${splitmap}.log

	    if needsUpdate $splitmap $refbfa $bfqfile
	    then
		log=${splitmaplog}
		logs[$lognr]=$log
		lognr=$(( $lognr + 1 ))
	    
		rundate=`date`
		echo "$MAQBIN map -n 3 $splitmap $refbfa $bfqfile : (${rundate})" >> ${log}

		arglist=$arglist"$splitmap $refbfa $bfqfile ${log} "
		splitmapupdate="yes"

		registerFile $splitmap temp
	    fi
	    
	    splittagfiles[$splitnr]=$splittagfile
	    splitmap[$splitnr]=$splitmap
	    splitnr=$(( $splitnr + 1 ))
	done

	if [ "$splitmapupdate" = "yes" ]
	then
	    echo $arglist | xargs -P $NPROC -n 4 sh -c "$MAQBIN map -n 3 \$1 \$2 \$3 2>> \$4" -
	
	    for mylog in ${logs[@]}
	    do
		cat $mylog >> ${map}.log
		rm $mylog
	    done

	    updates="yes"
	    mapupdate="yes"
	fi

	if [ "$mapupdate" = "yes" ] 
	then
	    $MAQBIN mapmerge $map ${splitmap[@]}
	    registerFile $map result
	fi
    else

    # check if bfq does not exist, or is older than fastq

	bfqfile=${tagfile%fastq}bfq

	if needsUpdate $bfqfile $tagfile
	then
	    $MAQBIN fastq2bfq $tagfile $bfqfile
	    registerFile $bfqfile temp	    
	    updates=yes
	fi

    # map
    # size-sep to avoid length cutoff?
	if needsUpdate $map $refbfa $bfqfile
	then
	    log=${map}.log
	    rundate=`date`
	    echo "map -n 3 $map $refbfa $bfqfile : (${rundate})" >>${log}
	    $MAQBIN map -n 3 $map $refbfa $bfqfile 2>>${log}
	    updates=yes
	fi
    fi
elif [ "$ALIGN" == "bowtie" ]
then
    log=${map}.log

    rundate=`date`
    echo "bowtie --best -M 1 -n 3 -p $NCPU $refhash $tagfile $map : (${rundate})" >> ${log}
    # -M 1 to account for multimappers..
    $BOWTIEBINDIR/bowtie --best -M 1 -n 3 -p $NCPU $refhash $tagfile $map
    registerFile $map result

fi
# no else for method choice -- this will have been cleared previously.

# convert binary maq map file to text format
if [ "$ALIGN" == "maq" ]
then
    if needsUpdate ${map}.txt $map
    then
	$MAQBIN mapview ${map} > ${map}.txt
	registerFile ${map}.txt temp
	updates=yes
    fi
fi

single=q${alqt}
if [ "$ALIGN" == "bowtie" ]
then 
    single=single
fi

mapsingle=${map}.${single}

#if [ "$ALIGN" == "bowtie" ]
#then
#    log=${map}.log
#
#    rundate=`date`
#    echo "bowtie --best -n 3 -m 1 -p $NCPU $refhash $tagfile $map : (${rundate})" >> ${log}
#    $BOWTIEBINDIR/bowtie --best -M 1 -n 3 -p $NCPU $refhash $tagfile $mapsingle
#fi

if [ "$ALIGN" == "bowtie" ]
then
	awk '($7==0) {print;}' < ${map} > ${mapsingle}.txt 
	registerFile ${mapsingle}.txt temp
	awk '($7>=2) {print;}' < ${map} > ${map}.multi.txt
	registerFile ${map}.multi.txt temp
fi

if [ "$ALIGN" == "maq" ]
then
    if needsUpdate ${mapsingle}.txt ${map}.txt
    then
	awk '($7>='$alqt') {print;}' < ${map}.txt > ${mapsingle}.txt
	registerFile ${mapsingle}.txt temp
	awk '($7==0) {print;}' < ${map}.txt > ${map}.multi.txt 
	registerFile ${map}.multi.txt temp
    fi
fi

# prepare for gff conversions (sizes needed for GBrowse GFF3 import)

refsizes=${reffile%fasta}sizes
if needsUpdate $refsizes $reffile $BINDIR/count_fasta_seq_len.pl
then
    $BINDIR/count_fasta_seq_len.pl < $reffile > $refsizes
    registerFile $refsizes temp
    updates=yes
fi

lib=${tagfile%%.*} # first part of tagfile name is library name
mapgff=${map}.gff

if [ "$ALIGN" == "maq" ]
then
    if needsUpdate $mapgff $map ${map}.txt $BINDIR/maqmaptxt2gff.pl
    then
	$BINDIR/maqmaptxt2gff.pl -s $refsizes -l ${lib} -f ${ALIGN}match${lib} < ${map}.txt > ${mapgff};
	registerFile $mapgff result
	updates=yes
    fi
elif [ "$ALIGN" == "bowtie" ]
then
    if needsUpdate $mapgff $map $BINDIR/bowtie_to_gff.pl
    then
	$BINDIR/bowtie_to_gff.pl -s $refsizes -l ${lib} -f ${ALIGN}match${lib} < ${map} > ${mapgff};
	registerFile $mapgff result
	updates=yes
    fi
fi

mapsinglegff=${mapsingle}.gff
if [ "$single" = "yes" ]
then
    if [ "$ALIGN" == "maq" ]
    then
	if needsUpdate $mapsinglegff $map ${mapsingle}.txt $BINDIR/maqmaptxt2gff.pl
	then
	    $BINDIR/maqmaptxt2gff.pl -s $refsizes -l ${lib} -f ${ALIGN}${single}${lib} < ${mapsingle}.txt > ${mapsinglegff};
	    registerFile ${mapsinglegff} result
	    updates=yes
	fi
    elif [ "$ALIGN" == "bowtie" ]
    then
	if needsUpdate $mapsinglegff $mapsingle $BINDIR/bowtie_to_gff.pl
	then
	    $BINDIR/bowtie_to_gff.pl -s $refsizes -l ${lib} -f ${ALIGN}${single}${lib} < ${mapsingle} > ${mapsinglegff};
	    registerFile ${mapsinglegff} result
	    updates=yes
	fi
    fi
fi

# prepare for assigning more genome info to mapping
# most is executed on a per chromosome basis, and final results merged into a complete file

fastalist=${reffile%fasta}fastalist

if needsUpdate $fastalist $reffile 
then
    grep \> $reffile |sed -e 's/>//;'|cut -f 1 -d' ' > $fastalist
    registerFile ${fastalist} temp
fi

# prepare per-sequence part (chromosome or such) prerequisite files

for seq in `cat $fastalist` ;
do
    # update refseqgff split if there is a new refseqgff available

    refseqgff=${reffile%.fasta}.${seq}.gff

    grepsafeseq=`echo ${seq}|sed -e 's/|/\\\|/g'`

    if needsUpdate ${refseqgff} ${refgff}
    then
	grep "${grepsafeseq}" ${refgff} > ${refseqgff}
	registerFile ${refseqgff} temp
    fi

    # update ref fasta split if there is a new reffile available

    refseq=${reffile%.fasta}.${seq}.fasta 

    if needsUpdate ${refseq} ${reffile} $BINDIR/fasta_header_grep.pl
    then
	$BINDIR/fasta_header_grep.pl "${grepsafeseq}" ${reffile} > ${refseq}
	registerFile ${refseq} temp
    fi

done

# produce main gff files

maptab=${map}.counts.tab
mapsingletab=${mapsingle}.counts.tab

mapuorftab=${maptab%.tab}.uorf.tab
mapsingleuorftab=${mapsingletab%.tab}.uorf.tab

# update flags for agglomeration of multi sequence unit files 
map_update=no
mapsingle_update=no
crunch_update=no
crunchsingle_update=no
uorf_update=no
uorfsingle_update=no

for seq in `cat $fastalist` ;
do
    # split the map data (gff format) into per sequence unit files, full set of mappers

    seqmapgff=${mapgff%.gff}.${seq}.gff
    
#    grepsafeseq=`echo ${seq}|sed -e 's/|/\\\|/g'`

    if needsUpdate $seqmapgff $mapgff
    then
	grep ${seq} ${mapgff} > ${seqmapgff}
	registerFile $seqmapgff temp
    fi

    # split the map data (gff format) into per sequence unit files, single mappers

    seqmapsinglegff=${mapsinglegff%.gff}.${seq}.gff

    if [ "$single" = "yes" ]
    then
	if needsUpdate $seqmapsinglegff $mapsinglegff
	then
	    grep ${seq} ${mapsinglegff} > ${seqmapsinglegff}
	    registerFile $seqmapsinglegff temp
	fi
    fi

    # the initial mapped read to gene connection in the form of gene.gffs, all mappers

    refseqgff=${reffile%.fasta}.${seq}.gff
    refseq=${reffile%.fasta}.${seq}.fasta
    
    genegff=${mapgff%.gff}.${seq}.gene.gff

    if needsUpdate $genegff $refseqgff $refseq $seqmapgff "$BINDIR/per_gene_splicesite_info.pl"
    then
	log=${map}.log
	rundate=`date`
	echo "$seq : per_gene_splicesite_info ($rundate)" >> $log
	$BINDIR/per_gene_splicesite_info.pl -g ${refseqgff} -s ${refseq} -t ${seqmapgff} -f ${ALIGN}match${lib} -o ${genegff} -c ${genetype} -U ${utrcutoff} >> $log
	registerFile $genegff result
    fi

    # initial gene-gffs, single mappers

    genesinglegff=${mapsinglegff%.gff}.${seq}.gene.gff
    if [ "$single" = "yes" ]
    then
	if needsUpdate $genesinglegff $refseqgff $refseq $seqmapsinglegff "$BINDIR/per_gene_splicesite_info.pl"
	then
	    log=${map}.log
	    rundate=`date`
	    echo "$seq : per_gene_splicesite_info ($rundate) : single mappers : " >> $log
	    $BINDIR/per_gene_splicesite_info.pl -g ${refseqgff} -s ${refseq} -t ${seqmapsinglegff} -f ${ALIGN}${single}${lib} -o ${genesinglegff} -c ${genetype} -U ${utrcutoff} >> $log
	    registerFile $genesinglegff result
	fi
    fi

    # antisense tags

    astaggff=${mapgff%.gff}.${seq}.astags.gff
    if needsUpdate $astaggff $refseqgff $refseq $seqmapgff "$BINDIR/list_antisense_splicesites.pl"
    then
	aslog=${map}.antisense.log
	rundate=`date`
	echo "$seq : list_antisense_splicesites ($rundate)" >> $aslog
	$BINDIR/list_antisense_splicesites.pl -g ${refseqgff} -s ${refseq} -t ${seqmapgff} -f ${ALIGN}match${lib} -o ${astaggff} -c ${genetype} >> $aslog
	registerFile $astaggff result
    fi
    # why not crunch these when we crunch the singlemappers?

    # antisense tags, single mappers

    astagsinglegff=${mapsinglegff%.gff}.${seq}.astags.gff
    if needsUpdate $astagsinglegff $refseqgff $refseq $seqmapsinglegff "$BINDIR/list_antisense_splicesites.pl"
    then
	aslogsingle=${mapsingle}.antisense.log
	rundate=`date`
	echo "$seq : list_antisense_splicesites ($rundate)" >> $aslogsingle
	$BINDIR/list_antisense_splicesites.pl -g ${refseqgff} -s ${refseq} -t ${seqmapsinglegff} -f ${ALIGN}${single}${lib} -o ${astagsinglegff} -c ${genetype} >> $aslogsingle
	registerFile $astagsinglegff temp
    fi

    astagsinglecrunchgff=${mapsinglegff%.gff}.${seq}.astags.crunch.gff
    if needsUpdate $astagsinglecrunchgff $mapsingle $astagsinglegff $BINDIR/crunch_gff_with_counts.pl
    then
	$BINDIR/crunch_gff_with_counts.pl -g ${astagsinglegff} -s ${refseq} -f ${ALIGN}${single}${lib} -N $libsize -c ${genetype} -o ${astagsinglecrunchgff}
	regsiterFile $astagsinglecrunchgff result
#	crunchsingle_update=yes
	updates=yes
    fi

    # counts tabs, reporting most of the per gene mapping data in a tabular format, for all mapping reads

    seqmaptab=${maptab%.tab}.${seq}.tab

    if needsUpdate $seqmaptab $genegff $refseq $BINDIR/tagged_gff_to_counts_tab.pl
    then
	$BINDIR/tagged_gff_to_counts_tab.pl -g ${genegff} -s ${refseq} -f ${ALIGN}match${lib} -N $libsize -c ${genetype} -o ${seqmaptab}
	registerFile $seqmaptab temp
	map_update=yes
	updates=yes
    fi
   
    # counts tabs, single mappers

    seqmapsingletab=${mapsingletab%.tab}.${seq}.tab
    if [ "$single" = "yes" ] 
    then
	if needsUpdate $seqmapsingletab $genesinglegff $refseq $BINDIR/tagged_gff_to_counts_tab.pl
	then
	    $BINDIR/tagged_gff_to_counts_tab.pl -g ${genesinglegff} -s ${refseq} -f ${ALIGN}${single}${lib} -N $libsize -c ${genetype} -o ${seqmapsingletab}
	    registerFile $genesingletab temp
	    mapsingle_update=yes
	    updates=yes
	fi
    fi

    # [TODO] not sure that we're actually using this now? useful for loading completed gene features onto gbrowse? fix multi-utr-len entry?!
    # all mappers
    seqmapcountsgff=${seqmapgff%gff}counts.gff
    if needsUpdate ${seqmapcountsgff} $genegff $refseq $BINDIR/tagged_gff_to_counts.pl 
    then
	$BINDIR/tagged_gff_to_counts.pl -g ${genegff} -s ${refseq} -f ${ALIGN}match${lib} -N $libsize -c ${genetype} -o ${seqmapcountsgff}
	registerFile $seqmapcountsgff result
    fi

    # and a separate single mapper run
    seqmapsinglecountsgff=${seqmapsinglegff%gff}counts.gff
    if needsUpdate ${seqmapsinglecountsgff} $genesinglegff $refseq $BINDIR/tagged_gff_to_counts.pl 
    then
	$BINDIR/tagged_gff_to_counts.pl -g ${genesinglegff} -s ${refseq} -f ${ALIGN}${single}${lib} -N $libsize -c ${genetype} -o ${seqmapsinglecountsgff}
	registerFile $seqmapsinglecountsgff result
    fi

    # crunch gff annotations down into single entries for each unique splicesite; in particular useful for gbrowse loading

    # all mappers first
    mapcrunchgff=${map}.crunch.gff
    seqmapcrunchgff=${mapcrunchgff%.crunch.gff}.${seq}.crunch.gff

    if needsUpdate $seqmapcrunchgff $map $genegff $BINDIR/crunch_gff_with_counts.pl
    then
	$BINDIR/crunch_gff_with_counts.pl -g ${genegff} -s ${refseq} -f ${ALIGN}match${lib} -N $libsize -c ${genetype} -o ${seqmapcrunchgff}
	registerFile $seqmapcrunchgff temp
	crunch_update=yes
	updates=yes
    fi

    # then crunch single mapper files
    if [ "$single" = "yes" ] 
    then
	mapsinglecrunchgff=${mapsingle}.crunch.gff

	seqmapsinglecrunchgff=${mapsinglecrunchgff%.crunch.gff}.${seq}.crunch.gff
	
	if needsUpdate $seqmapsinglecrunchgff $mapsingle $genesinglegff $BINDIR/crunch_gff_with_counts.pl
	then
	    $BINDIR/crunch_gff_with_counts.pl -g ${genesinglegff} -s ${refseq} -f ${ALIGN}${single}${lib} -N $libsize -c ${genetype} -o ${seqmapsinglecrunchgff}
	    registerFile $seqmapsinglecrunchgff temp
	    crunchsingle_update=yes
	    updates=yes
	fi
    fi
    
    # uorf analysis, using splicemodel derived uorf finding code. 

    if [ "${runuorf}" = "yes" ]
    then
	# First all mappers
	seqmapuorftab=${seqmaptab%.tab}.uorf.tab
			
	if needsUpdate $seqmapuorftab $map $genegff $BINDIR/tagged_gff_to_counts_uorf.pl
	then
	    $BINDIR/tagged_gff_to_counts_uorf.pl -g ${genegff} -s ${refseq} -f ${ALIGN}match${lib} -N $libsize -c ${genetype}  -l $uorfutrcutofflen -o ${seqmapuorftab}
	    registerFile $seqmapuorftab temp
	    uorf_update=yes
	    updates=yes
	fi

	# then single mappers if requested
 	if [ "$single" = "yes" ] 
	then
	    seqmapsingleuorftab=${seqmapsingletab%.tab}.uorf.tab
	
	    if needsUpdate $seqmapsingleuorftab $mapsingle $genesinglegff $BINDIR/tagged_gff_to_counts_uorf.pl
	    then
		$BINDIR/tagged_gff_to_counts_uorf.pl -g ${genesinglegff} -s ${refseq} -f ${ALIGN}${single}${lib} -N $libsize -c ${genetype} -l $uorfutrcutofflen -o ${seqmapsingleuorftab}
		registerFile $seqmapsingleuorftab temp
		uorfsingle_update=yes
		updates=yes
	    fi
	fi
    fi
done

# check update flags for agglomerate multi sequence unit files, and update the joined files as needed

# join the split crunch GFFs

: <<'POD_FUNC'

=item UpdateGFF(basefilename,suffix,searchstring)

    Find lines including searchstring in a file named almost like basefilename, but with each sequence name in the global fastalist inserted before .suffix in the basefilename.

    E.g. C<UpdateGFF $mapcrunchgff crunch.gff ${ALIGN}match>

=cut

POD_FUNC

function UpdateGFF()
{
    local this_gff=$1
    local suffix=$2
    local grepstring=$3

    cat /dev/null > $this_gff
    for seq in `cat $fastalist` ;
    do
	local this_seqgff=${this_gff%%.${suffix}}.${seq}.$suffix
	grep $grepstring ${this_seqgff} >> ${this_gff}
    done
}

if [ "${crunch_update}" = "yes" ] 
then
    mapcrunchgff=${map}.crunch.gff
    UpdateGFF $mapcrunchgff crunch.gff ${ALIGN}match
    registerFile $mapcrunchgff result
    updates=yes
fi

if [ "${crunchsingle_update}" = "yes" ] 
then
    
    mapsinglecrunchgff=${mapsingle}.crunch.gff    
    UpdateGFF $mapsinglecrunchgff crunch.gff ${ALIGN}${single}${lib}
    registerFile $mapsinglecrunchgff result
    updates=yes
fi

# join the tabs

: <<'POD_FUNC'

=item UpdateTab(tabfilename,suffix)

Update table tabfilename, by inserting sequence names from fasta list directly followed by the contents of per-sequence-sub-tabfilename with the sequnece name from the the global fastalist inserted before .suffix in the basefilename. E.g. C<UpdateTab ${mapsingleuorftab} uorf.tab>

=cut

POD_FUNC

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
    seqmaptab=${maptab%.tab}.${seq}.tab

    if needsUpdate $maptab $seqmaptab
    then
	map_update=yes
    fi

    # check single mapper counts tab seqtabs    
    if [ "$single" = "yes" ]
    then
	seqmapsingletab=${mapsingletab%.tab}.${seq}.tab

	if needsUpdate $mapsingletab $seqmapsingletab
	then
	    mapsingle_update=yes
	fi
    fi
    
    # check uorf tab seqtabs
    if [ "${runuorf}" = "yes" ]
    then
	seqmapuorftab=${seqmaptab%.tab}.uorf.tab
	
	if needsUpdate $mapuorftab $seqmapuorftab
	then
	    uorf_update=yes
	fi

	seqmapsingleuorftab=${seqmapsingletab%.tab}.uorf.tab
	if [ "$single" = "yes" ]
	then
	    if needsUpdate $mapsingleuorftab $seqmapsingleuorftab
	    then
		uorfsingle_update=yes
	    fi
	fi
    fi
done

if [ "${map_update}" = "yes" ] 
then
    UpdateTab ${maptab} tab
    registerFile ${maptab} result
fi

if [ "${mapsingle_update}" = "yes" ] 
then
    UpdateTab ${mapsingletab} tab
    registerFile ${mapsingletab} result
fi

# join the uORF tabs

if [ "${uorf_update}" = "yes" ]
then
    UpdateTab $mapuorftab uorf.tab
    registerFile ${mapuorftab} result
fi
  
if [ "$uorfsingle_update" = "yes" ]
then
    UpdateTab ${mapsingleuorftab} uorf.tab
    registerFile ${mapsingleuorftab} result
fi

: <<'POD_FUNC'

=item ExtractSites(optdesc,optstr)

E.g. C<ExtractSites "internal.once" "-E -1">

=cut

POD_FUNC

function ExtractSites()
{
    optdesc=$1
    optstr=$2

    splicesitesout=${map}.${optdesc}.splicesites.out
    
    splicesitesout_update=no

    for seq in `cat $fastalist`
    do
	genegff=${map}.${seq}.gene.gff
	refseq=${reffile%.fasta}.${seq}.fasta
	seqsplicesitesout=${map}.${seq}.${optdesc}.splicesites.out
		
	if needsUpdate $splicesitesout $reffile $genegff $BINDIR/extract_splicesites.pl $SEQLOGOBIN
	then
	    $BINDIR/extract_splicesites.pl -g ${genegff} -s ${refseq} -f ${ALIGN}match${lib} -c ${genetype} -o $seqsplicesitesout -b $beforeag -a $afterag $optstr
	    registerFile $seqsplicesitesout temp
	    splicesitesout_update=yes
	fi
    done

    if [ "$splicesitesout_update" = "yes" ] 
    then
	echo -n > $splicesitesout
	for seq in `cat $fastalist`
	do
	    seqsplicesitesout=${map}.${seq}.${optdesc}.splicesites.out
	    cat $seqsplicesitesout >> $splicesitesout
	done    
	registerFile $splicesitesout temp

	startpos=$(( $beforeag + 2 ))
	optdescnodot=${optdesc//./ }
	
	$SEQLOGOBIN -f $splicesitesout -h 8 -w 26 -n -Y -c -k1 -s -${startpos} -d 1 -t "Spliced leader Addition Site logo ${optdescnodot} ${lib}" -F PDF -o $splicesitesout
	registerFile ${splicesitesout}.pdf result
    fi
    
    splicesitesstatsummary=$splicesitesout.stat.summary
    if needsUpdate $splicesitesstatsummary $splicesitesout $BINDIR/splicesite_stats.pl
    then
	$BINDIR/splicesite_stats.pl < $splicesitesout > $splicesitesstatsummary
#	perl -e '%accdinuc; my $sum=0; print "Dinuc\tNum\tFraction\n"; while ($str=<STDIN>) { chomp $str; my $dinuc=substr $str,-8,2; $accdinuc{$dinuc}++; $sum++} for my $dinuc (sort { $accdinuc{$b}<=>$accdinuc{$a}} keys %accdinuc) { print $dinuc,"\t", $accdinuc{$dinuc}, "\t", sprintf("%.2f",$accdinuc{$dinuc}/$sum), "\n"; }' 
	registerFile $splicesitesstatsummary result
    fi

	# max p
	# reverse

#	while($thisseq =~ /([tcyn]+)([agrnx]?[tcynx]*){0,$max_poly_py_mismatches}([tcynx]+)/g) {
#	    print 

}

if [ "$plotsplicesites" == "yes" ] 
then
    ExtractSites "once" "-1"
    ExtractSites "major.once" "-m -1"
    ExtractSites "minor.once" "-M -1"
    ExtractSites "internal.once" "-E -1"
    ExtractSites "external.once" "-I -1"
fi

: <<'POD_FUNC'

=item CreateSubTable(basetab, tabname, awkscript[, comment_for_summary])

Function to make filtered sub-tables from the big ones

createsubtable basetab_filename subtabname converting_awk_script [comment_for_summary]

E.g. C<CreateSubTable $mapsingletab "genes_with_major_UTR_gt2k" '($9 > 2000) {print}' "${single}-mappers">

=cut

POD_FUNC

# function to make filtered sub-tables from the big ones
function CreateSubTable()
{
    # createsubtable basetab_filename subtabname converting_awk_script [comment_for_summary]
    local basetab=$1
    local tabname=$2
    local awkscript=$3

    local comment=""
    if [ $# -gt 3 ]
    then
	comment=" $4"
    fi

    subtab=${basetab%.tab}.${tabname}.tab
    if needsUpdate $subtab $maptab
    then
	echo -n > $subtab
	for seq in `cat $fastalist`
	do
	    seqbasetab=${basetab%.tab}.${seq}.tab
#	    echo try awk \"$awkscript\" $seqbasetab to $subtab
	    awk "$awkscript" $seqbasetab >> $subtab
	done
	registerFile $subtab temp
	updates=yes
    fi
     
    echo -n "${lib} ${tabname}${comment}" >> $summary
    wc -l $subtab |awk '{ print "\t",$1 }' >> $summary
}

: <<'POD_FUNC'

=item OneHeader(basetab_filename subtabname header_starts_with_word)

Get rid of all but one header line after joining.

E.g. C<CreateSubTable $maptab "genes_with_major_UTR_lt2k" 'BEGIN ...$10}'> and then 
            C<OneHeader $maptab "genes_with_major_UTR_lt2k" "Name">

=cut

POD_FUNC

# get rid of all but one header line after joining..
function OneHeader()
{
    # oneheader basetab_filename subtabname header_starts_with_word
    local basetab=$1 #map/mapsingle
    local tabname=$2
    local headerstart=$3

    subtab=${basetab%.tab}.${tabname}.tab
    subtabtmp=${subtab}.tmp
    grep ^${headerstart} $subtab |head -1> $subtabtmp 
    grep -v ^${headerstart} $subtab >> $subtabtmp
    mv $subtabtmp $subtab
}

: <<'POD_FUNC'

=item RUTRlenstats(tabname liblabel)

Write custom R script for calculating mean and median UTR lengths and plotting their distributions and run it.
Generates pdf plots.

E.g. C<RUTRlenstats "genes_with_major_UTR_lt2k" "major splice site below 2 kbp">

=cut

POD_FUNC

# R script for calculating mean and median UTR lengths and plotting their distributions
function RUTRlenstats()
{
    tabname=$1
    liblabel=$2
    subtab=${map}.${tabname}.tab

    cat > ${map}.utrlen.R <<-RUTRLEN
	lib<-read.table("${subtab}", header = TRUE, sep = "\\t");

	ctcountutrlen <- cor.test(lib\$Count, lib\$fpUTRlen)
	corrstr<-sprintf("The estimated correlation between 5pUTRlen and count, $lib $liblabel: %.3f, p = %.3e", ctcountutrlen\$estimate, ctcountutrlen\$p.value)
	cat(corrstr,file="${summary}",sep="\\n", append=TRUE)

	libutrlen <- rep(lib\$fpUTRlen,lib\$Count)

	meanutrlen <- mean(libutrlen)
	meanstr <- sprintf("Mean UTR len, $lib $liblabel: %.2f",meanutrlen)
	cat(meanstr,file="${summary}", sep="\\n", append=TRUE)

	medianutrlen <- median(libutrlen)
	medianstr <- sprintf("Median UTR len, $lib $liblabel: %.1f",medianutrlen)
	cat(medianstr,file="${summary}", sep="\\n", append=TRUE)

	pdf("${subtab%.tab}.UTRlen.pdf")
	hist(libutrlen,100,main="${lib} ${liblabel} 5'UTR length, transcript counts")
	dev.off()

	pdf("${subtab%.tab}.UTRlen.genes.pdf")
	hist(lib\$fpUTRlen,100,main="${lib} ${liblabel} 5'UTR length, gene counts")
	dev.off()

	RUTRLEN

    R CMD BATCH --no-save --slave ${map}.utrlen.R /dev/null

    registerFile ${map}.utrlen.R temp
    registerFile ${subtab%.tab}.UTRlen.genes.pdf result
    registerFile ${subtab%.tab}.UTRlen.pdf result
}

summary=${maptab%.tab}.summary

if needsUpdate $summary $map.txt ${mapcrunchgff} $maptab ${mapuorftab} ${mapsingletab} ${mapsingleuorftab} ${mapsingle}.txt $BINDIR/sites_per_gene_summary.pl $BINDIR/median.pl
then
    echo -n > $summary

    # produce a nice summary stats file..
    echo "${lib} Number of reads: $libsize" >> $summary

    echo -n "${lib} Number of reads aligned ok: " >> $summary
    wc -l $map.txt |awk '{ print "\t",$1 }' >> $summary

    # percent aligned ok.. 

    echo -n "${lib} number of $ALIGN singlemappers, (maq: q >= ${alqt}, bowtie: one reportable): " >> $summary
    wc -l ${mapsingle}.txt |awk '{ print "\t",$1 }' >> $summary

    echo -n "${lib} Number of $ALIGN multimappers, (maq: q == 0, bowtie more reportable): " >> $summary
    wc -l ${map}.multi.txt |awk '{ print "\t",$1 }' >> $summary

    # use a perl script for the median calculations 

    echo -n "${lib} ${ALIGN} single mappers alignment quality: " >> $summary
    awk '($7>0) {print $7;}' ${map}.txt | $BINDIR/median.pl >> $summary

    echo -n "${lib} ${ALIGN} single mappers mismatches: " >> $summary
    awk '($7>0) {print $10;}' ${map}.txt | $BINDIR/median.pl >> $summary
    
    echo -n "${lib} $ALIGN all mapped reads mismatches: " >> $summary
    awk '{print $10;}' ${map}.txt | $BINDIR/median.pl >> $summary
    
    echo -n "${lib} " >> $summary; 
    echo -n "Number of splicesites " >> $summary
    grep -c ${ALIGN}match ${mapcrunchgff} >> $summary

    echo -n "${lib} ">> $summary
    echo -n "Genes with counts">> $summary
    awk 'BEGIN {tick = 0} ($4 > 0) {tick=tick+1} END {print " ",tick}' $maptab >> $summary;

    #echo -n "${lib} Genes_with_lt2k_UTR" >> $summary
    #awk 'BEGIN {tick = 0} ($4>0 && $11 < 2000) {tick=tick+1} END {print " ",tick}' $maptab >> $summary

    echo -n "${lib} Genes_with_gt2k_lt5k_UTR" >> $summary
    awk 'BEGIN {tick = 0} ($4>0 && $11 >= 2000 && $11<5000) {tick=tick+1} END {print " ",tick}' $maptab >> $summary

    CreateSubTable $maptab "genes_with_major_site_internal" '($4 > 0 && $9 < 0) {print}'
    CreateSubTable $maptab "genes_with_only_internal_sites" '($4 > 0 && $12 < 0) {print}'

    CreateSubTable $maptab "genes_with_shortest_UTR_gt2k" '($11 > 2000) {print}'
    CreateSubTable $maptab "genes_with_major_UTR_gt2k" '($9 > 2000) {print}'
    
    CreateSubTable $maptab "genes_with_major_UTR_lt2k" 'BEGIN {print "Name\tPos\tfpUTRlen\tCount"} ($4 > 0 && $9 <= 2000 && $9 >= 0 ) {print $1"\t"$2"\t"$9"\t"$10}'
        
    OneHeader $maptab "genes_with_major_UTR_lt2k" "Name"

    if [ "$single" = "yes" ] 
    then
        CreateSubTable $mapsingletab "genes_with_major_site_internal" '($4 > 0 && $9 < 0) {print}' "${single}-mappers"
	CreateSubTable $mapsingletab "genes_with_only_internal_sites" '($4 > 0 && $12 < 0) {print}' "${single}-mappers"

	CreateSubTable $mapsingletab "genes_with_shortest_UTR_gt2k" '($11 > 2000) {print}' "${single}-mappers"
	CreateSubTable $mapsingletab "genes_with_major_UTR_gt2k" '($9 > 2000) {print}' "${single}-mappers"
    
	CreateSubTable $mapsingletab "genes_with_major_UTR_lt2k" 'BEGIN {print "Name\tPos\tfpUTRlen\tCount"} ($4 > 0 && $9 <= 2000 && $9 >= 0 ) {print $1"\t"$2"\t"$9"\t"$10}' "${single}-mappers"
        
	OneHeader $mapsingletab "genes_with_major_UTR_lt2k" "Name"
    fi  
    
    # count uORFs - no need for a separate tab, yet?
    if [ "$runuorf" = "yes" ]
    then
	echo -n "$lib genes with at least one uORF: " >> $summary
	awk 'BEGIN { sum=0 } ($14 >0) { sum=sum+1; } END { print sum }' $mapuorftab >> $summary

	if [ "$single" = "yes" ] 
	then
	    echo -n "$lib ${single} genes with at least one uORF: " >> $summary
	    awk 'BEGIN { sum=0 } ($14 >0) { sum=sum+1; } END { print sum }' $mapsingleuorftab >> $summary
	fi
    fi


#    echo -n "$lib genes with upstream starts :" >> $summary
#    awk 'BEGIN { sum=0 } ($14 >0) { sum=sum+1; } END { print sum }'

    # process tab to produce one with all UTRlens with counts expanded on separate lines (or get from crunch gff?) 
    # note RAW counts on the individual UTRs! correct for line separated tab seems to be the most straightforward

    allposutrs=${map}.all_pos_UTRs_lt2k.tab
    perl -e 'print "Name\tPos\tfpUTRlen\tCount\n"; while(<STDIN>) { chomp; @r=split(/\t+/); if ($r[3] > 0) { @counts=split(/,/,$r[5]); @UTRlen=split(/,/,$r[7]); for ($i=0; $i<@counts;$i++) { if($UTRlen[$i]>0 && $UTRlen[$i]<2000) { $normcount=$counts[$i]/'$libsize'*1000000; print "$r[0]\t$r[1]\t$UTRlen[$i]\t".sprintf("%.0f",$normcount)."\n"; } } } }' < $maptab > $allposutrs

    registerFile $allposutrs result

    RUTRlenstats "genes_with_major_UTR_lt2k" "major splice site below 2 kbp"
    RUTRlenstats "all_pos_UTRs_lt2k" "all UTRS between 0 nt and 2 knt"

    # perl script for the median; a heredoc might be better, but its a little cumbersome to edit 

    echo "$lib Number of sites/gene counting UTRs between 0 and 2000 nt " >> $summary
    $BINDIR/sites_per_gene_summary.pl < $maptab >> $summary

    # expression level summary: tags/gene measures, plot and a little tab with the 30 most highly expressed ones
    # including count 0 genes for some kind of completeness. some tuning can be achieved using the gene feature type.
 
    topcounts=${maptab%tab}topcounts

    awk 'BEGIN { print "Name","Count","Desc"; } ($5 >= 0) { print $1,$5,$13; }' $maptab > $topcounts

    cat > ${map}.tags.R <<-RTAGS
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
    
    R CMD BATCH --no-save --slave ${map}.tags.R /dev/null
    registerFile ${map}.tags.R temp
    registerFile ${topcounts} temp
    
    showtoptagged=40
    echo >> $summary
    echo "Genes with the $showtoptagged highest tag counts" >> $summary    
    head -1 $topcounts >> $summary
    grep -v ^Name $topcounts |sort -k2,2nr | head -$showtoptagged >> $summary

    registerFile $summary result
    updates=yes
fi

if [ "$updates" = "no" ]
then
    log=${map}.log
    rundate=`date`
    echo "Project was already up to date. Nothing to do. Please update date on appropriate input/program/intermediate results files to force results updates. ($rundate)" | tee -a $log
else
    log=${map}.log
    rundate=`date`
    echo "Project has been brought up to date. ($rundate)" | tee -a $log
fi
exit

: <<'POD_EOF'

=back

=cut

POD_EOF

# quality checks

# number of unaligned?

# number of rather ok single-mappers?

# average align quality?

# join size split reads over a certain size...
# average (trimmed) read size
#grep -c ^@ GDG-*/*GDG*.fastq > libsizecounts 
#cat libsizecounts |sed -e 's/_/:/g' |cut -d: -f 2,8 |sed -e 's/:/ /' |sort -k1,1n |grep -v Adapter |cut -d" " -f2 |column -t -c4 |perl -e '$count=0;$partsum=0; while(<STDIN>) { chomp; $count++; $partsum+=$_; if($count % 4 == 0) { print $count/4-1, "\t", $partsum, "\n"; $partsum = 0; } }' > libsizecounts.sum
#R
#<<-"RSCRIPT"
#libc<-read.table("libsizecounts.sum");
#libcraw<-rep(libc$V1, libc$V2)
#mean(libcraw)
#median(libcraw)
#RSCRIPT

# randomised trials -- belongs in a wrapper script, not here

#for lib in 1 2 3 4 ; do for count in 5000 25000 50000 100000 150000 200000 250000 300000 350000 400000 450000 500000 ; do echo GDG$lib - $count; ../install/bin/randomly_pick_fastq.pl -n $count -N `grep -c ^@ ../fasteris/GDG-${lib}/Size24andup.fastq` -o GDG-${lib}.Size24andup.rnd${count}.fastq /Users/daniel/fasteris/GDG-$lib/Size24andup.fastq ; done ; done

