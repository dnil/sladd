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

# MAQBIN is the MAQ installation directory (path to the maq binary)
if [ -z "$MAQBIN" ]
then
	#MAQBIN=~/bin/maq
	MAQBIN=~/install/maq-0.7.1/maq
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

. pipelinefunk.sh

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

refbfa=${reffile%fasta}bfa

if needsUpdate $refbfa $reffile
then
    $MAQBIN fasta2bfa $reffile $refbfa
    updates=yes
fi

# count the number of reads in tagfile
# [TODO] save result for repetitive script use?

libsize=`grep -c ^@ $tagfile`

# maq performs optimally around 2M reads, according to its manual pages
# if the number of reads is considerably larger than that (3M), split before bfq generation, 
# map separately and mapmerge the alignments.

maqmap=${tagfile%.fastq}_vs_${reffile%.fasta}.maqmap

if [ $libsize -gt 3000000 ] 
then
    tagfileprefix=${tagfile}.split
    maqmapupdate="no"
    
    # split into 8M/4=2M line files
    registerFile $tagfile.done.split temp
    if needsUpdate $tagfile.done.split $tagfile
    then
	split -l 8000000 $tagfile ${tagfileprefix}
	echo ${tagfileprefix}[a-z][a-z] > $tagfile.done.split
    fi

    splitnr=0

    # check if bfq does not exist, or is older than fastq
    for splittagfile in ${tagfileprefix}[a-z][a-z]
    do
	bfqfile=${splittagfile}.bfq
	if needsUpdate $bfqfile $splittagfile
	then	    
	    $MAQBIN fastq2bfq $splittagfile $bfqfile
	    updates=yes
	fi

	splitmaqmap=${splittagfile}_vs_${reffile%.fasta}.maqmap

	if needsUpdate $splitmaqmap $refbfa $bfqfile
	then
	    log=${maqmap}.log
	    rundate=`date`
	    echo "map -n 3 $splitmaqmap $refbfa $bfqfile : (${rundate})" >> ${log}
	    $MAQBIN map -n 3 $splitmaqmap $refbfa $bfqfile 2>> ${log}
	    updates=yes
	    maqmapupdate="yes"
	fi
	
	splittagfiles[$splitnr]=$splittagfile
	splitmaqmap[$splitnr]=$splitmaqmap
	splitnr=$(( $splitnr + 1 ))
    done

    if [ "$maqmapupdate" = "yes" ] 
    then
	$MAQBIN mapmerge $maqmap ${splitmaqmap[@]}
    fi
else

    # check if bfq does not exist, or is older than fastq

    bfqfile=${tagfile%fastq}bfq

    if needsUpdate $bfqfile $tagfile
    then
	$MAQBIN fastq2bfq $tagfile $bfqfile
	updates=yes
    fi

    # map
    # size-sep to avoid length cutoff?
    if needsUpdate $maqmap $refbfa $bfqfile
    then
	log=${maqmap}.log
	rundate=`date`
	echo "map -n 3 $maqmap $refbfa $bfqfile : (${rundate})" >>${log}
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
    fi

    # update ref fasta split if there is a new reffile available

    refseq=${reffile%.fasta}.${seq}.fasta 

    if needsUpdate ${refseq} ${reffile} $BINDIR/fasta_header_grep.pl
    then
	$BINDIR/fasta_header_grep.pl "${grepsafeseq}" ${reffile} > ${refseq}
    fi

done

# produce main gff files

maqmaptab=${maqmap}.counts.tab
maqmapsingletab=${maqmapsingle}.counts.tab

maqmapuorftab=${maqmaptab%.tab}.uorf.tab
maqmapsingleuorftab=${maqmapsingletab%.tab}.uorf.tab

# update flags for agglomeration of multi sequence unit files 
maqmap_update=no
maqmapsingle_update=no
crunch_update=no
crunchsingle_update=no
uorf_update=no
uorfsingle_update=no

for seq in `cat $fastalist` ;
do
    # split the maqmap data (gff format) into per sequence unit files, full set of mappers

    seqmaqmapgff=${maqmapgff%.gff}.${seq}.gff
    
#    grepsafeseq=`echo ${seq}|sed -e 's/|/\\\|/g'`

    if needsUpdate $seqmaqmapgff $maqmapgff
    then
	grep ${seq} ${maqmapgff} > ${seqmaqmapgff}
    fi

    # split the maqmap data (gff format) into per sequence unit files, single mappers

    seqmaqmapsinglegff=${maqmapsinglegff%.gff}.${seq}.gff

    if [ "$single" = "yes" ]
    then
	if needsUpdate $seqmaqmapsinglegff $maqmapsinglegff
	then
	    grep ${seq} ${maqmapsinglegff} > ${seqmaqmapsinglegff}
	fi
    fi

    # the initial mapped read to gene connection in the form of gene.gffs, all mappers

    refseqgff=${reffile%.fasta}.${seq}.gff
    refseq=${reffile%.fasta}.${seq}.fasta
    
    genegff=${maqmapgff%.gff}.${seq}.gene.gff

    if needsUpdate $genegff $refseqgff $refseq $seqmaqmapgff "$BINDIR/per_gene_splicesite_info.pl"
    then
	log=${maqmap}.log
	rundate=`date`
	echo "$seq : per_gene_splicesite_info ($rundate)" >> $log
	$BINDIR/per_gene_splicesite_info.pl -g ${refseqgff} -s ${refseq} -t ${seqmaqmapgff} -f maqmatch${lib} -o ${genegff} -c ${genetype} -U ${utrcutoff} >> $log
    fi

    # initial gene-gffs, single mappers

    genesinglegff=${maqmapsinglegff%.gff}.${seq}.gene.gff
    if [ "$single" = "yes" ]
    then
	if needsUpdate $genesinglegff $refseqgff $refseq $seqmaqmapsinglegff "$BINDIR/per_gene_splicesite_info.pl"
	then
	    log=${maqmap}.log
	    rundate=`date`
	    echo "$seq : per_gene_splicesite_info ($rundate) : single mappers : " >> $log
	    $BINDIR/per_gene_splicesite_info.pl -g ${refseqgff} -s ${refseq} -t ${seqmaqmapsinglegff} -f maqq${alqt}${lib} -o ${genesinglegff} -c ${genetype} -U ${utrcutoff} >> $log
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
    fi

    # antisense tags, single mappers

    astagsinglegff=${maqmapsinglegff%.gff}.${seq}.astags.gff
    if needsUpdate $astagsinglegff $refseqgff $refseq $seqmaqmapsinglegff "$BINDIR/list_antisense_splicesites.pl"
    then
	aslogsingle=${maqmapsingle}.antisense.log
	rundate=`date`
	echo "$seq : list_antisense_splicesites ($rundate)" >> $aslogsingle
	$BINDIR/list_antisense_splicesites.pl -g ${refseqgff} -s ${refseq} -t ${seqmaqmapsinglegff} -f maqq${alqt}${lib} -o ${astagsinglegff} -c ${genetype} >> $aslogsingle
    fi

    astagsinglecrunchgff=${maqmapsinglegff%.gff}.${seq}.astags.crunch.gff
    if needsUpdate $astagsinglecrunchgff $maqmapsingle $astagsinglegff $BINDIR/crunch_gff_with_counts.pl
    then
	$BINDIR/crunch_gff_with_counts.pl -g ${astagsinglegff} -s ${refseq} -f maqq${alqt}${lib} -N $libsize -c ${genetype} -o ${astagsinglecrunchgff}
#	crunchsingle_update=yes
	updates=yes
    fi

    # counts tabs, reporting most of the per gene mapping data in a tabular format, for all mapping reads

    seqmaqmaptab=${maqmaptab%.tab}.${seq}.tab

    if needsUpdate $seqmaqmaptab $genegff $refseq $BINDIR/tagged_gff_to_counts_tab.pl
    then
	$BINDIR/tagged_gff_to_counts_tab.pl -g ${genegff} -s ${refseq} -f maqmatch${lib} -N $libsize -c ${genetype} -o ${seqmaqmaptab}
	maqmap_update=yes
	updates=yes
    fi
   
    # counts tabs, single mappers

    seqmaqmapsingletab=${maqmapsingletab%.tab}.${seq}.tab
    if [ "$single" = "yes" ] 
    then
	if needsUpdate $seqmaqmapsingletab $genesinglegff $refseq $BINDIR/tagged_gff_to_counts_tab.pl
	then
	    $BINDIR/tagged_gff_to_counts_tab.pl -g ${genesinglegff} -s ${refseq} -f maqq${alqt}${lib} -N $libsize -c ${genetype} -o ${seqmaqmapsingletab}
	    maqmapsingle_update=yes
	    updates=yes
	fi
    fi

    # [TODO] not sure that we're actually using this now? useful for loading completed gene features onto gbrowse? fix multi-utr-len entry?!
    # all mappers
    seqmaqmapcountsgff=${seqmaqmapgff%gff}counts.gff
    if needsUpdate ${seqmaqmapcountsgff} $genegff $refseq $BINDIR/tagged_gff_to_counts.pl 
    then
	$BINDIR/tagged_gff_to_counts.pl -g ${genegff} -s ${refseq} -f maqmatch${lib} -N $libsize -c ${genetype} -o ${seqmaqmapcountsgff}
    fi

    # and a separate single mapper run
    seqmaqmapsinglecountsgff=${seqmaqmapsinglegff%gff}counts.gff
    if needsUpdate ${seqmaqmapsinglecountsgff} $genesinglegff $refseq $BINDIR/tagged_gff_to_counts.pl 
    then
	$BINDIR/tagged_gff_to_counts.pl -g ${genesinglegff} -s ${refseq} -f maqq${alqt}${lib} -N $libsize -c ${genetype} -o ${seqmaqmapsinglecountsgff}
    fi    

    # crunch gff annotations down into single entries for each unique splicesite; in particular useful for gbrowse loading

    # all mappers first
    maqmapcrunchgff=${maqmap}.crunch.gff
    seqmaqmapcrunchgff=${maqmapcrunchgff%.crunch.gff}.${seq}.crunch.gff

    if needsUpdate $seqmaqmapcrunchgff $maqmap $genegff $BINDIR/crunch_gff_with_counts.pl
    then
	$BINDIR/crunch_gff_with_counts.pl -g ${genegff} -s ${refseq} -f maqmatch${lib} -N $libsize -c ${genetype} -o ${seqmaqmapcrunchgff}
	crunch_update=yes
	updates=yes
    fi

    # then crunch single mapper files
    if [ "$single" = "yes" ] 
    then
	maqmapsinglecrunchgff=${maqmapsingle}.crunch.gff

	seqmaqmapsinglecrunchgff=${maqmapsinglecrunchgff%.crunch.gff}.${seq}.crunch.gff
	
	if needsUpdate $seqmaqmapsinglecrunchgff $maqmapsingle $genesinglegff $BINDIR/crunch_gff_with_counts.pl
	then
	    $BINDIR/crunch_gff_with_counts.pl -g ${genesinglegff} -s ${refseq} -f maqq${alqt}${lib} -N $libsize -c ${genetype} -o ${seqmaqmapsinglecrunchgff}
	    crunchsingle_update=yes
	    updates=yes
	fi
    fi
    
    # uorf analysis, using splicemodel derived uorf finding code. 

    if [ "${runuorf}" = "yes" ]
    then
	# First all mappers
	seqmaqmapuorftab=${seqmaqmaptab%.tab}.uorf.tab
			
	if needsUpdate $seqmaqmapuorftab $maqmap $genegff $BINDIR/tagged_gff_to_counts_uorf.pl
	then
	    $BINDIR/tagged_gff_to_counts_uorf.pl -g ${genegff} -s ${refseq} -f maqmatch${lib} -N $libsize -c ${genetype}  -l $uorfutrcutofflen -o ${seqmaqmapuorftab}
	    uorf_update=yes
	    updates=yes
	fi

	# then single mappers if requested
 	if [ "$single" = "yes" ] 
	then
	    seqmaqmapsingleuorftab=${seqmaqmapsingletab%.tab}.uorf.tab
	
	    if needsUpdate $seqmaqmapsingleuorftab $maqmapsingle $genesinglegff $BINDIR/tagged_gff_to_counts_uorf.pl
	    then
		$BINDIR/tagged_gff_to_counts_uorf.pl -g ${genesinglegff} -s ${refseq} -f maqq${alqt}${lib} -N $libsize -c ${genetype} -l $uorfutrcutofflen -o ${seqmaqmapsingleuorftab}
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

    E.g. C<UpdateGFF $maqmapcrunchgff crunch.gff maqmatch>

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
    maqmapcrunchgff=${maqmap}.crunch.gff
    UpdateGFF $maqmapcrunchgff crunch.gff maqmatch
    updates=yes
fi

if [ "${crunchsingle_update}" = "yes" ] 
then
    
    maqmapsinglecrunchgff=${maqmapsingle}.crunch.gff    
    UpdateGFF $maqmapsinglecrunchgff crunch.gff maqq${alqt}${lib}

    updates=yes
fi

# join the tabs

: <<'POD_FUNC'

=item UpdateTab(tabfilename,suffix)

Update table tabfilename, by inserting sequence names from fasta list directly followed by the contents of per-sequence-sub-tabfilename with the sequnece name from the the global fastalist inserted before .suffix in the basefilename. E.g. C<UpdateTab ${maqmapsingleuorftab} uorf.tab>

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
    
    # check uorf tab seqtabs
    if [ "${runuorf}" = "yes" ]
    then
	seqmaqmapuorftab=${seqmaqmaptab%.tab}.uorf.tab
	
	if needsUpdate $maqmapuorftab $seqmaqmapuorftab
	then
	    uorf_update=yes
	fi

	seqmaqmapsingleuorftab=${seqmaqmapsingletab%.tab}.uorf.tab
	if [ "$single" = "yes" ]
	then
	    if needsUpdate $maqmapsingleuorftab $seqmaqmapsingleuorftab
	    then
		uorfsingle_update=yes
	    fi
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

# join the uORF tabs

if [ "${uorf_update}" = "yes" ]
then
    UpdateTab $maqmapuorftab uorf.tab
fi
  
if [ "$uorfsingle_update" = "yes" ]
then
    UpdateTab ${maqmapsingleuorftab} uorf.tab
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

    splicesitesout=${maqmap}.${optdesc}.splicesites.out
    
    splicesitesout_update=no

    for seq in `cat $fastalist`
    do
	genegff=${maqmap}.${seq}.gene.gff
	refseq=${reffile%.fasta}.${seq}.fasta
	seqsplicesitesout=${maqmap}.${seq}.${optdesc}.splicesites.out
		
	if needsUpdate $splicesitesout $reffile $genegff $BINDIR/extract_splicesites.pl $SEQLOGOBIN
	then
	    $BINDIR/extract_splicesites.pl -g ${genegff} -s ${refseq} -f maqmatch${lib} -c ${genetype} -o $seqsplicesitesout -b $beforeag -a $afterag $optstr
	    splicesitesout_update=yes
	fi
    done

    if [ "$splicesitesout_update" = "yes" ] 
    then
	echo -n > $splicesitesout
	for seq in `cat $fastalist`
	do
	    seqsplicesitesout=${maqmap}.${seq}.${optdesc}.splicesites.out
	    cat $seqsplicesitesout >> $splicesitesout
	done    

	startpos=$(( $beforeag + 2 ))
	optdescnodot=${optdesc//./ }
	
	$SEQLOGOBIN -f $splicesitesout -h 8 -w 26 -n -Y -c -k1 -s -${startpos} -d 1 -t "Spliced leader Addition Site logo ${optdescnodot} ${lib}" -F PDF -o $splicesitesout
    fi
    
    splicesitesstatsummary=$splicesitesout.stat.summary
    if needsUpdate $splicesitesstatsummary $splicesitesout $BINDIR/splicesite_stats.pl
    then
	$BINDIR/splicesite_stats.pl < $splicesitesout > $splicesitesstatsummary
#	perl -e '%accdinuc; my $sum=0; print "Dinuc\tNum\tFraction\n"; while ($str=<STDIN>) { chomp $str; my $dinuc=substr $str,-8,2; $accdinuc{$dinuc}++; $sum++} for my $dinuc (sort { $accdinuc{$b}<=>$accdinuc{$a}} keys %accdinuc) { print $dinuc,"\t", $accdinuc{$dinuc}, "\t", sprintf("%.2f",$accdinuc{$dinuc}/$sum), "\n"; }' 
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

E.g. C<CreateSubTable $maqmapsingletab "genes_with_major_UTR_gt2k" '($9 > 2000) {print}' "q${alqt}-mappers">

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
    if needsUpdate $subtab $maqmaptab
    then
	echo -n > $subtab
	for seq in `cat $fastalist`
	do
	    seqbasetab=${basetab%.tab}.${seq}.tab
#	    echo try awk \"$awkscript\" $seqbasetab to $subtab
	    awk "$awkscript" $seqbasetab >> $subtab
	done
	
	updates=yes
    fi
     
    echo -n "${lib} ${tabname}${comment}" >> $summary
    wc -l $subtab |awk '{ print "\t",$1 }' >> $summary
}

: <<'POD_FUNC'

=item OneHeader(basetab_filename subtabname header_starts_with_word)

Get rid of all but one header line after joining.

E.g. C<CreateSubTable $maqmaptab "genes_with_major_UTR_lt2k" 'BEGIN ...$10}'> and then 
            C<OneHeader $maqmaptab "genes_with_major_UTR_lt2k" "Name">

=cut

POD_FUNC

# get rid of all but one header line after joining..
function OneHeader()
{
    # oneheader basetab_filename subtabname header_starts_with_word
    local basetab=$1 #maqmap/maqmapsingle
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
    subtab=${maqmap}.${tabname}.tab

    cat > ${maqmap}.utrlen.R <<-RUTRLEN
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

    R CMD BATCH --no-save --slave ${maqmap}.utrlen.R /dev/null
}

summary=${maqmaptab%.tab}.summary

if needsUpdate $summary $maqmap.txt ${maqmapcrunchgff} $maqmaptab ${maqmapuorftab} ${maqmapsingletab} ${maqmapsingleuorftab} ${maqmapsingle}.txt $BINDIR/sites_per_gene_summary.pl $BINDIR/median.pl
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
    
    echo -n "${lib} " >> $summary; 
    echo -n "Number of splicesites " >> $summary
    grep -c maqmatch ${maqmapcrunchgff} >> $summary

    echo -n "${lib} ">> $summary
    echo -n "Genes with counts">> $summary
    awk 'BEGIN {tick = 0} ($4 > 0) {tick=tick+1} END {print " ",tick}' $maqmaptab >> $summary;

    #echo -n "${lib} Genes_with_lt2k_UTR" >> $summary
    #awk 'BEGIN {tick = 0} ($4>0 && $11 < 2000) {tick=tick+1} END {print " ",tick}' $maqmaptab >> $summary

    echo -n "${lib} Genes_with_gt2k_lt5k_UTR" >> $summary
    awk 'BEGIN {tick = 0} ($4>0 && $11 >= 2000 && $11<5000) {tick=tick+1} END {print " ",tick}' $maqmaptab >> $summary

    CreateSubTable $maqmaptab "genes_with_major_site_internal" '($4 > 0 && $9 < 0) {print}'
    CreateSubTable $maqmaptab "genes_with_only_internal_sites" '($4 > 0 && $12 < 0) {print}'

    CreateSubTable $maqmaptab "genes_with_shortest_UTR_gt2k" '($11 > 2000) {print}'
    CreateSubTable $maqmaptab "genes_with_major_UTR_gt2k" '($9 > 2000) {print}'
    
    CreateSubTable $maqmaptab "genes_with_major_UTR_lt2k" 'BEGIN {print "Name\tPos\tfpUTRlen\tCount"} ($4 > 0 && $9 <= 2000 && $9 >= 0 ) {print $1"\t"$2"\t"$9"\t"$10}'
        
    OneHeader $maqmaptab "genes_with_major_UTR_lt2k" "Name"

    if [ "$single" = "yes" ] 
    then
        CreateSubTable $maqmapsingletab "genes_with_major_site_internal" '($4 > 0 && $9 < 0) {print}' "q${alqt}-mappers"
	CreateSubTable $maqmapsingletab "genes_with_only_internal_sites" '($4 > 0 && $12 < 0) {print}' "q${alqt}-mappers"

	CreateSubTable $maqmapsingletab "genes_with_shortest_UTR_gt2k" '($11 > 2000) {print}' "q${alqt}-mappers"
	CreateSubTable $maqmapsingletab "genes_with_major_UTR_gt2k" '($9 > 2000) {print}' "q${alqt}-mappers"
    
	CreateSubTable $maqmapsingletab "genes_with_major_UTR_lt2k" 'BEGIN {print "Name\tPos\tfpUTRlen\tCount"} ($4 > 0 && $9 <= 2000 && $9 >= 0 ) {print $1"\t"$2"\t"$9"\t"$10}' "q${alqt}-mappers"
        
	OneHeader $maqmapsingletab "genes_with_major_UTR_lt2k" "Name"
    fi  
    
    # count uORFs - no need for a separate tab, yet?
    if [ "$runuorf" = "yes" ]
    then
	echo -n "$lib genes with at least one uORF: " >> $summary
	awk 'BEGIN { sum=0 } ($14 >0) { sum=sum+1; } END { print sum }' $maqmapuorftab >> $summary

	if [ "$single" = "yes" ] 
	then
	    echo -n "$lib q${alqt} genes with at least one uORF: " >> $summary
	    awk 'BEGIN { sum=0 } ($14 >0) { sum=sum+1; } END { print sum }' $maqmapsingleuorftab >> $summary
	fi
    fi


#    echo -n "$lib genes with upstream starts :" >> $summary
#    awk 'BEGIN { sum=0 } ($14 >0) { sum=sum+1; } END { print sum }'

    # process tab to produce one with all UTRlens with counts expanded on separate lines (or get from crunch gff?) 
    # note RAW counts on the individual UTRs! correct for line separated tab seems to be the most straightforward

    allposutrs=${maqmap}.all_pos_UTRs_lt2k.tab
    perl -e 'print "Name\tPos\tfpUTRlen\tCount\n"; while(<STDIN>) { chomp; @r=split(/\t+/); if ($r[3] > 0) { @counts=split(/,/,$r[5]); @UTRlen=split(/,/,$r[7]); for ($i=0; $i<@counts;$i++) { if($UTRlen[$i]>0 && $UTRlen[$i]<2000) { $normcount=$counts[$i]/'$libsize'*1000000; print "$r[0]\t$r[1]\t$UTRlen[$i]\t".sprintf("%.0f",$normcount)."\n"; } } } }' < $maqmaptab > $allposutrs

    RUTRlenstats "genes_with_major_UTR_lt2k" "major splice site below 2 kbp"
    RUTRlenstats "all_pos_UTRs_lt2k" "all UTRS between 0 nt and 2 knt"

    # perl script for the median; a heredoc might be better, but its a little cumbersome to edit 

    echo "$lib Number of sites/gene counting UTRs between 0 and 2000 nt " >> $summary
    $BINDIR/sites_per_gene_summary.pl < $maqmaptab >> $summary

    # expression level summary: tags/gene measures, plot and a little tab with the 30 most highly expressed ones
    # including count 0 genes for some kind of completeness. some tuning can be achieved using the gene feature type.
 
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

