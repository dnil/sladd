#!/bin/bash
#
# Daniel Nilsson, 2009
# daniel.nilsson@izb.unibe.ch, daniel.k.nilsson@gmail.com
#
# POD documentation to follow throughout the file - use e.g. perldoc to read
#

: <<'POD_INIT'

=head1 NAME

run_comparison.sh - compare SLT libraries, to turn single library results into comparative results

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@izb.unibe.ch, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2009, 2010 held by Daniel Nilsson. The package is realesed for use under the Perl Artistic License.

=head1 SYNOPSIS

USAGE: C<run_compare.sh comparison.list [gene_feature_type(gene)]>

=head1 DESCRIPTION

Compare splice sites and expression levels between SLT libraries.

Works on the output from F<run_maq.sh>, in particular F<counts.tab> files, 
and a short ordered list file describing the libraries to compare.

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

The script expects:

=over 4

=item *

an ordered comparison list

=item *

F<counts.tab> files, named in the comparison list, as produced by the F<run_maq.sh> part of the pipe

=item *

tag fastq files for libraries referenced (for lib sizes)

=item *

several of the other "sladdlampa"-scripts in $BINDIR

=item *

indirectly, bioPerl in the perl library path (and a perl 5.8.x interpreter at /usr/bin)

=back

=head1 SHELL ENVIRONMENT VARIABLES

Several environment variables influence the pipeline behaviour.

=over 4

=item BINDIR [path (~/install/sladd)]          

Directory where the rest of the SLADD pipeline lives.

=item forceupdate [<yes|no> (no)]
                     
Force rerun of all analysis regardless of 
modification times.

=item usekegg [<yes|no> (no)] 

Set to C<yes> to perform KEGG pathway specific analysis set. Note that the default is C<no>.

=back

=head1 OPTIONS AND ARGUMENTS

=over 4

=cut

POD_INIT

# check input environment variables and set unset ones to default values

# BINDIR holds the scripts/binaries for this pipeline
if [ -z "$BINDIR" ]
then   	
	#BINDIR=~/src/sladd
	BINDIR=~/install/sladd
fi

export PERL5LIB=$PERL5LIB:$BINDIR

# MAQ alignment quality cutoff NOT RESPECTED YET
if [ -z "$alqt" ]
then
    alqt=30
fi

# single mapping filtered files to be used for tab generation NOT RESPECTED YET
if [ -z "$single" ]
then
    single=yes
fi

if [ -z "$usekegg" ]
then
    usekegg=no
fi

# uses pipelinefunk.sh for needsUpdate, registerFile etc.

. $BINDIR/pipelinefunk.sh

# CALLED will contain a complete pathname for this script
CALLED=$0

if [ $# -lt 1 ]
then
	perldoc $CALLED
        echo "USAGE: ${0##*/} tabs_to_compare.ordered.list [gene_gff_type] [kegg_organism]"
        exit 1
fi

: <<'POD_ARG'

=item first argument: tabs_to_compare.ordered.list

first argument : F<tabs_to_compare.ordered.list> 

File with ORDERED list of tabs to compare. See sladd_howto.pdf for an example.

Read tablistfile: columns as follow.

libname	tabfilename

=cut

POD_ARG

tablistfile=$1 # file with ORDERED list of tabs to compare

if [ ! -e "$tablistfile" ]
then
    echo "Tab list file $tablistfile not found."
#    echo "USAGE: ${0##*/} tabs_to_compare.ordered.list [gene_gff_type] [kegg_organism]"
    perldoc $CALLED
    exit
fi

: <<'POD_ARG'

=item second argument, optional: genetype

second argument, optional: genetype : GFF feature type to associate tags with (gene, CDS, mRNA or such).
Default: gene.

=cut

POD_ARG

genetype=gene

if [ $# -eq 2 ]
then
    genetype=$2
fi

: <<'POD_ARG'

=item third argument, optional

third argument: KEGG organism abbreviation for organism in question. Default: tbr.

=back

=cut

POD_ARG

# KEGG organism abbreviation for organism in question
org=tbr
if [ $# -eq 3 ]
then
    org=$3
fi

# read tablistfile: columns as follow
# libname	tabfilename

unset libsize
unset libname
unset libtabfilename

libnr=0;

# extract lib information

# TODO: size counting is slow and could well only be made conditional on an update of the tablist, 

for libtab in `cut -f1 $tablistfile`;
do
#    echo $libtab
    libtabfilename[${libnr}]=$libtab;

    tagfile=${libtab%%_vs_*}.fastq

    alignmentprogramtmp=${libtab%%.counts.tab}
    
    alignmentsingle=${alignmentprogramtmp##*.}
    echo $alignmentsingle
    if [[ "$alignmentsingle" == *map ]]
    then
	alignmentprogramtmp=${libtab%%map.counts.tab}
	alignmentprogram=${alignmentprogramtmp##*.}
	alignmentsingle="multi"
    else
	# q30 or single after last dot before .counts.tab
	alignmentprogramtmp=${libtab%%map*}
	alignmentprogram=${alignmentprogramtmp##*.}
    fi

    echo $alignmentprogram
    echo $alignmentsingle

    libsize[${libnr}]=`grep -c ^@ ${tagfile}`
    libname[${libnr}]=`grep $libtab $tablistfile |awk '{ print $2 }'` 

    if [ ! -e "$libtab" ]
    then
	echo "Tab file $libtab not found (for lib ${libname[${libnr}]}, ordered ${libnr})!"
	exit;
    else
	echo "Tab file $libtab INDEED found (for lib ${libname[${libnr}]}, ordered ${libnr})!"	
    fi

    libnr=$(( $libnr + 1 ))
done

#for row in `awk '{ print $1,"+",$2,"+",$3; }' $tablistfile`; do 
#    ${libname[${libnr}]}=`echo $row |cut -f1 -d+`
#    ${libtabfilename[${libnr}]}= `echo $row |cut -f2 -d+`
#    ${libdesc[${libnr}]}=`echo $row |cut -f3 -d+`    
#    libnr=$(( $libnr + 1 ))
#done

updates=no

#numlibs=${#libtabfilename}
numlibs=$libnr

summary=${tablistfile}.summary

libnr=0
declare -a comptabs

echo There were $numlibs libs.

for ((liba=0;liba<$numlibs;liba++))
do
    libboffset=$(( $liba + 1 ))
#    echo a is $liba -- libb offset $libboffset
    for ((libb=$libboffset;libb<$numlibs;libb++))
    do
	comptab=${libname[${liba}]}_vs_${libname[${libb}]}.${alignmentprogram}.map.${alignmentsingle}.counts.compare
	if [ "$alignmentsingle" == "multi" ] 
	then
	    comptab=${libname[${liba}]}_vs_${libname[${libb}]}.${alignmentprogram}.map.counts.compare
	fi

	if needsUpdate $comptab ${tablistfile} ${libtabfilename[${liba}]} ${libtabfilename[${libb}]} $BINDIR/compare_counts_tabs.pl
	then
	    echo Compare by: $BINDIR/compare_counts_tabs.pl -o $comptab -M ${libsize[${liba}]} -N ${libsize[${libb}]} ${libtabfilename[${liba}]} ${libtabfilename[${libb}]}
	    registerFile ${comptab} result
	    $BINDIR/compare_counts_tabs.pl -o $comptab -M ${libsize[${liba}]} -N ${libsize[${libb}]} ${libtabfilename[${liba}]} ${libtabfilename[${libb}]}
	    updates=yes

	    registerFile ${comptab}.reg temp
	    perl -ne 'my $section=0; my $print_section=3; while(my $r=<STDIN>) {chomp $r; if($r=~/^\*\*\*/) { $section++ } if($print_section == $section) {print $r,"\n";}}' <$comptab > ${comptab}.reg
	    perl -ne 'my $section=0; my $print_section=4; while(my $r=<STDIN>) {chomp $r; if($r=~/^\*\*\*/) { $section++ } if($print_section == $section) {print $r,"\n";}}' < $comptab |grep -v "\#\#\#" >> ${comptab}.reg 
	    
	else
	    echo Compare ${libtabfilename[${liba}]} vs ${libtabfilename[${libb}]}: not done, already up to date.
	fi
	
	comptabs[${libnr}]=$comptab

	libnr=$(( $libnr + 1 ))
    done

    # list genes with interesting alternative splicesites
    # we initially screened for genes with both internal and external sites
    altspliceanalysis=${libtabfilename[${liba}]%%tab}altsplicestat.tab
    if needsUpdate $altspliceanalysis ${libtabfilename[${liba}]} $tablistfile $BINDIR/altsplice_analysis.pl
    then
	# -i 900 -u 2000 -f 0.6
        # the -A option now gives all genes passing the major site maximum fraction threshold with at least two alternative sites, valid according to the filter criteria.
	registerFile $altspliceanalysis result
	$BINDIR/altsplice_analysis.pl -A -o $altspliceanalysis -M ${libsize[${liba}]} -i 50000 -u 2000 -c 0 -f 0.6 ${libtabfilename[${liba}]}
	updates=yes
    fi

done

# co-regulation?  generalise beyond three samples... 
# unique features, like expressed only in, or upregulated only in.. 
# genes that stay the same, but are on (in all/at least one sample) 

# update the comparison summary - not much in there at the moment.
if needsUpdate $summary ${comptabs[@]} $tablistfile
then
    rundate=`date`
    echo "Comparison with ${tablistfile} ($rundate)" > $summary
    registerFile $summary result

    for comptab in ${comptabs[@]}
    do
	# count total number of up/down reg genes
	echo -n "${comptab%%.map.counts.compare} " >> $summary
	grep -v "^\*\|^\=\|^\#" ${comptab}.reg|wc -l >> $summary 
	grep -v "^\*\|^\=\|^\#" ${comptab}.reg |awk 'BEGIN {up=0; down=0;} ($2 !="" && $5>$10) { up = up +1 } ($2 !="" && $5<$10) { down = down + 1 } END { print "significantly up: " up ", significantly down: " down "." }  ' >> $summary
	grep -v "^\*\|^\=\|^\#" ${comptab}.reg |awk 'BEGIN {up=0; down=0;} ($2 !="" && ($5+1)/($10+1) > 2 ) { up = up +1 } ($2 !="" && ($5+1)/($10+1) < 0.5) { down = down + 1 } END { print "Of the significantly changed genes, " up " are twofold up but " down " twofold down." }' >> $summary
    done

    updates=yes
fi

: <<'POD_FUNCTION'

=head1 FUNCTIONS

=over 4

=item updateNormtab(mynormtab, sign)

Function to update the norm tab, differs slightly for getting significant only (second arg).
Second argument sign is yes or no, with significant only output produced when set to yes.

E.g. C<updateNormtab $normsigntab yes>

=back

=cut

POD_FUNCTION

# function to update the norm tab, differs slightly for getting significant only (second arg)
function updateNormtab()
{
    local mynormtab=$1
    local sign=$2

    if [ "$sign" = "yes" ] 
    then	
        # get all the genes that have sign difference between any two samples
	for compreg in ${comptabs[@]/compare/compare.reg}
	do
	    registerFile ${compreg}.list temp
	    grep -v "^#\|\*" $compreg |awk '($2 !="") { print $1 }' > ${compreg}.list;
	done


	cat ${comptabs[@]/compare/compare.reg.list} |sort -n |uniq > ${tablistfile}.reg.list
	registerFile ${tablistfile}.reg.list temp

        # individual columns for the norm tab, only genes listed above that have sign diff btw any two samples
	for tab in ${libtabfilename[@]}
	do 
	    on_exp_list=${tab%%.counts.tab}.on_exp_list.tab

	    registerFile ${on_exp_list} temp
	    $BINDIR/grep_big_list.pl -w -f ${tablistfile}.reg.list $tab > ${on_exp_list}

	    registerFile ${on_exp_list%%.tab}.name.tab temp
	    awk '{print $1}' ${on_exp_list} > ${on_exp_list%%.tab}.name.tab
	    registerFile ${on_exp_list%%.tab}.norm.tab temp
	    awk '{print $5}' ${on_exp_list} > ${on_exp_list%%.tab}.norm.tab
	done
    else
	for tab in ${libtabfilename[@]}
	do 
	    registerFile ${tab%%counts.tab}name.tab temp
	    awk '($2 >= 0) {print $1}' ${tab} > ${tab%%counts.tab}name.tab
	    registerFile ${tab%%counts.tab}norm.tab temp
	    awk '($2 >= 0) {print $5}' ${tab} > ${tab%%counts.tab}norm.tab
	done
    fi

    # norm tab header
    registerFile $mynormtab result
    echo Name ${libname[@]} Desc| perl -ne 's/ /\t/g; print;' > ${mynormtab}
    
    if [ "$sign" = "yes" ] 
    then
        # gene descriptions
	firstfile=${libtabfilename[0]/counts.tab/on_exp_list.tab}
    else
	firstfile=${libtabfilename[0]}
    fi
    registerFile ${firstfile}.desc temp
    awk '($2 != "") { print $13}' $firstfile > ${firstfile}.desc
#    awk '($2 != "") { print $13}' $firstfile > ${firstfile}.dir
#    awk '($2 != "") { print $13}' $firstfile > ${firstfile}.start

    # make one big tab
    if [ "$sign" = "yes" ] 
    then
	$BINDIR/join_columns.pl -H ${libtabfilename[0]/counts.tab/on_exp_list.name.tab} ${libtabfilename[@]/counts.tab/on_exp_list.norm.tab} ${firstfile}.desc >> ${mynormtab}
    else
	$BINDIR/join_columns.pl -H ${libtabfilename[0]/counts.tab/name.tab} ${libtabfilename[@]/counts.tab/norm.tab} ${firstfile}.desc >> ${mynormtab}
    fi

    registerFile ${mynormtab}.stats.R temp
    # standard dev, average on each gene, using R
    cat > ${mynormtab}.stats.R <<-STATS
	normtab <- read.table("${mynormtab}", header=TRUE)
	normtab\$mean<-apply(normtab[2:(2-1+${numlibs})],1,mean)
	normtab\$sd<-apply(normtab[2:(2-1+${numlibs})],1,sd)
	write.table(normtab[c(1:(2-1+${numlibs}),(2-1+${numlibs}+2):(2-1+${numlibs}+3),(2-1+${numlibs}+1))], file="${mynormtab%%tab}stats.tab", sep="\t",quote=F)
STATS
    registerFile ${mynormtab%%tab}stats.tab result
    R CMD BATCH --no-save --slave ${mynormtab}.stats.R /dev/null
#    R CMD BATCH ${mynormtab}.stats.R 

    updates=yes
}

normsigntab=${tablistfile}.sign.norm.tab
if needsUpdate $normsigntab $tablistfile $libtabfilename[@] $comptabs[@] $BINDIR/grep_big_list.pl $BINDIR/join_columns.pl
then
    updateNormtab $normsigntab yes
fi

normtab=${tablistfile}.norm.tab
if needsUpdate $normtab $tablistfile $libtabfilename[@] $comptabs[@] $BINDIR/grep_big_list.pl $BINDIR/join_columns.pl
then
    updateNormtab $normtab no
fi

# scatterplot, all on all (also outputs additional summary statistics..)
normtabscatterplot=${normtab%%tab}scatterall.pdf
if needsUpdate $normtabscatterplot $normtab $normsigntab $tablistfile
then

    # use R to scatterplot in a numlibs by numlibs matrix, with labels on the diagonal, and all dots above, significants below 
    # one could also do pairwise significants as before
    registerFile ${normtab}.scatterplot.all.R temp
    cat > ${normtab}.scatterplot.all.R <<-SCATTERALL
	normtab <- read.table("${normtab}", header=TRUE)
        normsigntab<-read.table("${normsigntab}", header=TRUE)
	pdf("$normtabscatterplot")
	par(mar=c(2,2,0.2,0.2))
	
	nsquared<-${numlibs}*${numlibs}
	layout(matrix(c(1:nsquared),${numlibs},${numlibs},byrow=TRUE))

        for(i in 1:nsquared) {
	        colnr<-i%%${numlibs}
                if (colnr == 0) {
                      colnr<-${numlibs}
		}
                rownr=ceiling(i/${numlibs})
		if (colnr==rownr) {
			plot.new(); 
        		plot.window(c(0,1), c(0,1), log = "");
        		text(0.5,0.5,substitute(bold(thelabel), list(thelabel = names(normtab)[1+rownr])),adj=c(0.5,0.5),cex=2)
		} else {
			xlogged=log10(normtab[,(1+colnr)])
			ylogged=log10(normtab[,(1+rownr)])
			xlogged[xlogged==-Inf] <- 0
			ylogged[ylogged==-Inf] <- 0

			xsignlogged=log10(normsigntab[,(1+colnr)])
			ysignlogged=log10(normsigntab[,(1+rownr)])
			xsignlogged[xsignlogged==-Inf] <- 0
			ysignlogged[ysignlogged==-Inf] <- 0
                        if(colnr>rownr) {
				plot( xlogged, ylogged, pch=".",cex=1.5, xlab=names(normtab)[1+colnr], ylab=names(normtab)[1+rownr],xlim=c(0,4), ylim=c(0,4)) 
			} else {
				plot( xsignlogged, ysignlogged, pch=".",cex=1.5, xlab=names(normsigntab)[1+colnr], ylab=names(normsigntab)[1+rownr],xlim=c(0,4), ylim=c(0,4)) 
			}
			spearmancor <- cor.test(normtab[,1+colnr], normtab[,1+rownr], method="spearman")
			pearsoncor <- cor.test(normtab[,1+colnr], normtab[,1+rownr], method="pearson")
			linreg <- summary(lm(normtab[,1+colnr]~normtab[,1+rownr]))

			corstr <- sprintf("Overall TPM norm correlations %s to %s; Spearman rho=%.2f (p<=%.4f), Pearson r=%.2f (p<=%.4f), Linear reg model R2=%.3f\\n",names(normtab)[1+colnr], names(normtab)[1+rownr], spearmancor\$estimate, spearmancor\$p.value,pearsoncor\$estimate,pearsoncor\$p.value,linreg\$adj.r.squared)

			cat(corstr,file="${summary}", sep="\\n", append=TRUE)

		}
	}
	dev.off()
SCATTERALL

#R CMD BATCH --no-save --slave ${normtab}.scatterplot.all.R /dev/null
    R CMD BATCH ${normtab}.scatterplot.all.R

    updates=yes
fi

if [ "$usekegg" = "yes" ] 
then
   
    # get per organism KEGG pathway files if old or missing

    echo -n Check for KEGG organism updates...

    seconds_in_two_days=$(( 60 * 60 * 24 * 2))

    update_pathways=no

    org_kegg_list=${org}.list
    if needsUpdateTimeBased ${org}.list $seconds_in_two_days
    then
	wget -c ftp://ftp.genome.jp/pub/kegg/pathway/organisms/${org}/${org}.list
	update_pathways=yes
	updates=yes
    fi
    
    if needsUpdateTimeBased map_title.tab $seconds_in_two_days
    then
	wget -c ftp://ftp.genome.jp/pub/kegg/pathway/map_title.tab 
	update_pathways=yes
	updates=yes
    fi

    echo Check for KEGG organism updates done.

    # make a set of gene lists for each pathway in a subdir called pathway
    org_gene_list=${org}_kegg.list.gene_first
    
    pathwaylist=${org}.kegg.pathways

    if needsUpdate $org_gene_list $org_kegg_list
    then
	registerFile ${org_kegg_list}.${org} temp
	grep $org\: $org_kegg_list > ${org_kegg_list}.${org}

 	registerFile $org_gene_list temp
	awk '{ print $2,"\t",$1;}' ${org_kegg_list}.${org} |sed -e 's/'$org'://; s/path://;'  > $org_gene_list

	update_pathways=yes
	updates=yes
   fi 
    
    if needsUpdate $pathwaylist $org_gene_list
    then
	registerFile ${pathwaylist} temp
	awk '($2 != "") {print $2}' $org_gene_list |sort|uniq > ${pathwaylist}
	update_pathways=yes
	updates=yes
    fi

    if [ ! -d pathways ]
    then
	registerFile pathways temp 
	mkdir pathways
	update_pathways=yes
	updates=yes
    fi

    if [ "$update_pathways" = "yes" ] 
    then
	for pathway in `cat ${pathwaylist}` ; 
	do 
	    grep $pathway $org_gene_list | awk '{print $1}' > pathways/$pathway.${org}_genes ;
	done
	
	for pathway in `cat ${pathwaylist} | sed -e 's/'${org}'//g;'`;
	do
	    grep $pathway map_title.tab > pathways/${org}${pathway}.map_title
	done
	updates=yes
    fi

    pathwaysummary=${tablistfile}.pathway.summary.tab

    if needsUpdate $pathwaysummary ${tablistfile} ${libtabfilename[@]} 
    then
	for tabfile in ${libtabfilename[@]} 
	do
	    for pathway in `cat ${pathwaylist}` ; 
	    do
		grep -w -f pathways/${pathway}.${org}_genes $tabfile > pathways/${tabfile%%.counts.tab}.$pathway.counts.tab 	   
	    done
	done

	# sum up total counts per pathway
	registerFile $pathwaysummary result
	echo PathID Pathway ${libname[@]} > $pathwaysummary
	for pathway in `cat ${pathwaylist}` ;
	do 
	    cat pathways/${pathway}.map_title|perl -ne 'chomp; print;' 
	    for tabfile in ${libtabfilename[@]}
	    do 
		awk 'BEGIN {sum=0}; {sum = sum + $5} END {printf "\t%i",sum}' < pathways/${tabfile%%.counts.tab}.$pathway.counts.tab
	    done
	    echo 
	done >> $pathwaysummary

	updates=yes
    fi

    foldchangetab=${pathwaysummary%%.summary.tab}.foldavg.summary.tab
    if needsUpdate $foldchangetab $normtab $tablistfile
    then

	# calculate fold changes per gene relative to the first lib
	# (use joined compare files!)
	# and average per pathway
	for pathway in `cat ${pathwaylist}`
	do
	    pathwaynormtab=pathways/${normtab%%.norm.tab}.$pathway.norm.tab
#	    echo Name ${libname[@]} Desc| perl -ne 's/ /\t/g; print;' > ${pathwaynormtab}
	    registerFile $pathwaynormtab temp
	    if [ "$org" == "tcr" ]
	    then

		# KEGG decided not use abbreviated pathway ids of
		# tc. That won't work with grep -w. And non-w is a
		# clear no-no when dealing with gene names ending in
		# .10, .100 and so on.

		awk '{print "Tc00.1047053"$1}' <pathways/${pathway}.${org}_genes >pathways/${pathway}.${org}_genes.edit
		mv pathways/${pathway}.${org}_genes.edit pathways/${pathway}.${org}_genes
	    fi

	    startcol=3
	    endcol=$(( $numblibs + 1 ))
	    
#	    perl -ne 'chomp; $r = $_; @c = split(/\t/,$r); printf("$r"."\t%.2f"x.('${numlibs}'+1)."\n", '"$c["{$startcol..$endcol}"]/$c[2]"');'

	    # compute fold changes per and sum up and average over each fold change, 
	    # note: set denominators of 0 to 0.9 ("fold change over 0 is slightly more impressive than fold change over 1") to avoid div by 0.
	    registerFile ${pathwaynormtab%%.norm.tab}.fold_change_over_prev_stage.tab temp
	    perl -e 'my $nfolds=0;
                     my @c;
                     my @foldsum; my @colsum;
                     while (<STDIN>) { 
                       chomp; $r = $_; 
                       @c = split(/\t/,$r);
                       $patternstr=""; 
                       $nlibs='${numlibs}'; 
                       $counts= $nlibs-1; 
                       $patternstr .= "\t%.2f" x $counts; 
                       $patternstr .= "\n"; 
		       $firstcol=1;
                       for (my $col=$firstcol; $col<($nlibs+$firstcol); $col++) {
                         $den=$c[$col];
                         if($den==0) {
			 	     $den=0.9;
			 };
			 if($col == $nlibs) {
			   $nom=$c[$firstcol];
                         } else {
		           $nom = $c[$col+1];
                         }
                         $fold[$col-1] = $nom/$den;
			 $foldsum[$col-1]+=$fold[$col-1];
                         $colsum[$col-1]+=$c[$col];
                         $nfolds++; 
		       } 
		       print $r,sprintf $patternstr, @fold;
		     } 
		     $rtabs=scalar(@c)-2;
		     @foldavg = map { $_/$nfolds } @foldsum;
		     print "AverageFoldChange\t", join("\t", @colsum),"\t", join("\t", @foldavg),"\n";' < $pathwaynormtab > ${pathwaynormtab%%.norm.tab}.fold_change_over_prev_stage.tab
	done

       	# same as pathway normtab plus the extra divs..
	# construct header
        #  comparisons made 
	for ((libden=0;libden<$numlibs;libden++))
	do
	    if [ $libden == $(( $numlibs - 1 )) ]
	    then 
		libnom=0
	    else		
		libnom=$(( $libden + 1 )) 
	    fi
	    
	    libvsexp[$libden]=${libname[${libnom}]}"vs"${libname[${libden}]}
	done

#  "PathwayNr","Pathway","PathwayId",'$' "LS","SS","PC","SSvsLS","PCsSS","LSvsPC"),"\n"

	registerFile $foldchangetab result
	perl -e 'print join("\t", @ARGV)' "PathwayNr" "Pathway" "PathwayId" ${libname[@]} ${libvsexp[@]} "\n" > $foldchangetab

	for pathway in `cat ${pathwaylist}`
	do
	    cat pathways/${pathway}.map_title|perl -ne 'chomp; print; print "\t"'
	    pathwaynormtab=pathways/${normtab%%.norm.tab}.$pathway.norm.tab
	    grep ^AverageFoldChange ${pathwaynormtab%%.norm.tab}.fold_change_over_prev_stage.tab |sed -e 's/AverageFoldChange/'${pathway}'/g;'
	done >> $foldchangetab

	updates=yes
	# then 
    fi
fi

if [ "$updates" = "no" ]
then
    rundate=`date`
    echo "Project was already up to date. Nothing to do. Please update date on appropriate input/program/intermediate results files to force results updates. ($rundate)"
else
    rundate=`date`
    echo "Project has been brought up to date. ($rundate)"
fi
exit

: <<'POD_FILES'

=head1 FILES

Please see the sladd_howto.pdf for file descriptions.

=cut

POD_FILES


# scatterplots, pairwise

#function scatterplot() 
#{
#    local comptabfile=$1;
#    local libaname=$2
#    local libbname=$3
#    local regsignnorm=$1.reg.sign_norm_only
#    cat > ${comptabfile}.scatterplot.R <<-SCATTERPLOT
#    tab <- read.table(${regsignnorm}, header=TRUE)
#    pdf("${comptabfile}.scatterplot.pdf")
#    plot( log10(tab\$V1), log10(tab\$V2),pch=".",cex=1.5, xlab="${libaname}", ylab="${libbname}",main="$libaname vs $libbname TPM" , xlim=c(0,4), ylim=c(0,4))   
#    dev.off()
#
#SCATTERPLOT
#
#    R CMD BATCH --no-save --slave ${comptabfile}.scatterplot.R /dev/null
#}

#scatterplots=${comptabs[@]/compare/plot}

#for scatterplot in ${scatterplots[@]}
#do
   # actually there is only one comp tab to check dates on for each scatterplot

#    libaname=${scatterplot%%_vs_*}
#    libbnametmp=${scatterplot##*_vs_}
#    libbname=${libbnametmp%%.${alignmentprogram}.map.counts.compare}
#    comptab=${plot/plot/compare}

 #   if needsUpdate $scatterplot $comptab
 #   then
#	scatterplot $comptab $libaname $libbname
#    fi
#done

# extract all on all up/down summary from pairwise
# --- also do an ordered progression on these

# r-matrix from list order norm exp levels
#     plots on selected genes
#     overall stats from the progressions
#     profile clustering (via bioconductor?)

# extract all on all sign diff splice from pairwise

# pathways

# also do any of this for major transcript only??

