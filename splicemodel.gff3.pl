#!/usr/bin/perl -w -I. -I/Users/daniel/install/splicemodel
# Daniel Nilsson, 2002-09 -- 2002-10
#
# simple splicemodel to separate genes conforming to current opinion on tryp trans-splicing
# and find potential alternative trans-splice sites
#
# Lm: leBowitz et al 1993: cds-stop-utr-|polyA------polypy---AG-----ATG-----stop
#     polya - pyag precisely specified as 200-500 nt. =/
#     Distance constrained to 100-400 in Clayton, 2002 review.
# Tb: 
# 
# Tc: 
# 
# begin and end of transcripts? What happens at the ultimate ends of a directional cluster?
# Is there a requirement on "junk-sites"/"junk-genes" after the last gene? "Reason" for duplication ("due to" duplications)??
# There are also several other separate mechanisms involved:
# * polI promoters - hows the polyA/SL addition there? On the genes nearby?
# Transcription start & stop: few, and probably not so easy to catch -- but this will be a possible inroad
# add up features in common table? insert keywords in sequence? Keep separate tables??

# "most distal sites" as suggested by Clayton et al fits well with using perls greedy patterns - although the imposed distance 
# constraints become the edge limits rather than some actual biological feature..

# Please check length of utr seq - these are afaik fairly variable ~4-500 bp, but they surely influence the model, since the bioresults 
# deal with translation start and term rather than "transcript" start and term. (And *transcription* start and term remain complete unknowns.)
#
# 5'UTR; more like 50-150 nt, 3'UTR much longer -- up to about 3kb in the extreme (PABP) case.
#
# It is reasonable to assume that the shortest polypy-ag-distance is the wanted one, *but* there is certainly a possiblility that suitable nearness to the atg is preferable

# version 0.1 : print sequential list of possible sites
# version 0.2 : state-model, checking which are ok -- needs two passes for forward/reverse
# version 0.3 : multiple ok sites checked for "alternative splicing"

# 030228 NOTE : coordinate choice on 3p end fasta export are really off - they point to gene start instead... =)

my $version = "060210";

use seqfunk;

sub polypyagscan;
sub polypyag_cds_start_check;
sub polyasite_check_genes;
sub get_aaa_around;
sub extremes;
sub check_uorf;

sub output_uorfs;
sub write_fasta_utrs;
sub write_gff_utrs_for_cds;
sub write_uorf_statistics;
sub output_extremes_stats;

sub usage;

# Default values for model parameters (last adjusted 040313 to fit Tc 5p ESTs alignment data)
my $cds_5p_ag_atg_dist = 1850; # essentialy length of 5'utr - SL as seen by ribosome # tb at least 850
#my $min_poly_py_length = 15; # 16, normally # T
#my $max_poly_py_mismatches = 3; 
my $min_poly_py_length = 10; # 16, normally # T
my $max_poly_py_mismatches = 1;
my $min_qualified_pyrun = 8;

my $exp_stop_polya_ag_distance =  5040; # 2500; maximum length of 3'UTR + polya-dst_polypyag (distance cds stop to pyag)

my $exp_polypy_ag_distance = 261; # length btw polypy-stretch and the downstream splice-site -- 261 in Tb, 105 in Tc
my $polya_pyag_dist = 100; # distance btw the 5p of the downstream (polypy-connected) ag and upstream polya-site (essentialy the pseudo-intron length)
my $polya_aaa_scanning_window = 30; # half size of win (well, N, where winsize is 2N+1) in which to scan for tri-a for p-a site selection

my $DEBUG = 0;
my $WARNING = 1;

my $gff3outfh = *STDOUT;
my $uorfoutfh = *STDOUT;
my $fastaoutfh = *STDOUT;
my $reportoutfh = *STDOUT;

my $seqFileName = undef;
my $reportFileName = undef;
my $gffFileName = undef;
my $fastaFileName = undef;
my $uorfFileName = undef;
my $statisticsfile = undef;
my $gff_utr_file = undef;
my $rout = undef;
my $revcomp = 1;
my $fasta = 0;
my $gff_utrs = 0;
my $spliceok = 0;
my $report = 0;
my $long = 0;
my $short = 0;
my $check_uorf = 0;
my $no_aaa_check = 0;
my $downstream_gene_interferes_with_splicing = 0;
my $closest_a = 1;
my $aaa_scan = 0;
my $disable_dupe_checking = 0;
my $print_internal_alt_site = 0;
my $internal_alt_file_name = "";
my $splitstats = 0;

my $longest_uorf_free = 1;

# add "never include other cds in utr"?

# my $nr_of_ppy_found = 0;
my $nr_of_ppy_ag_found = 0;

while (my $arg = shift @ARGV) {
    if($arg eq "--debug") {
	$DEBUG = 1;
    } elsif ($arg eq "--infile") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No infile name given.");
	} else {
	    $seqFileName = $arg;
	}
    } elsif ($arg eq "--report") {
	$report = 1;
    } elsif ($arg eq "--reportfile") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No report file name given.");
	} else {
	    # reportfile implies report
	    $report = 1;
	    $reportFileName = $arg;
	}
    } elsif ($arg eq "--no-revcomp") {
	$revcomp = 0;
    } elsif ($arg eq "--no-warn") {
	$WARNING = 0;
    } elsif ($arg eq "--gff") {
	$gff = 1;
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No gff file name given.");
	} else {
	    $gffFileName = $arg;
	}
    } elsif ($arg eq "--fasta") {
	$fasta = 1;
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No fasta outfile name given.");
	} else {
	    $fastaFileName = $arg;
	}
    } elsif ( $arg eq '--gff-utrs' ) {
	$gff_utrs = 1;
	$arg = shift @ARGV;
	if ($arg eq '') {
	    usage("Error! No fasta outfile name given.");
	} else {
	    $gff_utr_file = $arg;
	}
    } elsif ($arg eq "--no-aaa-check") {
	$no_aaa_check = 1;
    } elsif ($arg eq "--polya-model") {
	$polyamodel = shift; 
	if ( $polyamodel eq 'fix' ) {
	    $no_aaa_check = 1;
	} elsif ( $polyamodel eq 'aaa') {
	    $no_aaa_check = 0;
	    $aaa_scan = 1;
	    $closest_a = 0;
       } elsif ( $polyamodel eq 'center' ) {
	   $closest_a = 1;
	   $aaa_scan = 0;
	   $no_aaa_check = 0;
       }
    } elsif ($arg eq '--coupled-polya') {
	$downstream_gene_interferes_with_splicing = 1;
    } elsif ($arg eq "--all-ppyag") {
	$disable_dupe_checking = 1;
    } elsif ($arg eq "--long") {
	$long = 1;
    } elsif ($arg eq "--short") {
	$short = 1;
    } elsif ($arg eq "--uorf") {
	$check_uorf = 1;
    } elsif ($arg eq "--uorf-file") {
	$arg = shift;
	if ($arg eq "") {
	    usage("Error! No uORF file name given.");
	} else {
	    $check_uorf = 1;
	    $uorfFileName = $arg;
	}
    } elsif ($arg eq "--print-int-alt-sites-file") {
	$arg = shift;
	if ($arg eq "") {
	    usage("Error! No internal alternative sites file name given.");
	} else {	    
	    $internal_alt_file_name = $arg;
	    $print_internal_alt_site = 1;
	}
    } elsif ($arg eq "--polypy-ag") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No polypy to downstream ag distance given.");
	} else {
	    $exp_polypy_ag_distance = $arg;
	}
    } elsif ($arg eq "--min-polypy-length") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No minimum polypyrimidine stretch length given.");
	} else {
	    $min_poly_py_length = $arg;
	}
    } elsif ($arg eq "--min-qual-polypy-run") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No minimum qualified polypyrimidine stretch length given.");
	} else {
	    $min_qualified_pyrun = $arg;
	}
    } elsif ($arg eq "--max-polypy-mismatches") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No maximum number of polypy mismatches given.");
	} else {
	    $max_poly_py_mismatches = $arg;
	}
    } elsif ($arg eq "--5pag-atg") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No 5pag-atg distance given.");
	} else {
	    $cds_5p_ag_atg_dist = $arg;
	}
    } elsif ($arg eq "--polya-pyag") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No pa-pyag distance given.");
	} else {
	    $polya_pyag_dist = $arg;
	}
    } elsif ($arg eq "--half-polya-win") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No half polya window length given.");
	} else {
	    $polya_aaa_scanning_window = $arg;
	}
    } elsif ($arg eq "--max-stop-ag-distance") { # FIX - obsolete name!
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No maximum distance btw stop codon and polypy-ag site given for polya-model.");
	} else {
	    $exp_stop_polya_ag_distance = $arg;
	}
    } elsif ($arg eq "--rout") {
	$rout = 1;
    } elsif ($arg eq "--splitstats") {
	$splitstats = 1;
    } elsif ($arg eq "--statistics") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No statistics file given.");
	} else {
	    $statisticsfile = $arg;
	}
    } else {
	usage("Unknown option: $arg");
    }
}

use Bio::Seq;
use Bio::SeqIO;

!defined($seqFileName) && usage("No input file name given (required)!");
-e $seqFileName or fatal("File not found $seqFileName.");
$rout && (defined($statisticsfile) || fatal("Cannot output R-format statistics file when no statistics file name is given."));
my $gbkin = Bio::SeqIO->new('-format' => 'genbank', -file => $seqFileName); # stdin-problem? check seqio doc..

# open gff3 outfile or STDOUT
if($gff && !defined($gffFileName)) {
    $DEBUG && debug("Writing gff to STDOUT."); # default behaviour
} elsif($gff) {
    local *OUTFILE;
    open(OUTFILE, ">$gffFileName") or fatal("Could not open $gffFileName for writing.");
    $gff3outfh = *OUTFILE;
}

# open uORF outfile or STDOUT
if($check_uorf && !defined($uorfFileName)) {
    $DEBUG && debug("Writing uORF gff to STDOUT."); # default behaviour
} elsif($check_uorf) {
    local *OUTFILE;
    open(OUTFILE, ">$uorfFileName") or fatal("Could not open $uorfFileName for writing.");
    $uorfoutfh = *OUTFILE;
}

# open fasta outfile or STDOUT
if($fasta && !defined($fastaFileName)) {
    $DEBUG && debug("Writing fasta to STDOUT."); # default behaviour
} elsif($fasta) {
    local *OUTFILE;
    open(OUTFILE, ">$fastaFileName") or fatal("Could not open $fastaFileName for writing.");
    $fastaoutfh = *OUTFILE;
}

if(!defined($reportFileName)) {
    $DEBUG && debug("Writing report to STDOUT."); # default behaviour
} elsif($report) {
    local *OUTFILE;
    open(OUTFILE, ">$reportFileName") or fatal("Could not open $reportFileName for writing.");
    $reportoutfh = *OUTFILE;
}

my $statisticsfh; # set to *STDOUT per default?
if(defined($statisticsfile) && ($splitstats == 0)) {
    $DEBUG && debug("Writing statistics to $statisticsfile.");
    open STATISTICS, ">$statisticsfile" or fatal("Could not open $statisticsfile for writing.");
    $statisticsfh = *STATISTICS;
}

if($print_internal_alt_site) {
    $DEBUG && debug("Writing internal alternative sites to file $internal_alt_file_name.");
    open ALTSITE, ">$internal_alt_file_name";    
}

my $seq;

$seq = $gbkin->next_seq();

$DEBUG && debug("Scanning for polypy");
polypyagscan(\$seq);

# make a feature list

my @cfeat; # considered features
my $nCDS = 0;
my $npPyAg = 0;

$DEBUG && debug("Construct initial feature list..");
FEAT:    foreach my $feat ( $seq->all_SeqFeatures() ) {
     if( ! defined($feat) ) {
	 # Skip feature-less entries, such as contig scaffolds..
	 next FEAT;
     } 
     $revcomp || $feat->strand==0 || $feat->strand==1 || next FEAT;

     ($feat->strand == 0) && ($feat->strand = 1);
     # get cds feature
     if ($feat->primary_tag eq "pPyAg") {
	 push @cfeat, $feat;
	 $npPyAg++;
     } elsif($feat->primary_tag eq "CDS") {
	 $nCDS++;
	 push @cfeat, $feat;
     }
}

my @dfeat; # features in discordance with model
my @ofeat; # feat in agreement with model 

if($nCDS == 0) {
    $WARNING && warning("WARNING: no CDS for this sequence ".($seq->id)." - has $npPyAg PyAG sites." );
}

if($npPyAg == 0) {
    $WARNING && warning("WARNING: no pPyAG foud for sequence ".($seq->id).", so the $nCDS CDSs are without sites.");
}

my @opafeat;
my @dpafeat;
my @pyfeat;
my @usedppy;
my @usedpappy;

# check the polypyag - cds start distances
$DEBUG && debug("Check cds start site..");
polypyag_cds_start_check(); 

# fill in extremes (longest and shortest utrs) on request
my %shortest5p;
my %shortest5ppy;
my %longest5p;
my %longest5ppy;

my %shortest3p;
my %shortest3ppy;
my %longest3p;
my %longest3ppy;

# ehhr.. Why not do this before? Simply rewrite pyag scan to be gene centric, and implement there directly..

if ($long || $short) { 
    $DEBUG && debug("Extreme length filtering 5p..");
    extremes_5p();
    extremes_clean_cds_5p();
    defined($statisticsfile) && output_extremes_stats();
}

# check polypyag - potential polya site
$DEBUG && debug("Check polya site..");
polyasite_check_genes();

if ($long || $short) { 
    $DEBUG && debug("Extreme length filtering 3p..");
    extremes_3p();
    extremes_clean_cds_3p();
    defined($statisticsfile) && output_extremes_stats();
}

# output results
$fasta && write_fasta_utrs();
$gff_utrs && write_gff_utrs_for_cds();

# check for upstream ORFs
my $uorffree;
my $utrs_uorf_checked;
my $uorf;
my @uorfcount;

if($check_uorf) {
    ($uorffree, $uorfolap, $utrs_uorf_checked, @uorfcount) = check_uorf();
}

$check_uorf && output_uorfs();

# $check_uorf && $rout && write_uorf_statistics();

if ($gff) { # GFF-output # longest/shortest here as well?
    map { print $gff3outfh gff3_string($$_)."\n" } @ofeat;     # output all ok cds feats
    map { print $gff3outfh gff3_string($$_)."\n" } @pyfeat;    # output all polypy features (with "ok" annotation added)
}

if ($report) { # model results summary + gff partitioned into different categories # report == 1
    if($cds_5p_ag_atg_dist) {
	print $reportoutfh "\n*** In summary, ".scalar(@dfeat)." cds features disagree with the polypy ag model, and ".scalar(@ofeat)." agree. \n";
	print $reportoutfh "\n*** Model parameters: exp_polypy_ag_distance=$exp_polypy_ag_distance max_poly_py_mismatches=$max_poly_py_mismatches".
	    " min_poly_py_length=$min_poly_py_length cds_5p_ag_atg_dist=$cds_5p_ag_atg_dist\n\n";
	print $reportoutfh "\n*** In total, $nr_of_ppy_ag_found ppyag sites were found over the entire sequence (forward AND reverse).\n";
	if ($check_uorf) {
	    print $reportoutfh "*** $uorffree 5pUTRs out of $utrs_uorf_checked checked were free from uORFs.\n";
	    # print "*** n cds checked have uorf free utrs.." -> use shortest...
	    print $reportoutfh "*** $uorfolap 5pUTRs show uORFs overlapping annotated CDS start.\n";

	    print $reportoutfh "*** The number 5pUTRs with 1, 2, 3, etc uAUGs was as follows: ";
	    for($i = 1; $i<@uorfcount; $i++) {
		defined($uorfcount[$i]) && print $reportoutfh $uorfcount[$i]; 		       
		print $reportoutfh ", ";
	    }
	    print $reportoutfh "\n";   
	}
	print $reportoutfh "\n*** The following cds features disagree with current model: \n";
	map { print $reportoutfh gff3_string($$_)."\n" } @dfeat;

	print $reportoutfh "\n*** The following cds features agree with current model: \n";
	map { print $reportoutfh gff3_string($$_)."\n" } @ofeat;
    }
    if($exp_stop_polya_ag_distance) {
	print $reportoutfh "\n*** In summary, ".scalar(@dpafeat)." cds features disagree with the polya-polypy model, and ".scalar(@opafeat)." agree. \n";
	print $reportoutfh "\n*** Model parameters: exp_stop_polya_ag_distance=$exp_stop_polya_ag_distance exp_polypy_ag_distance=$exp_polypy_ag_distance".
	    " max_poly_py_mismatches=$max_poly_py_mismatches min_poly_py_length=$min_poly_py_length cds_5p_ag_atg_dist=$cds_5p_ag_atg_dist\n\n";

	print $reportoutfh "\n*** The following cds features disagree with current model: \n";
	map { print $reportoutfh gff3_string($$_)."\n" } @dpafeat;

	print $reportoutfh "\n*** The following cds features agree with current model: \n";
	map { print $reportoutfh gff3_string($$_)."\n" } @opafeat;
    }
}

###### end main ###### ###### end main ###### ###### end main ###### ###### end main ###### ###### end main ######

sub write_uorf_statistics {

    # r-header
    print "length\tstart\tend\tterminated\toverlapping\n";
    
    # select a set of ppy feats
    my $selectedppyfeats;
    
    if ($long) {
	$selectedppyfeats = [values %longest5ppy];
    } elsif($short) {
	$selectedppyfeats = [values %shortest5ppy];

    } else {
	$selectedppyfeats = \@usedppy;
    }
    
    my %uORFfree;

    foreach my $polypyfeat (@$selectedppyfeats) {

	my ($cdsatgpos) = $$polypyfeat->each_tag_value('cds-atg');
	my $sladdpos;
	if ($long) {
	    $sladdpos = $longest5p{$cdsatgpos};
	} elsif($short) {
	    $sladdpos = $shortest5p{$cdsatgpos};
	} else {
	    ( ($$polypyfeat->strand == 1) && ($sladdpos = $$polypyfeat->end) ) || ($sladdpos = $$polypyfeat->start);
	    $DEBUG && debug("Setting sladdpos = $sladdpos, cdsatgpos $cdsatgpos.");
	}

	if($$polypyfeat->has_tag('uORF')) {
	    	    
	    my @uORFtag = $$polypyfeat->each_tag_value('uORF');    
	    
	    foreach my $uORFt ( @uORFtag ) {
		# print "DEBUG: $uORFt\n";
 
		%uORFtag_values = ( split(/(?:\s+|=)/, $uORFt) );
			 	
		print "$uORFtag_values{length}\t$uORFtag_values{start}\t$uORFtag_values{end}\t$uORFtag_values{terminated}\t";
		defined($uORFtag_values{overlapping}) && print "$uORFtag_values{overlapping}";
		print "\n";
		# Remember; several $$polypyfeat->add_tag_value('uORF', "start=$augp end=$endp terminated=0 overlapping=0 length=".($uORFend-$uORFstart));
	    }	    
	} else {
	    # no uORF on sl
	    #print "uORF-free \n";
	    if( $longest_uorf_free ) {
		my ($uORFfree) = $$polypyfeat->each_tag_value('uORFfree'); # length of TLR
		print "uORFfree cdsatg=$cdsatgpos ppystart=", $$polypyfeat->start, " ppyend=", $$polypyfeat->end, " $uORFfree\n";
		if (!defined($uORFfree{$cdsatgpos})) {
		    $uORFfree{$cdsatgpos} = $uORFfree; # length
		} elsif ( $uORFfree > $uORFfree{$cdsatgpos} ) {
		    $uORFfree{$cdsatgpos} = $uORFfree;
		}		
	    }
	}
    }
    
    if( $longest_uorf_free ) {
	foreach my $cdsatgpos (keys %uORFfree ) {
	    print "$cdsatgpos\t$uORFfree{$cdsatgpos}\n";
	}
    }
}

sub write_fasta_utrs {
    # push @ofeat, @usedppy;
    # @ofeat = sort { $a->strand <=> $b->strand || ($a->strand == 1 && $a->start <=> $b->start) || ($a->strand == -1 && $b->end <=> $a->end)} @ofeat;
    
    # foreach used ppy, write the sequence to the next cds

    my $selectedppyfeats;
    my %utrpos;
    
    if ($long) {
	$selectedppyfeats = [values %longest5ppy];
    } elsif($short) {
	$selectedppyfeats = [values %shortest5ppy];
    } else {
	$selectedppyfeats = \@usedppy;
    }

    foreach my $polypyfeat (@$selectedppyfeats) {
	my ($cdsatgpos) = $$polypyfeat->each_tag_value('cds-atg');
	my $sladdpos;
	if ($long) {
	    $sladdpos = $longest5p{$cdsatgpos};
	} elsif($short) {
	    $sladdpos = $shortest5p{$cdsatgpos};
	} else {
	    ( ($$polypyfeat->strand == 1) && ($sladdpos = $$polypyfeat->end) ) || ($sladdpos = $$polypyfeat->start);
	    $DEBUG && debug("Setting sladdpos = $sladdpos, cdsatgpos $cdsatgpos.");
	}

	$DEBUG && debug("Writing 5putr [ppy at ".$$polypyfeat->start."-".$$polypyfeat->end."] for direction ".$$polypyfeat->strand.".");
	my $fiveputr;
	if($$polypyfeat->strand == 1) {
	    $DEBUG && debug("Getting subseq from ".$sladdpos." to ".$cdsatgpos.".");
	    if(!defined($utrpos{$sladdpos."+".$cdsatgpos})) {
		$DEBUG && debug("getting the utr...");
		$fiveputr = $seq->subseq($sladdpos,$cdsatgpos);
		print $fastaoutfh ">5pUTR_".$sladdpos."-".$cdsatgpos;
		print $fastaoutfh "\n$fiveputr\n";
		$DEBUG && debug("writing $sladdpos - $cdsatgpos");
		$utrpos{$sladdpos."+".$cdsatgpos} = 1;
	    } else {
		$DEBUG && debug("Duplicate 5'UTR - ignoring it.");
	    }
	} else {
	    $DEBUG && debug("Getting subseq from ".$cdsatgpos." to ".$sladdpos.".");
	    if(!defined($utrpos{$sladdpos."-".$cdsatgpos})) {
		$DEBUG && debug("getting the rc utr...");
		$fiveputr = revcomp($seq->subseq($cdsatgpos,$sladdpos));
		print $fastaoutfh ">5pUTR_".$cdsatgpos."-".$sladdpos;
		$DEBUG && debug("writing $cdsatgpos - $sladdpos");
		print $fastaoutfh "\n$fiveputr\n";		
		$utrpos{$sladdpos."-".$cdsatgpos} = 1;
	    } else {
		$DEBUG && debug("Duplicate 5'UTR - ignoring it.");
	    }
	}
    }

    if ($long) {
	$selectedppyfeats = [values %longest3ppy];
    } elsif($short) {
	$selectedppyfeats = [values %shortest3ppy];
    } else {
	$selectedppyfeats = \@usedpappy;
    }

    foreach my $polypyfeat (@$selectedppyfeats) {
	my ($cdsendpos) = $$polypyfeat->each_tag_value('cds-end');
	my ($polyapos) = $$polypyfeat->each_tag_value('polya');
	
	$DEBUG && debug("Writing 3putr [".$$polypyfeat->start."-".$$polypyfeat->end."] for direction ".$$polypyfeat->strand.".");

	my $threeputr;
	print $fastaoutfh ">3pUTR_";
	if($$polypyfeat->strand == 1) {
	    $DEBUG && debug("Getting subseq from ".$cdsendpos." to ".$polyapos.".");
	    $threeputr = $seq->subseq($cdsendpos, $polyapos);
	    print $fastaoutfh $cdsendpos."-".$polyapos;
	} else {
	    $DEBUG && debug("Getting subseq from ".$polyapos." to ".$cdsendpos.".");
	    $threeputr = revcomp($seq->subseq($polyapos, $cdsendpos));
	    print $fastaoutfh $polyapos."-".$cdsendpos;
	}
	print $fastaoutfh "\n$threeputr\n";
    }
}

sub write_gff_utrs_for_cds {
    open GFF_UTR_OUT, ">$gff_utr_file";

    my $seq_name = $seq->id;
    
    foreach $accepted_cds_feat ( @ofeat ) {

	my $strand = $$accepted_cds_feat->strand;
	my $cdsatgpos = ($strand>0)?$$accepted_cds_feat->start:$$accepted_cds_feat->end;

	foreach my $sladdpos ( $$accepted_cds_feat->each_tag_value('SLadd') ) {
	    my $strandchar = ($strand == 1) ? '+' : '-';
	    my $comment = ""; # $$polypyfeat->all_tags();

	    if($strand == 1) {
		# the UTR start is one nt closer to the AUG than the AG-add site. Also, obviously, the UTR ends on the base just before the AUG.
		# This actually means we get GFFs where the end is smaller than start when facing a 0-len UTR, that has its start pos on the 0: A in aug and at the base before; -1. 
		# That makes length == 0, obviously start - end + 1 == -1 - 0 + 1 == 0.
		print GFF_UTR_OUT "$seq_name\tpapya_$version\t5pUTR\t".($sladdpos+1)."\t".($cdsatgpos-1)."\t.\t$strandchar\t.\t$comment\n"; # don't include the initial and end A's on the "UTR"
	    } else {
		print GFF_UTR_OUT "$seq_name\tpapya_$version\t5pUTR\t".($cdsatgpos+1)."\t".($sladdpos-1)."\t.\t$strandchar\t.\t$comment\n"; # don't include the initial and end A's on the "UTR"
	    }
	}
    }
    
    foreach $accepted_cds_feat ( @opafeat ) {

	my $strand = $$accepted_cds_feat->strand;
	my $cdsendpos = ($strand>0)?$$accepted_cds_feat->end:$$accepted_cds_feat->start;

	foreach my $polyapos ( $$accepted_cds_feat->each_tag_value('polya') ) {
	    my $strandchar = ($strand == 1) ? '+' : '-';
	    my $comment = ""; # $$polypyfeat->all_tags();
	    if($strand == 1) {
		print GFF_UTR_OUT "$seq_name\tpapya_$version\t3pUTR\t".($cdsendpos+1)."\t$polyapos\t.\t$strandchar\t.\t$comment\n"; 
	    } else {
		# we don't want to include the last base of the stopcodon..
		print GFF_UTR_OUT "$seq_name\tpapya_$version\t3pUTR\t$polyapos\t".($cdsendpos-1)."\t.\t$strandchar\t.\t$comment\n";
	    }
	}	
    }

    close GFF_UTR_OUT;
}

sub extremes_5p {

   foreach my $polypyfeat (@usedppy) {
        # getting one cds only! ok, then thats enough for identifying "the" cds.
	my (@cdsatgpos) = $$polypyfeat->each_tag_value('cds-atg'); 

	for (my $i=0; $i < @cdsatgpos; $i++) { # a polypyfeat can be the nearest polyasite for several cds
	    if($$polypyfeat->strand == 1) {
		if(!defined($shortest5p{$cdsatgpos[$i]}) || $shortest5p{$cdsatgpos[$i]} < $$polypyfeat->end) {
		    $shortest5p{$cdsatgpos[$i]} = $$polypyfeat->end;
		    $shortest5ppy{$cdsatgpos[$i]} = $polypyfeat;
		}
		if(!defined($longest5p{$cdsatgpos[$i]}) || $longest5p{$cdsatgpos[$i]} > $$polypyfeat->end) {
		    $longest5p{$cdsatgpos[$i]} = $$polypyfeat->end;
		    $longest5ppy{$cdsatgpos[$i]} = $polypyfeat;
		}
	    } else {
		if(!defined($shortest5p{$cdsatgpos[$i]}) || $shortest5p{$cdsatgpos[$i]} > $$polypyfeat->start) {
		    $shortest5p{$cdsatgpos[$i]} = $$polypyfeat->start;
		    $shortest5ppy{$cdsatgpos[$i]} = $polypyfeat;
		}
		if(!defined($longest5p{$cdsatgpos[$i]}) || $longest5p{$cdsatgpos[$i]} < $$polypyfeat->start) {
		    $longest5p{$cdsatgpos[$i]} = $$polypyfeat->start;
		    $longest5ppy{$cdsatgpos[$i]} = $polypyfeat;
		}
	    }
	}
    }
}
sub extremes_3p {

    foreach my $polypyfeat (@usedpappy) {
	my (@cdsendpos) = $$polypyfeat->each_tag_value('cds-end');
	my (@polyapos) = $$polypyfeat->each_tag_value('polya');

	for (my $i=0; $i < @cdsendpos; $i++) { # a polypyfeat can be the nearest polyasite for several cds
	    if($$polypyfeat->strand == 1) {
		if(!defined($shortest3p{$cdsendpos[$i]}) || $shortest3p{$cdsendpos[$i]} > $polyapos[$i]) {
		    $shortest3p{$cdsendpos[$i]} = $polyapos[$i];
		    $shortest3ppy{$cdsendpos[$i]} = $polypyfeat;
		}		
		if(!defined($longest3p{$cdsendpos[$i]}) || $longest3p{$cdsendpos[$i]} < $polyapos[$i]) {
		    $longest3p{$cdsendpos[$i]} = $polyapos[$i];
		    $longest3ppy{$cdsendpos[$i]} = $polypyfeat;
		}
	    } else {
		if(!defined($shortest3p{$cdsendpos[$i]}) || $shortest3p{$cdsendpos[$i]} < $polyapos[$i]) {
		    $shortest3p{$cdsendpos[$i]} = $polyapos[$i];
		    $shortest3ppy{$cdsendpos[$i]} = $polypyfeat;
		}
		if(!defined($longest3p{$cdsendpos[$i]}) || $longest3p{$cdsendpos[$i]} > $polyapos[$i]) {
		    $longest3p{$cdsendpos[$i]} = $polyapos[$i];
		    $longest3ppy{$cdsendpos[$i]} = $polypyfeat;
		}
	    }
	}
    }
}

sub extremes_clean_cds_5p {

    if($short) {
	# shortest5p

	# remove all ppy feat SLaddForCDS and cds-atg tags..
	foreach my $polypyfeat ( @usedppy ) {
	    if( $$polypyfeat->has_tag('SLaddForCDS') ) {
		$$polypyfeat->remove_tag('SLaddForCDS');
		$$polypyfeat->remove_tag('cds-atg');
	    }
	}

	foreach my $accepted_cds_feat ( @ofeat ) {

	    # get the cds ATG pos (start or end depending on strand)
	    my $accepted_cds_feat_atg_pos;
	    
	    if( $$accepted_cds_feat->strand == 1 ) {
		$accepted_cds_feat_atg_pos = $$accepted_cds_feat->start;
	    } else {
		$accepted_cds_feat_atg_pos = $$accepted_cds_feat->end;
	    }

	    # obtain the nearest ppy-ag
	    my $shortest_ppy_feat = $shortest5ppy{$accepted_cds_feat_atg_pos};
	    my $shortest_ppy_pos = $shortest5p{$accepted_cds_feat_atg_pos};

	    # remove all old SLadd tags
	    if($$accepted_cds_feat->has_tag('SLadd')) {
		$$accepted_cds_feat->remove_tag('SLadd');
	    }

	    # polypymodel tags? check with the polya part that also uses this!
	    if($$accepted_cds_feat->has_tag('Polypymodel') ) {
		$$accepted_cds_feat->remove_tag('Polypymodel');
	    }

	    # re-add those for the nearest ppy-ag site
	    $$accepted_cds_feat->add_tag_value('SLadd', $shortest_ppy_pos);  # 'Possible spliced leader addition site at '
	    $$accepted_cds_feat->add_tag_value('Polypymodel',"Agreement: distance from 5p_ag to cds atg is "
					       .abs($accepted_cds_feat_atg_pos - $shortest_ppy_pos));

	    # remove old tags on the now de-referenced ppyag feats? or simply clean or additionaly qualify the "extreme" one?
	    $$shortest_ppy_feat->add_tag_value('SLaddForCDS',"This ppyag was implied as the trans-splicesite for the downstream CDS ".gff3_string($$accepted_cds_feat));
	    $$shortest_ppy_feat->add_tag_value('cds-atg',$accepted_cds_feat_atg_pos);

	}
    }

    if ($long) {
	# longest5p
	foreach my $accepted_cds_feat ( @ofeat ) {

	    # get the cds ATG pos (start or end depending on strand)
	    my $accepted_cds_feat_atg_pos;
	    
	    if( $$accepted_cds_feat->strand == 1 ) {
		$accepted_cds_feat_atg_pos = $$accepted_cds_feat->start;
	    } else {
		$accepted_cds_feat_atg_pos = $$accepted_cds_feat->end;
	    }

	    # obtain the nearest ppy-ag
	    my $longest_ppy_feat = $longest5ppy{$accepted_cds_feat_atg_pos};
	    my $longest_ppy_pos = $longest5p{$accepted_cds_feat_atg_pos};

	    # remove all old SLadd tags
	    if($$accepted_cds_feat->has_tag('SLadd')) {
		$$accepted_cds_feat->remove_tag('SLadd');
	    }

	    # polypymodel tags? check with the polya part that also uses this!
	    if( $$accepted_cds_feat->has_tag('Polypymodel') ) {
		$$accepted_cds_feat->remove_tag('Polypymodel');
	    }

	    # re-add those for the nearest ppy-ag site
	    $$accepted_cds_feat->add_tag_value('SLadd', $longest_ppy_pos); # 'Possible spliced leader addition site at '
	    $$accepted_cds_feat->add_tag_value('Polypymodel',"Agreement: distance from 5p_ag to cds atg is "
					       .abs($accepted_cds_feat_atg_pos - $longest_ppy_pos));
	}
    }
}

sub extremes_clean_cds_3p {	

    if($short) {
        # shortest3p
	foreach my $accepted_cds_feat ( @opafeat ) {

	    # get the cds STOP pos (start or end depending on strand)
	    my $accepted_cds_feat_end_pos;
	    
	    if( $$accepted_cds_feat->strand == 1 ) {
		$accepted_cds_feat_end_pos = $$accepted_cds_feat->end;
	    } else {
		$accepted_cds_feat_end_pos = $$accepted_cds_feat->start;

	    }

	    # obtain the nearest ppy-ag
	    my $shortest_pa_pyag_feat = $shortest3ppy{$accepted_cds_feat_end_pos};
	    my $shortest_polya_pos = $shortest3p{$accepted_cds_feat_end_pos};
	    
	    # remove all old polya tags - Polypymodel tag already removed in 5p step
	    if( $$accepted_cds_feat->has_tag('polya') ) {
		$$accepted_cds_feat->remove_tag('polya');
	    }
	    
           # quick debug
#	    defined($accepted_cds_feat_end_pos) || print STDERR "undefined acc_cds_feat_end for ".$$accepted_cds_feat->gff3_string."\n";
#	    defined($shortest_polya_pos) ||  print STDERR "undefined shortest polyapos for ".$$accepted_cds_feat->gff3_string."\n";

	    # re-add those for the nearest ppy-ag site 
	    $$accepted_cds_feat->add_tag_value('Polypymodel',"Agreement: distance from downstream 5p_ag to rev cds end is "
					       .abs( $accepted_cds_feat_end_pos - $shortest_polya_pos ));
	    $$accepted_cds_feat->add_tag_value('polya', $shortest_polya_pos);
	}
    }

    if ($long) {
        # longest3p 
	foreach my $accepted_cds_feat ( @opafeat ) {

	    # get the cds ATG pos (start or end depending on strand)
	    my $accepted_cds_feat_end_pos;
	    
	    if( $$accepted_cds_feat->strand == 1 ) {
		$accepted_cds_feat_end_pos = $$accepted_cds_feat->end;
	    } else {
		$accepted_cds_feat_end_pos = $$accepted_cds_feat->start;
	    }

	    # obtain the nearest ppy-ag
	    my $longest_pa_pyag_feat = $longest3ppy{$accepted_cds_feat_end_pos};
	    my $longest_polya_pos = $longest3p{$accepted_cds_feat_end_pos};

	    # remove all old polya tags - Polypymodel tag already removed in 5p step
	    if($$accepted_cds_feat->remove_tag('polya')) {
		$$accepted_cds_feat->remove_tag('polya');
	    }

	    # re-add those for the nearest ppy-ag site
	    $$accepted_cds_feat->add_tag_value('Polypymodel',"Agreement: distance from downstream 5p_ag to rev cds end is "
					       .abs( $accepted_cds_feat_end_pos - $longest_polya_pos ));
	    $$accepted_cds_feat->add_tag_value('polya', $longest_polya_pos);
	}
    }
}

sub check_uorf {
    my $selectedppyfeats;

    my $uorffree = 0;
    my $uorfolap = 0;
    my $utrs_checked = 0;    
    my @uaugcount;

    my (@uorfstart, @uorfend);

    if ($long) {
	$selectedppyfeats = [values %longest5ppy];
    } elsif($short) {
	$selectedppyfeats = [values %shortest5ppy];
    } else {
	$selectedppyfeats = \@usedppy;
    }

    foreach my $polypyfeat (@$selectedppyfeats) {
	my (@cdsatgpos) = $$polypyfeat->each_tag_value('cds-atg');
	for (my $i=0; $i < @cdsatgpos; $i++) { # a polypyfeat can be the nearest splice-site for several cds
	    
	    if ($long) {
		$sladdpos = $longest5p{$cdsatgpos[$i]};
	    } elsif($short) {
		$sladdpos = $shortest5p{$cdsatgpos[$i]};
	    } else {
		( ($$polypyfeat->strand == 1) && ($sladdpos = $$polypyfeat->end) ) || ($sladdpos = $$polypyfeat->start);
	    }
	    
	    $DEBUG && debug("Checking 5putr for uorf [".$$polypyfeat->start."-".$$polypyfeat->end."] for direction ".$$polypyfeat->strand.".");
	    
	    my $direction = $$polypyfeat->strand;
	    if($direction == 1) {
		$first = $sladdpos -1;
		$last = $cdsatgpos[$i] - 1;
	    } else {
		$first = $cdsatgpos[$i] + 1;
		$last = $sladdpos +1;
	    }

	    $DEBUG && debug("Getting subseq from  $first to $last.");

	    my $fiveputr = lc $seq->subseq($first,$last);
	    ($direction == -1) && ($fiveputr = revcomp($fiveputr));
	    
	    $utrs_checked++;
	    my $olapped = 0;
	    
	    if($fiveputr =~ m/atg/) {
		my $augs = 0;
		my $aug_position = 0;
		my @augpos;
		
		while( $fiveputr =~ m/(.*?)atg/g ) {
		    $augs++;
		    $aug_position += length($1);
		    push @augpos, $aug_position;
		    
		    $DEBUG && debug("ppyfeat ".$$polypyfeat->start."-".$$polypyfeat->end." (".$$polypyfeat->strand."): uORF start codon found (nr $augs at "
				    .$aug_position." on ".(length($fiveputr))." bp utr).");
		    $aug_position += 3; # for the part of the pattern not included in match..
		}
		$uaugcount[$augs]++;
		
		foreach my $augp ( @augpos ) {
		    my $terminated = 0;
		    my $ds = substr($fiveputr,$augp);
		    
		    # print STDERR "DEBUG: $ds\n";
		  TRI: for(my $i = 0; 3 * $i < length($ds); $i++) {
		      my $tri = substr($ds, 3*$i, 3);
		      if ($tri eq "taa" || $tri eq "tag" || $tri eq "tga") {
			  $endp = $augp + 3*$i;
			  # SOooo, what about "uORFs" ending in the same stop codon? Use a "longest uORF"-approach?
			  # stop codon in frame for this uORF?
			  $DEBUG && debug("uORF at $augp ends at $endp (".($endp-$augp)." bp long)");
			  my ($uORFstart, $uORFend);
			  if($direction == 1) {
			      $uORFstart = $first + $augp;
			      $uORFend = $first + $endp + 2; # -1, then add 3 to include stop codon 
			  } else {
			      $uORFstart = $last - $endp - 2; # +1, then deduct 3 to include stop codon
			      $uORFend = $last - $augp;
			  }

			  # screen for dupes - should actually be done a lot earlier..
			  for (my $olduorfnr=0; $olduorfnr < @uorfstart; $olduorfnr++) {
			      if ($uORFstart == $uorfstart[$olduorfnr] && $uORFend == $uorfend[$olduorfnr]) {
				  $DEBUG && print "Not adding duplicate uORF between $uORFstart and $uORFend.\n";
				  $terminated = 1;
				  last TRI;
			      }
			  }
			  
			  push @uorfstart, $uORFstart;
			  push @uorfend, $uORFend;
			  
			  $seq->add_SeqFeature(new Bio::SeqFeature::Generic ( -start => $uORFstart, -end => $uORFend, -strand => $direction,
									      -primary=>'uORF',
									      -source=>"papya_$version",
									      -tag => { note => "uORF terminated on UTR, length ".($uORFend-$uORFstart) }));
			  $$polypyfeat->add_tag_value('uORF', "start=$augp end=$endp terminated=1 length=".($uORFend-$uORFstart));   # adding uORF tag
			  $terminated = 1;
			  last TRI;
		      }
		  }

		    if( !$terminated && ((length($fiveputr) - $augp) % 3 == 0) ) { # how about 0 % 3 == 0?
			# inframe with actual gene?
			my ($uORFstart, $uORFend);
			if($direction == 1) {   
			    $uORFstart = $first + $augp;
			    $uORFend = $last;
			} else {
			    $uORFstart = $first ;
			    $uORFend = $last - $augp;
			}

			# screen for dupes - should actually be done a lot earlier..
			my $add= 1;
		        for (my $olduorfnr=0; $olduorfnr < @uorfstart; $olduorfnr++) {
			    if ($uORFstart == $uorfstart[$olduorfnr] && $uORFend == $uorfend[$olduorfnr]) {
				$DEBUG && print "Not adding duplicate uORF between $uORFstart and $uORFend.\n";
				$add = 0;
			    }
			}

			if($add) {
			    push @uorfstart, $uORFstart;
			    push @uorfend, $uORFend;

			    $seq->add_SeqFeature(new Bio::SeqFeature::Generic ( -start => $uORFstart, -end => $uORFend, -strand => $direction, 
										-primary=>'uORF',
										-source=>"papya_$version", 
										-tag => { note => "uAUG is in frame with annotated gene, length ".($uORFend-$uORFstart) }));
			    $$polypyfeat->add_tag_value('uORF', "start=$augp terminated=0 overlapping=0 length=".($uORFend-$uORFstart));

			    $DEBUG && debug("uORF at $augp is in frame with annotated gene, and unterminated before annotated start.");
			}
		    } elsif(!$terminated) {
			my ($uORFstart, $uORFend);
			if($direction == 1) {
			    $uORFstart = $first + $augp;
			    $uORFend = $last;
			} else {
			    $uORFstart = $first;
			    $uORFend = $last - $augp;
			}

			# endp really ok???

			if(!defined($endp)) {
			    $endp = $augp; # if distance to olapping actual start < 3bp (one codon).
			}

			# screen for dupes - should actually be done a lot earlier..
			my $add = 1;
		        for (my $olduorfnr=0; $olduorfnr < @uorfstart; $olduorfnr++) {
			    if ($uORFstart == $uorfstart[$olduorfnr] && $uORFend == $uorfend[$olduorfnr]) {
				$DEBUG && print "Not adding duplicate uORF between $uORFstart and $uORFend.\n";
				$add = 0;
			    }
			}

			if($add) {
			    push @uorfstart, $uORFstart;
			    push @uorfend, $uORFend;

			    $seq->add_SeqFeature(new Bio::SeqFeature::Generic ( -start => $uORFstart, -end => $uORFend, -strand => $direction, 
										-primary=>'uORF', 
										-source=>"papya_$version", 
										-tag => { note => "uORF overlaps start of annotated gene." }));
			    $$polypyfeat->add_tag_value('uORF', "start=$augp end=$endp terminated=0 overlapping=1 length=".($uORFend-$uORFstart));
			    $olapped = 1;
			    $DEBUG && debug("uORF at $augp is overlapping the start of the annotated gene.");
			}
		    }
		}
	    } else {
		$$polypyfeat->add_tag_value('uORFfree', length($fiveputr));
		$DEBUG && debug("ppyfeat ".$$polypyfeat->start."-".$$polypyfeat->end." (".$$polypyfeat->strand."): No uORF start codons found (on ".(length($fiveputr))." bp utr).");
		$uorffree++;
	    }

	    $olapped && $uorfolap++;
	}
    }

    return ($uorffree, $uorfolap, $utrs_checked, @uaugcount);
}

sub output_uorfs {
    foreach my $feat ($seq->all_SeqFeatures()) {
	if ($feat->primary_tag eq "uORF") {
	    print $uorfoutfh gff3_string($feat)."\n";
	}
    }
}

sub split_statistics_file {
    my $part = shift;

     if($short) {
	 $part = $part.".short";
     } elsif ($long) {
	 $part = $part.".long";
     }
    if ($rout) {
	$part .= ".r";
    }

    if(defined($statisticsfile)) {
	if(defined($statisticsfh)) {
	    close $statisticsfh;
	}
	$DEBUG && debug("Writing statistics to $statisticsfile.$part.");
	open STATISTICS, ">$statisticsfile.$part" or fatal("Could not open $statisticsfile.$part for writing.");
	$statisticsfh = *STATISTICS;
    }
    return $statisticsfh;
}

sub output_extremes_stats {
    
     # 5p

     my %selectedppyfeats;
     if ($long) {
   	%selectedppyfeats = %longest5ppy;
     } elsif($short) {
    	%selectedppyfeats = %shortest5ppy;
     } else {
# 	foreach my $usedppy ( @usedppy ) {
#	    my (@cdsatgpos) = ($usedppy->each_tag_value('cds-atg'));
# 	    foreach $cdsatgpos ( @cdsatgpos ) {
# 		push @$selectedppyfeat{$cdsatgpos} , $usedppy; # of ppy for that cdsatgpos...
# 	    }
#	}
    }
     if ($splitstats == 1) {
	 $statisticsfh = split_statistics_file("5p");
     }

     if(defined($rout)) {	 
  	print $statisticsfh "pyagstart\tpyagstop\tcdsstart\tcdsstop\tpyagdist\tpolypylen\tpyagatgdist\n";
     }     

     foreach my $accepted_cds_feat ( @ofeat ) {

  	my $strand = $$accepted_cds_feat->strand;
  	my $cdsatgpos = ($strand>0)?$$accepted_cds_feat->start:$$accepted_cds_feat->end;
	
  	foreach my $sladdpos ( $$accepted_cds_feat->each_tag_value('SLadd') ) {

  	    my $polypyfeat = $selectedppyfeats{$cdsatgpos};

  	    my $strandchar = ($strand == 1) ? '+' : '-';

 	    my ($polyagdist) = $$polypyfeat->each_tag_value('pyagdist');
 	    my ($polypylen) = $$polypyfeat->each_tag_value('polypylen');

	    if( $strand > 0 ) {
		my $pyagatgdist = $cdsatgpos - $sladdpos -1 ; # utr-len == (atg-1) - (ag+1) + 1

		if(defined($rout)) {
		    print $statisticsfh $$polypyfeat->start."\t".$$polypyfeat->end."\t".$$accepted_cds_feat->start."\t".$$accepted_cds_feat->end."\t"
			.$polyagdist."\t".$polypylen."\t".$pyagatgdist."\n";
		} else {
		    print $statisticsfh "ppyag\t".$$polypyfeat->start."\t".$$polypyfeat->end."\t".$$accepted_cds_feat->start."\t".$$accepted_cds_feat->end."\t"
			.$polyagdist."\t".$polypylen."\t".$pyagatgdist."\t".gff3_string($$accepted_cds_feat)."\n";
		}
	    } else {
		my $pyagatgdist = ($sladdpos - $cdsatgpos - 1); # utrlen == (ag-1) - (atg+1) + 1 

		if(defined($rout)) {
		    print $statisticsfh $$polypyfeat->start."\t".$$polypyfeat->end."\t".$$accepted_cds_feat->start."\t".$$accepted_cds_feat->end."\t"
			.$polyagdist."\t".$polypylen."\t".$pyagatgdist."\n";
		} else {
		    print $statisticsfh "ppyag\t".$$polypyfeat->start."\t".$$polypyfeat->end."\t".$$accepted_cds_feat->start."\t".$$accepted_cds_feat->end."\t"
			.$polyagdist."\t".$polypylen."\t".$pyagatgdist."\t".gff3_string($$accepted_cds_feat)."\n";
		}		
	    }
	}
    }

     # 3p
     if ($splitstats == 1) {
	 $statisticsfh = split_statistics_file("3p");
     }

     if(defined($rout)) {
	 print $statisticsfh "pyagstart\tpyagstop\tcdsstart\tcdsstop\tpyagdist\tpolypylen\tpana\tpolyapyagdist\ttputrlen\n";
     }
     
     if ($long) {
	 %selectedppyfeats = %longest3ppy;
     } elsif($short) {
	 %selectedppyfeats = %shortest3ppy;
     } else {
     }

     foreach my $accepted_cds_feat ( @opafeat ) {
	 my $strand = $$accepted_cds_feat->strand;

	 # get the cds STOP pos (start or end depending on strand)
	 my $accepted_cds_feat_end_pos;
	    
	 if( $$accepted_cds_feat->strand == 1 ) {
	     $accepted_cds_feat_end_pos = $$accepted_cds_feat->end;
	 } else {
	     $accepted_cds_feat_end_pos = $$accepted_cds_feat->start;
	 }
	 
	 # obtain the nearest ppy-ag
	 
	 if(! exists($selectedppyfeats{$accepted_cds_feat_end_pos})) {
	     fatal("Oops! Apparently ppy feat for accepted cds end pos $accepted_cds_feat_end_pos (strand $strand) is not defined!");
	 }

	 my $polypyfeat = $selectedppyfeats{$accepted_cds_feat_end_pos};
	 
	 my ($cdsendpos) = $$polypyfeat->each_tag_value('cds-end');
	 my ($polyapos) = $$polypyfeat->each_tag_value('polya');
	 my ($pyagdist) = $$polypyfeat->each_tag_value('pyagdist');
	 my ($polypylen) = $$polypyfeat->each_tag_value('polypylen');
	 my ($na) = $$polypyfeat->each_tag_value('polyasite_na');
# FORWARD 3ps
#					my ($tpUTRlen) = ($possible_polyasite - $cfeat[$i]->end - 1); # distance between two inclusive coords..
#					my ($papyagdist) = ($lastpolypyagstart - $possible_polyasite);
# REV 3p 
#					my ($tpUTRlen) = ($cfeat[$i]->start - $possible_polyasite -1 );
#					my ($papyagdist) = ($possible_polyasite - $lastpolypyagstart); # start is the pY end -- the feat end for the rev strand
#	 my ($tpUTRlen) = ($polyapos - $cdsendpos - 1); # distance between two inclusive coords..

	 # strand dep
	 my $papyagdist;
	 my $tpUTRlen;
	 if( $strand > 0 ) {
	     $tpUTRlen = $polyapos - $cdsendpos - 1; # distance between two inclusive coords..
	     $papyagdist = $$polypyfeat->start - $polyapos;
	 } else {
	     $tpUTRlen = $cdsendpos - $polyapos - 1; # cds-end is already start pos for rev strand
	     $papyagdist = $polyapos - $$polypyfeat->end;
	 }

	 if(defined($rout)) {
	     print $statisticsfh $$polypyfeat->start."\t".$$polypyfeat->end."\t".$$accepted_cds_feat->start."\t".$$accepted_cds_feat->end."\t"
		 .$pyagdist."\t".$polypylen."\t".$na."\t".$papyagdist."\t".$tpUTRlen."\n";
	 } else {
	     print $statisticsfh "pappy\t".$$polypyfeat->start."\t".$$polypyfeat->end."\t".$$accepted_cds_feat->start."\t".$$accepted_cds_feat->end
		 ."\t".$polyagdist."\t".$na."\t".$polypylen."\t".$tpUTRlen."\t".gff3_string($accepted_cds_feat)."\n";
	 }
     }
 }


sub polypyag_cds_start_check_genes {

    # and why not make that coditionally nearest only as well.. should make --short runs a LOT faster
    my $seqlen = length($seq->seq());
}

sub polypyag_cds_start_check {

    my $lastpolypyfeat;
    my $seqlen = length($seq->seq());

    $DEBUG && debug("First sort...");
    @cfeat = sort { $a->strand <=> $b->strand || ($a->strand == 1 && $a->start <=> $b->start) || ($a->strand == -1 && $b->end <=> $a->end)} @cfeat;

    $DEBUG && debug("First sort done.");
    # $DEBUG && map { debug($_->gff3_string) } @cfeat;

    if ($splitstats == 1) {    
	$statisticsfh = split_statistics_file("5p");
    }

    if(!($long || $short) && defined($rout)) {
	print $statisticsfh "pyagstart\tpyagstop\tcdsstart\tcdsstop\tpyagdist\tpolypylen\tpyagatgdist\n";
    }

    foreach $direction ( (-1, 1) ) {

        # if user has enabled no-revcomp, don't do -1
	($direction == 1) || $revcomp || next;

	my $lastpolypyagend = (($direction == 1) ? 1 : $seqlen);
	if($cds_5p_ag_atg_dist) {
	    my $featc = 0;
	    my $feats = @cfeat;
FEAT:	    foreach my $feat (@cfeat) {
		$featc++;
		$feat->strand == $direction || next;

		if ($feat->primary_tag eq "pPyAg") {
		    ($lastpolypyagend = (($direction == 1) ? $feat->end : $feat->start));
		    $lastpolypyfeat = \$feat;
		    push @pyfeat, \$feat;
OTHERFEAT:		    for(my $i = $featc; $i < $feats; $i++) {
			if($cfeat[$i]->strand == $direction && $cfeat[$i]->primary_tag eq "CDS") {
			    
			    # be careful about joins - although most of these will be from pseudogenes..
			    if ($cfeat[$i]->location->isa('Bio::Location::SplitLocationI')) {
				# print "DEBUG: Ok, found join..\n";
				$cfeat_start = $cfeat[$i]->location->min_start;
				$cfeat_end = $cfeat[$i]->location->max_end;
			    } else {
				$cfeat_start = $cfeat[$i]->start;
				$cfeat_end = $cfeat[$i]->end;
			    }

			    if($direction == 1 && $cfeat_start - $lastpolypyagend < $cds_5p_ag_atg_dist) {

				# should we markup the particular ppyag that are used with a particular cds, for making a proper statistic?
				if($cfeat_start < $lastpolypyagend) {
				    if($cfeat_end >= $lastpolypyagend ) {
					# that would be indicative of a polypy-ag element spanning an annotated cds start -- does this kill good candidates, further away?
					# yes. so only do a next OTHERFEAT!
					$DEBUG && $WARNING && warning("Possible alternative splice site; direction $direction, lastpolypyagend is $lastpolypyagend, feat-end "
							    .$cfeat_end." feat-start ".$cfeat_start." NOT INCLUDED in continued analysis.");
					if($print_internal_alt_site) { 			    
					  #  PRINT INTERNAL ALTERNATIVE SITE
					    print ALTSITE "direction $direction, lastpolypyagend is $lastpolypyagend, feat-end ".$cfeat_end." feat-start ".$cfeat_start."\n";
					}
					next OTHERFEAT;
				    } else {
					$WARNING && debug("Bug or pyag spanned ORF? featstart ".$cfeat_start."< lastpolypyagend $lastpolypyagend AND featend "
					      .$cfeat_end." < lastpolypyagend $lastpolypyagend"); # assert
					next OTHERFEAT;
				    }
				}

				if($lastpolypyagend != 1) {
                                    # $DEBUG && debug("Agreement: distance 5p_ag to atg is ".($cfeat_start - $lastpolypyagend)." for:\n".$cfeat[$i]->gff3_string);
				    if(!$cfeat[$i]->has_tag('SLadd')) {
					push @ofeat, \$cfeat[$i];
				    } elsif($DEBUG) {
                                        # debug("Feature not added to ok list since its already there.");
				    }

				    # couldn't the stat part be put in a separate function?! Ugly hack, like the rest of this program.
				    if(!($long || $short) && defined($statisticsfile)) { # type, start, stop, pyagdist, polypylen,  pyag-atgdist
					# printing for ppyag implied as transsplicesite for the downstream CDS.
					my ($polyagdist) = $$lastpolypyfeat->each_tag_value('pyagdist');
					my ($polypylen) = $$lastpolypyfeat->each_tag_value('polypylen');
					if(defined($rout)) {
					    print $statisticsfh $$lastpolypyfeat->start."\t".$$lastpolypyfeat->end."\t".$cfeat_start."\t".$cfeat_end."\t"
						.$polyagdist."\t".$polypylen."\t".($cfeat_start - $lastpolypyagend -1)."\n";
					} else {
					    print $statisticsfh "ppyag\t".$$lastpolypyfeat->start."\t".$$lastpolypyfeat->end."\t".$cfeat_start."\t".$cfeat_end."\t"
						.$polyagdist."\t".$polypylen."\t".($cfeat_start - $lastpolypyagend -1)."\t".gff3_string($cfeat[$i])."\n";
					}
				    }

				    $$lastpolypyfeat->add_tag_value('SLaddForCDS',"This ppyag was implied as the trans-splicesite for the downstream CDS ".gff3_string($cfeat[$i]));
				    $$lastpolypyfeat->add_tag_value('cds-atg',$cfeat_start);
				    
				    $cfeat[$i]->add_tag_value('Polypymodel',"Agreement: distance from 5p_ag to cds atg is ".($cfeat_start - $lastpolypyagend -1));
				    $cfeat[$i]->add_tag_value('SLadd', $lastpolypyagend);

				    push @usedppy, $lastpolypyfeat;
				} else {
				    $WARNING && warning("Warning: edge effect detected (direction $direction) in splice model!");
				}
				
				# next FEAT; # next pPyAg - this site will only be marked as a sladd-site for one cds.

			    } elsif ($direction == -1 && $lastpolypyagend - $cfeat_end < $cds_5p_ag_atg_dist) {

				if($lastpolypyagend <= $cfeat_end) { # i e the negative cds_5p_ag_atg_dist
				    if($lastpolypyagend >= $cfeat_start) {

					# that would be indicative of a polypy-ag element spanning an annotated cds start -- does this kill good candidates, further away?
					# yes. so only do a next OTHERFEAT!
					$DEBUG && $WARNING && warning("Possible alternative splice site; direction $direction, lastpolypyagend is $lastpolypyagend, feat-end "
							    .$cfeat_end." feat-start ".$cfeat_start." NOT INCLUDED in continued analysis.");

					# PRINT INTERNAL ALTERNATIVE SITE
					if($print_internal_alt_site) {		
					    # PRINT INTERNAL ALTERNATIVE SITE
					    print ALTSITE "direction $direction, lastpolypyagend is $lastpolypyagend, feat-end ".$cfeat_end." feat-start ".$cfeat_start."\n";
					}

					next OTHERFEAT;
				    } else {
					$WARNING && debug("Bug or pyag spanned ORF? On rev strand, featend ".$cfeat_end." >= lastpolypyagend $lastpolypyagend AND featstart "
							  .$cfeat_start." > lastpolypyagend $lastpolypyagend"); # assert
					next OTHERFEAT;
				    }
				}

				if(!$cfeat[$i]->has_tag('SLadd')) {
				    push @ofeat, \$cfeat[$i];
				} elsif($DEBUG) {
                                   # debug("Feature not added to ok list since its already there.");
				}

				if($lastpolypyagend != $seqlen) {
                                    # $DEBUG && debug("Agreement: distance 5p_ag ($lastpolypyagend) to atg (".$cfeat_end.") is "
				    # .($lastpolypyagend - $cfeat_end)." for:\n".$cfeat[$i]->gff3_string);

				    $cfeat[$i]->add_tag_value('Polypymodel',"Agreement: distance to nearest previous 5p_ag to cds atg is ".($lastpolypyagend - $cfeat_end -1));
				    $$lastpolypyfeat->add_tag_value('cds-atg',$cfeat_end);
				    $$lastpolypyfeat->add_tag_value('SLaddForCDS',"This ppyag was implied as the trans-splicesite for the downstream CDS ".gff3_string($cfeat[$i]));
				    $cfeat[$i]->add_tag_value('SLadd', $lastpolypyagend); # "Possible spliced leader addition site at "
				    
				    if(!($long || $short) && defined($statisticsfile)) { # type, start, stop, pyagdist, polypylen, pyag-atgdist
					# printing for ppyag implied as transsplicesite for the downstream CDS.
					my ($polyagdist) = $$lastpolypyfeat->each_tag_value('pyagdist');
					my ($polypylen) = $$lastpolypyfeat->each_tag_value('polypylen');
					if(defined($rout)) {
					    print $statisticsfh $$lastpolypyfeat->start."\t".$$lastpolypyfeat->end."\t".$cfeat_start."\t".$cfeat_end."\t"
						.$polyagdist."\t".$polypylen."\t".($lastpolypyagend - $cfeat_end -1)."\n";
					} else {
					    print $statisticsfh "ppyag\t".$$lastpolypyfeat->start."\t".$$lastpolypyfeat->end."\t".$cfeat_start."\t".$cfeat_end."\t"
						.$polyagdist."\t".$polypylen."\t".($lastpolypyagend - $cfeat_end -1)."\t".gff3_string($cfeat[$i])."\n";
					}
				    }

				    push @usedppy, $lastpolypyfeat;
				} else {
				    $WARNING && warning("Warning: edge effect detected (direction $direction) in splice model!");
				}

				# next FEAT; # next pPyAg - this site will only be marked as a sladd-site for one cds.
			    } else { 
				# drop to next polypyfeat, since the start (well, end in case of dir..) of this cdsfeat was outside reach of the cds_5p_ag_atg_dist.
				# please, test the result by diffing genes accepted w/ and w/o this
				# ok, at least the gff output is identical!
				next FEAT;
			    }
			}
		    }
		} elsif($feat->primary_tag eq "CDS") {
		    if(! $feat->has_tag('SLadd') ) {
#			$DEBUG && debug("The following CDS is NOT in agreement with splice model (lacks close enough polypy-ag)!\n".$feat->gff3_string);
			push @dfeat, \$feat;
		    }
		}
	    }
	}
    }
}

sub polyasite_check_genes {
   
    # todo: when to call "not accepted cds"

    # split the cfeats into pPyAg sites and cds:es

    my $seqlen = length($seq->seq());
    
    my @cds_feat;
    my @pyag_feat;

    $DEBUG && debug("split ".scalar(@cfeat)." cfeat array into cds and pyag arrays");
    
    for ( my $f = 0; $f < @cfeat ; $f++) {
	if( $cfeat[$f]->primary_tag eq 'CDS' ) {
	    push @cds_feat, \$cfeat[$f];
	} elsif ( $cfeat[$f]->primary_tag eq 'pPyAg' ) {
	    push @pyag_feat, \$cfeat[$f];
	}
    }

    $DEBUG && debug("DEBUG: done splitting into cds ".scalar(@cds_feat)." element array and and pyag ".scalar(@pyag_feat)." arrays.\n");

    # sort cds with smallest end first on forward, largest start first on reverse
    @cds_feat = sort { ${$b}->strand <=> ${$a}->strand || (${$a}->strand == 1 && ${$a}->end <=> ${$b}->end) || (${$a}->strand == -1 && ${$b}->start <=> ${$a}->start)} @cds_feat;
    # sort pyag with smallest start first on forward, largest end first on reverse
    @pyag_feat = sort { ${$b}->strand <=> ${$a}->strand || (${$a}->strand == 1 && ${$a}->start <=> ${$b}->start) || (${$a}->strand == -1 && ${$b}->end <=> ${$a}->end)} @pyag_feat;

    $DEBUG && debug("\nSecond sort done.\n");

#    for (my $i = 0; $i < @pyag_feat; $i++ ) {
#	print STDERR ${$pyag_feat[$i]}->start." ".${$pyag_feat[$i]}->end." ".${$pyag_feat[$i]}->strand."\n";
#    }

#    for (my $i = 0; $i < @cds_feat; $i++ ) {
#	print STDERR ${$cds_feat[$i]}->start." ".${$cds_feat[$i]}->end." ".${$cds_feat[$i]}->strand."\n";
#    }

   
    if ($splitstats == 1) {
	$statisticsfh = split_statistics_file("3p");
    }

    if(!($long || $short) && defined($rout)) {
	print $statisticsfh "pyagstart\tpyagstop\tcdsstart\tcdsstop\tpyagdist\tpolypylen\tpana\tpolyapyagdist\ttputrlen\n";
    }
    
    my $last_used_pyag_index = 0;
    for( my $i = 0; $i < @cds_feat; $i++ ) {
	
	# remember to care for joins (eg introns)

	if (${$cds_feat[$i]}->location->isa('Bio::Location::SplitLocationI')) {
	    $cds_feat_start = ${$cds_feat[$i]}->location->min_start;
	    $cds_feat_end = ${$cds_feat[$i]}->location->max_end;
	} else {
	    $cds_feat_start = ${$cds_feat[$i]}->start;
	    $cds_feat_end = ${$cds_feat[$i]}->end;
	}

        # forward

	if(${$cds_feat[$i]}->strand == 1) {

	    # if last on this strand..
	    if($i == @cds_feat-1  || ${$cds_feat[$i+1]}->strand == -1) {
		$next_cds_start = $seqlen;
	    } else {
		if (${$cds_feat[$i+1]}->location->isa('Bio::Location::SplitLocationI')) {
		    $next_cds_start = ${$cds_feat[$i+1]}->location->min_start;
		} else {
		    $next_cds_start = ${$cds_feat[$i+1]}->start;
		}
	    }
	    
	    for($j = $last_used_pyag_index; $j < @pyag_feat; $j++) {
		if(${$pyag_feat[$j]}->strand == -1) {
		    # done all feats on this strand now!
		    $DEBUG && debug("No more pyag feats on forward strand!");
		    last;
		}
		my $pyag_feat_start = ${$pyag_feat[$j]}->start;

		if($downstream_gene_interferes_with_splicing && !${$pyag_feat[$j]}->has_tag('SLaddForCDS')) {
		    # only look at signals implicated as actual sladd-sites

		    if($pyag_feat_start > $next_cds_start) {
			# site is downstream of next gene start. begin looking at next gene.
			$DEBUG && debug("Site is *not* SLadd for any gene, but downstream of next gene start codon ($pyag_feat_start > $next_cds_start). "
					."Begin looking at next gene.");
			last;
		    } else {
			$DEBUG && debug("Not SLadd (start $pyag_feat_start) - next!");
			next;
		    }
		}
		
		if( $pyag_feat_start - $cds_feat_end > $exp_stop_polya_ag_distance ) {
		    # outside accepted 3p utr length - begin checking next gene!
		    $DEBUG && debug("cds end $cds_feat_end and pyag start $pyag_feat_start outside accepted 3p utr length - begin checking next gene!");
		    last;
		}

		if($pyag_feat_start < $cds_feat_end) {
		    $DEBUG && debug("Ignoring cds-internal or cds-preceeding py-ag (start "
				    .$pyag_feat_start."< cdsend $cds_feat_end) in polyasite scan.");

                    # its now safe to move the leftmost pyag scan start index forward to this position 
		    # ia not examine this and the previous sites for the next cds
		    # unless ofcourse this cds is completely contained in another ORF...
		    
		    $last_used_pyag_index = $j+1;

		    next;
		}
		
		my ($possible_polyasite, $na) = (-1,0);
		if ( $aaa_scan ) {
		    ($possible_polyasite, $na) = get_aaa_around ( $pyag_feat_start - $polya_pyag_dist, 1, $cds_feat_end, $pyag_feat_start );
		} elsif ( $closest_a ) {
		    ($possible_polyasite, $na) = get_closest_a_around_center ( $pyag_feat_start - $polya_pyag_dist, 1, $cfeat_end, $pyag_feat_start );
		}

		if($possible_polyasite < $cds_feat_end) {
		    # polya site predicted within gene 
		    $DEBUG && debug("Ignoring cds-internal predicted polyasite ($possible_polyasite, pyag start $pyag_feat_start < cdsend $cds_feat_end) in polyasite scan.");
		    next;
		}
	    
		# assert
		if($possible_polyasite > $pyag_feat_start) {
		    $DEBUG && fatal("Predicted pyag feat internal poly adenylation site - this shouldn't happen! Check prediction routine!");
		    exit(-1);
		}

 		if($possible_polyasite > $next_cds_start && $downstream_gene_interferes_with_splicing) {
 		    $DEBUG && debug("Predicted polyadenylation site is past start of ds gene - ignoring this, marking last as the ultimate polyasite and skipping to next gene.");
 		    # note that this check can be too late if downstream genes overlap so that the $i+1 gene is contained completely in a further downstream ending one 
		    # only using this to avoid "unneeded" computation.. wonder if it introduces more than it saves?
 		    last;
 		}
		
		# all checks passed ok - this is a putative polyadenylation site for this cds!

		if(!${$cds_feat[$i]}->has_tag('polya')) { # the first time a certain cds is connected with a polyasite, a ref to it is added to the ok-polya-feat list
		    push @opafeat, $cds_feat[$i];
		}

		$DEBUG && debug("This ppyag ($pyag_feat_start) was implied as a possible polyadenylation position ($possible_polyasite) determining site for the upstream CDS:\n"
				.gff3_string(${$cds_feat[$i]})."\n");

		${$pyag_feat[$j]}->add_tag_value('paForCDS',"This ppyag was implied as a possible polyadenylation position determining site for the upstream CDS "
					       .gff3_string(${$cds_feat[$i]}));
		${$pyag_feat[$j]}->add_tag_value('polya', $possible_polyasite);
		${$pyag_feat[$j]}->add_tag_value('polyasite_na', $na);
		${$pyag_feat[$j]}->add_tag_value('cds-end',$cds_feat_end);
		# distance between two inclusive coordinates
		${$cds_feat[$i]}->add_tag_value('Polypymodel',"Agreement: distance from downstream 5p_ag to cds end is ".($pyag_feat_start - $cds_feat_end));
		${$cds_feat[$i]}->add_tag_value('polya', $possible_polyasite);
	   	    
		if(!($long || $short) && defined($statisticsfile)) { # type, start, stop, pyagdist, polypylen, polya-pyag
		    # printing for ppyag implied as determinant for polya-site for the upstream CDS
		    my ($pyagdist) = ${$pyag_feat[$j]}->each_tag_value('pyagdist');
		    my ($polypylen) = ${$pyag_feat[$j]}->each_tag_value('polypylen');
		    my ($tpUTRlen) = ($possible_polyasite - $cds_feat_end - 1); # distance between two inclusive coords..
		    my ($papyagdist) = ($pyag_feat_start - $possible_polyasite);
		    
		    if(defined($rout)) {
			print $statisticsfh ${$pyag_feat[$j]}->start."\t".${$pyag_feat[$j]}->end."\t".$cds_feat_start."\t".$cds_feat_end."\t".$pyagdist
			    ."\t".$polypylen."\t".$na."\t".$papyagdist."\t".$tpUTRlen."\n";
		    } else {
			print $statisticsfh "pappy\t".${$pyag_feat[$j]}->start."\t".${$pyag_feat[$j]}->end."\t".$cds_feat_start."\t".$cds_feat_end."\t".$polyagdist
			    ."\t".$polypylen."\t".$tpUTRlen."\t".gff3_string(${$cds_feat[$i]})."\n";
		    }
		}

		push @usedpappy, $pyag_feat[$j]; # ref
	    }

	} else {

	    # reverse

	    if( $i == @cds_feat-1 ) {
		$next_cds_end = 1;
	    } else {
		if (${$cds_feat[$i+1]}->location->isa('Bio::Location::SplitLocationI')) {
		    $next_cds_end = ${$cds_feat[$i+1]}->location->max_end; # more like the cds start.. 
		} else {
		    $next_cds_end = ${$cds_feat[$i+1]}->end;
		}
	    }       	

	    for($j = $last_used_pyag_index; $j < @pyag_feat; $j++) {
		my $pyag_feat_end = ${$pyag_feat[$j]}->end; # well, really more like the pyag start, but with it being the numerically higher of the coordinate pair members, I'll keep the name "end"     
		if( ${$pyag_feat[$j]}->strand == 1) {
		    # still on the forw strand in pyag sites.. bomb ahead!
		    $last_used_pyag_index = $j+1;
		    next;
		}

		if($downstream_gene_interferes_with_splicing && !${$pyag_feat[$j]}->has_tag('SLaddForCDS')) {
		    # only look at signals implicated as actual sladd-sites
		    if($pyag_feat_end < $next_cds_end) {
			# site is sladd for some gene, but downstream of next gene start codon. begin looking at next gene.
			$DEBUG && debug("Site is SLadd for some gene, but downstream of next gene start codon ($pyag_feat_end < $next_cds_end). begin looking at next gene.");
			last;
		    } else {
			$DEBUG && debug("Not SLadd (end $pyag_feat_end) - next!");
			next;
		    }
		}

		if( $cds_feat_start - $pyag_feat_end > $exp_stop_polya_ag_distance ) {
		    # outside accepted 3p utr length - begin checking next gene!
		    $DEBUG && debug("cds start $cds_feat_start and pyag end $pyag_feat_end outside accepted 3p utr length - begin checking next gene!");
		    last;
		}

		if($pyag_feat_end > $cds_feat_start) {
		    $DEBUG && debug("Ignoring cds-internal or cds-preceeding py-ag (end ".$pyag_feat_end." > cdsstart $cds_feat_start) in rev strand polyasite scan.");
                    # its now safe to move the pyag scan start index forward to this position 
		    # ia not examine this and the previous sites for the next cds
		    # unless ofcourse this cds is completely contained in another ORF...		    
		    $last_used_pyag_index = $j+1;
		    next;
		}

		my ($possible_polyasite, $na) = (-1, 0);
		if( $aaa_scan ) {
		    ($possible_polyasite, $na) = get_aaa_around ( $pyag_feat_end + $polya_pyag_dist, ${$cds_feat[$i]}->strand, $cds_feat_start, $pyag_feat_end );
		} elsif ( $closest_a )  { 
		    ($possible_polyasite, $na) = get_closest_a_around_center ( $pyag_feat_end + $polya_pyag_dist, ${$cds_feat[$i]}->strand, $cds_feat_start, $pyag_feat_end );
		}

		if($possible_polyasite > $cds_feat_start) {
		    # polya site predicted within gene 
		    $DEBUG && debug("Ignoring cds-internal predicted polyasite ($possible_polyasite, pyag end $pyag_feat_end, cdsstart $cds_feat_start) in rev strand polyasite scan.");
		    next;
		}

		# assert
		if($possible_polyasite < $pyag_feat_end) {
		    $DEBUG && fatal("Predicted pyag feat internal poly adenylation site - this shouldn't happen! Check prediction routine!");
		    exit(-1);
		}

 		if($downstream_gene_interferes_with_splicing && $possible_polyasite < $next_cds_end ) {
 		    $DEBUG && debug("Predicted polyadenylation site is past start of ds gene - ignoring this, marking last as the ultimate polyasite and skipping to next gene.");
 		    # note that this check can be too late if downstream genes overlap so that the $i+1 gene is contained completely in a further downstream ending one 
		    # only using this to avoid "unneeded" computation.. wonder if it introduces more than it saves?
 		    last;
 		}

		# all checks passed ok - this is a putative polyadenylation site for this cds!

		$DEBUG && debug("This ppyag ($pyag_feat_end) was implied as a possible polyadenylation position ($possible_polyasite) determining site for the upstream CDS:\n"
				.gff3_string(${$cds_feat[$i]})."\n");

		if(!${$cds_feat[$i]}->has_tag('polya')) { # the first time a certain cds is connected with a polyasite, a ref to it is added to the ok-polya-feat list 
		    push @opafeat, $cds_feat[$i];
		}

		${$pyag_feat[$j]}->add_tag_value('paForCDS',"This ppyag was implied as a potential polyadenylation position directing site for the upstream CDS ".
						gff3_string(${$cds_feat[$i]}));
		${$pyag_feat[$j]}->add_tag_value('cds-end',$cds_feat_start); # DEBUG: test moving 
		${$pyag_feat[$j]}->add_tag_value('polya', $possible_polyasite);
		${$pyag_feat[$j]}->add_tag_value('polyasite_na', $na);
		${$cds_feat[$i]}->add_tag_value('Polypymodel',"Agreement: distance from downstream 5p_ag to rev cds end is ".( $cds_feat_start - $pyag_feat_end ));
		${$cds_feat[$i]}->add_tag_value('polya', $possible_polyasite);
		
		if(!($long || $short) && defined($statisticsfile)) { # type, start, stop, pyagdist, polypylen, polya-pyag
		    # printing for ppyag implied as determinant for polya-site for the upstream CDS
		    my ($pyagdist) = ${$pyag_feat[$j]}->each_tag_value('pyagdist');
		    my ($polypylen) = ${$pyag_feat[$j]}->each_tag_value('polypylen');
		    my ($tpUTRlen) = ($cds_feat_start - $possible_polyasite -1 );
		    my ($papyagdist) = ($possible_polyasite - $pyag_feat_end); # start is the pY end -- the feat end for the rev strand
		    
		    if(defined($rout)) {
			print $statisticsfh ${$pyag_feat[$j]}->start."\t".${$pyag_feat[$j]}->end."\t".$cds_feat_start."\t".$cds_feat_end."\t".$pyagdist."\t"
			    .$polypylen."\t".$na."\t".$papyagdist."\t".$tpUTRlen."\n";
		    } else {
			print $statisticsfh "pappy\t".${$pyag_feat[$j]}->start."\t".${$pyag_feat[$j]}->end."\t".$cds_feat_start."\t".$cds_feat_end."\t".$polyagdist
			    ."\t".$polypylen."\t".$tpUTRlen."\t".gff3_string(${$cds_feat[$i]})."\n";
		    }
		}
		
		push @usedpappy, $pyag_feat[$j]; # ref..
	    }
	}

	# we have now been through all possible pyag feats for this cds (in the case of "short", we wouldn't have to check more than one, if it weren't for the fact
	# that the "extremes" filter is run after this routine..
	
	if(! ${$cds_feat[$i]}->has_tag('polya') ) {
	    $DEBUG && debug("The following CDS is NOT in agreement with polya-model (lacks close enough polypy-ag; may have a huge 3'-UTR?)!\n".gff3_string(${$cds_feat[$i]}));
	    push @dpafeat, $cds_feat[$i];
	}
    }
}

# sub polyasite_check {

#     my $lastpolypyfeat;
#     my $seqlen = length($seq->seq());

#     # resort: costly, but it saves quite some extra errorprone logic in the loop..
#     @cfeat = sort { $a->strand <=> $b->strand || ($a->strand == 1 && $b->end <=> $a->end) || ($a->strand == -1 && $a->start <=> $b->start)} @cfeat;

#     $DEBUG && debug("\nSecond sort done.\n");
# #    $DEBUG && map { debug($_->gff3_string) } @cfeat;

#     if ($splitstats == 1) {
# 	$statisticsfh = split_statistics_file("3p");
#     }

#     if(!($long || $short) && defined($rout)) {
# 	print $statisticsfh "pyagstart\tpyagstop\tcdsstart\tcdsstop\tpyagdist\tpolypylen\tpana\tpolyapyagdist\ttputrlen\n";
#     }

#     # polya-site feasibility model

#     foreach $direction ( (-1, 1) ) {
# 	($direction == 1) || $revcomp || next;
# 	if($exp_stop_polya_ag_distance) {
#             # ehr, this is not really used at the moment, right?!
# 	    $lastpolypyagend = (($direction == 1) ? $seqlen : 1); # ..but remember to check those edge-guys.. 
# 	    # keep a feature counter to be able to "scan" for next ds cds from a polypy 
# 	    my $featc = -1;
# 	    my $feats = @cfeat;
# FEAT:	    foreach my $feat (@cfeat) {
# 		$featc++;
# 		$feat->strand == $direction || next FEAT;
# 		if ($feat->primary_tag eq "pPyAg") { # should really be represented as a join btw the polY and the AG..
# 		    $lastpolypyfeat = \$feat;
# 		    $lastpolypyagend = ($direction == 1) ? $feat->end : $feat->start;
# 		    $lastpolypyagstart = ($direction == 1) ? $feat->start : $feat->end;

# OTHERFEAT:	    for(my $i = $featc; $i < $feats; $i++) {

                            
#                             if($cfeat[$i]->primary_tag eq 'CDS') {
# 			      # joined feature (as in gene with intron) 
# 			      if ($cfeat[$i]->location->isa('Bio::Location::SplitLocationI')) {
# 				  $cfeat_start = $cfeat[$i]->location->min_start;
# 				  $cfeat_end = $cfeat[$i]->location->max_end;
# 			      } else {
# 				  $cfeat_start = $cfeat[$i]->start;
# 				  $cfeat_end = $cfeat[$i]->end;
# 			      }
# 			  }

# 			if($cfeat[$i]->strand == $direction && $cfeat[$i]->primary_tag eq 'CDS') {

# 			    if($direction == 1 && $lastpolypyagend - $cfeat_end  < $exp_stop_polya_ag_distance) { 
# 				$lastpolypyagend >= $cfeat_end or debug("Warning: possible bug OR cds internal polya-signal with valid ds site: featend ".$cfeat_end
# 									     ."> lastpolypyagend $lastpolypyagend in polya model direction $direction.");
				
# 				    # $DEBUG && debug("Disagree: The predicted 3'UTR is shorter (".( $lastpolypyagend - $cfeat_end )."bp) than the polya-pyag-distance ("
# 				    # .$polya_pyag_dist.") leaving a negatively sized 3'UTR..");
# 				# ..but that certainly poses a scanning problem, since the scan is designed to pick up the last polypy, right...		
# #				    next OTHERFEAT;
# #				}

# 				if( $lastpolypyagstart < $cfeat_end ) {
# 				    next OTHERFEAT;
# 				}

# 				if($lastpolypyagend != $seqlen) {
#                                     # $DEBUG && debug("Agreement: distance cds-stop to 5p_ag is ".($lastpolypyagend - $cfeat_end)." for:\n".$cfeat[$i]->gff3_string);
# 				    # no - lets measuse from the 5p of the pPyAg rather than the polypyagend!

# 				    # attempt to locate a suitable polya-site
# 				    my ($possible_polyasite, $na) = (-1,0);
# 				    if ( $aaa_scan ) {
# 					($possible_polyasite, $na) = get_aaa_around ( $lastpolypyagstart - $polya_pyag_dist, $cfeat[$i]->strand, $cfeat_end, $lastpolypyagstart);
# 				    } elsif ( $closest_a ) {
# 					($possible_polyasite, $na) = get_closest_a_around_center ( $lastpolypyagstart - $polya_pyag_dist, $cfeat[$i]->strand, $cfeat_end, $lastpolypyagstart);
# 				    }

# 				    if( $na == 0) {
# 					# score one "against" the CC aaa hypothesis
# #					$possible_polyasite = $lastpolypyagstart - $polya_pyag_dist;

# 					if( ($possible_polyasite < $cfeat_end)
# 					    || ($possible_polyasite > $lastpolypyagstart)) {
# 					    next OTHERFEAT;
# 					}
					
# 					# print STDERR "The AAA hyp was not useful for cds ".$cfeat_start." ".$cfeat_end." ".$cfeat[$i]->strand." ppyagend $lastpolypyagend.\n";
#                                         # if ( $lastpolypyagstart - $cfeat_end < $polya_pyag_dist ) {
# 				    } 
#                                     # else {
#                                     #   print STDERR "Using AAA to set polya to ".$possible_polyasite." for cds ".$cfeat_start." "
#                                     #     .$cfeat_end." ".$cfeat[$i]->strand." ppyagend $lastpolypyagend.\n";
#                                     # }
				    
# 				    if( ($possible_polyasite < $cfeat_end)
# 					|| ($possible_polyasite > $lastpolypyagstart)) {
# 					$WARNING && debug("Possible internal or out of bounds polya-site discarded with aaa check routine ($possible_polyasite, cds end on forw strand "
# 							  .($cfeat_start).", lastpolypyagstart $lastpolypyagstart)" );

# 					next OTHERFEAT;
# 				    }

# 				    if(!$cfeat[$i]->has_tag('polya')) { # the first time a certain cds is connected with a polyasite, a ref to it is added to the ok-polya-feat list
# 					push @opafeat, \$cfeat[$i];
# 				    } 

# 				    $$lastpolypyfeat->add_tag_value('paForCDS',"This ppyag was implied as a possible polyadenylation position determining site for the upstream CDS "
# 								    .$cfeat[$i]->gff3_string);
# 				    $$lastpolypyfeat->add_tag_value('polya', $possible_polyasite);
# 				    $$lastpolypyfeat->add_tag_value('polyasite_na', $na);
# 				    $$lastpolypyfeat->add_tag_value('cds-end',$cfeat_end);
# 				    $cfeat[$i]->add_tag_value('Polypymodel',"Agreement: distance from downstream 5p_ag to cds end is ".($lastpolypyagstart - $cfeat_end));
# 				    $cfeat[$i]->add_tag_value('polya', $possible_polyasite);

# 				    if(!($long || $short) && defined($statisticsfile)) { # type, start, stop, pyagdist, polypylen, polya-pyag
# 					# printing for ppyag implied as determinant for polya-site for the upstream CDS
# 					my ($pyagdist) = $$lastpolypyfeat->each_tag_value('pyagdist');
# 					my ($polypylen) = $$lastpolypyfeat->each_tag_value('polypylen');
# 					my ($tpUTRlen) = ($possible_polyasite - $cfeat_end - 1); # distance between two inclusive coords..
# 					my ($papyagdist) = ($lastpolypyagstart - $possible_polyasite);
					
# 					if(defined($rout)) {
# 					    print $statisticsfh $$lastpolypyfeat->start."\t".$$lastpolypyfeat->end."\t".$cfeat_start."\t".$cfeat_end."\t".$pyagdist
# 						."\t".$polypylen."\t".$na."\t".$papyagdist."\t".$tpUTRlen."\n";
# 					} else {
# 					    print $statisticsfh "pappy\t".$$lastpolypyfeat->start."\t".$$lastpolypyfeat->end."\t".$cfeat_start."\t".$cfeat_end."\t".$polyagdist
# 						."\t".$polypylen."\t".$tpUTRlen."\t".$cfeat[$i]->gff3_string."\n";
# 					}
# 				    }

# 				    push @usedpappy, $lastpolypyfeat;
# 				} else {
# 				    $WARNING && warning("Warning: edge effect detected (direction $direction) in polya model!");
# 				}

# #				next FEAT; # next pPyAg - this site will only be marked as a polyasite for one cds.

# 			    } elsif($direction == -1 && $cfeat_start - $lastpolypyagend < $exp_stop_polya_ag_distance) {

# 				$lastpolypyagend <= $cfeat_start or debug("Warning: possible bug OR cds internal polya-signal with valid ds site:  featstart "
# 									  .$cfeat_start."< lastpolypyagend $lastpolypyagend in polya model direction $direction.");

# 				if($lastpolypyagend != 1) {
#                                     # $DEBUG && debug("Agreement: distance rev cds-stop to 5p_ag is ".( $cfeat_start - $lastpolypyagend )." for:\n".$cfeat[$i]->gff3_string);
				    
# 				    if( $lastpolypyagstart > $cfeat_start ) {
# 					next OTHERFEAT;
# 				    }

				    
# 				    my ($possible_polyasite, $na) = (-1, 0);
# 				    if( $aaa_scan ) {
# 					($possible_polyasite, $na) = get_aaa_around ( $lastpolypyagstart + $polya_pyag_dist, $cfeat[$i]->strand, $cfeat_start, $lastpolypyagstart );
# 				    } elsif ( $closest_a )  { 
# 					($possible_polyasite, $na) = get_closest_a_around_center ( $lastpolypyagstart + $polya_pyag_dist, $cfeat[$i]->strand, $cfeat_start, $lastpolypyagstart );
# 				    }
# 				    if( $possible_polyasite == -1 ) {
# 					# score one "against" the CC aaa hypothesis
# 					$possible_polyasite = $lastpolypyagstart + $polya_pyag_dist;
# 					# print STDERR "The AAA hyp was not useful for cds ".$cfeat_start." ".$cfeat_end
# 					# ." ".$cfeat[$i]->strand." ppyagend $lastpolypyagend.\n";
# 					# if( $possible_polyasite <  $cfeat_start ) {
# #					if ( $cfeat_start - $lastpolypyagstart < $polya_pyag_dist ) {					    
# 					    # $DEBUG && debug("Disagree: The predicted 3'UTR is shorter (".( $cfeat_start - $lastpolypyagend )."bp) than the polya-pyag-distance ("
# 					    # .$polya_pyag_dist."), leaving a negatively sized 3'UTR for ".$cfeat[$i]->gff3_string);
# 					#    next OTHERFEAT;
# 					# }
# 					if( ($possible_polyasite > $cfeat_start)
# 					    || ($possible_polyasite < $lastpolypyagstart) ) { 
# 					    # $WARNING && debug("Possible internal polya-site discarded ($possible_polyasite, cds end on rev strand ".($cfeat_start).")" );
# 					    next OTHERFEAT;
# 					}
# 				    } else {
#                                         # print STDERR "Using AAA to set polya to ".$possible_polyasite." for cds ".$cfeat_start." "
# 					# .$cfeat_end." ".$cfeat[$i]->strand." ppyagend $lastpolypyagend.\n";
# 				    }

# 				    if( ($possible_polyasite > $cfeat_start)
# 					|| ($possible_polyasite < $lastpolypyagstart) ) { 
# 					$WARNING && debug("Possible internal or out of bounds polya-site discarded with aaa check routine ($possible_polyasite, cds end on rev strand "
# 							  .($cfeat_start).", lastpolypyagstart $lastpolypyagstart)");
# 					next OTHERFEAT;
# 				    }

# 				    if(!$cfeat[$i]->has_tag('polya')) { # the first time a certain cds is connected with a polyasite, a ref to it is added to the ok-polya-feat list 
# 					push @opafeat, \$cfeat[$i];
# 				    }

# 				    $$lastpolypyfeat->add_tag_value('paForCDS',"This ppyag was implied as a potential polyadenylation position directing site for the upstream CDS ".
# 								    $cfeat[$i]->gff3_string);
# 				    $$lastpolypyfeat->add_tag_value('cds-end',$cfeat_start); # DEBUG: test moving 
# 				    $$lastpolypyfeat->add_tag_value('polya', $possible_polyasite);
# 				    $$lastpolypyfeat->add_tag_value('polyasite_na', $na);
# 				    $cfeat[$i]->add_tag_value('Polypymodel',"Agreement: distance from downstream 5p_ag to rev cds end is ".( $cfeat_start - $lastpolypyagend ));
# 				    $cfeat[$i]->add_tag_value('polya', $possible_polyasite);

# 				    if(!($long || $short) && defined($statisticsfile)) { # type, start, stop, pyagdist, polypylen, polya-pyag
# 					# printing for ppyag implied as determinant for polya-site for the upstream CDS
# 					my ($pyagdist) = $$lastpolypyfeat->each_tag_value('pyagdist');
# 					my ($polypylen) = $$lastpolypyfeat->each_tag_value('polypylen');
# 					my ($tpUTRlen) = ($cfeat_start - $possible_polyasite -1 );
# 					my ($papyagdist) = ($possible_polyasite - $lastpolypyagstart); # start is the pY end -- the feat end for the rev strand

# 					if(defined($rout)) {
# 					    print $statisticsfh $$lastpolypyfeat->start."\t".$$lastpolypyfeat->end."\t".$cfeat_start."\t".$cfeat_end."\t".$pyagdist."\t"
# 						.$polypylen."\t".$na."\t".$papyagdist."\t".$tpUTRlen."\n";
# 					} else {
# 					    print $statisticsfh "pappy\t".$$lastpolypyfeat->start."\t".$$lastpolypyfeat->end."\t".$cfeat_start."\t".$cfeat_end."\t".$polyagdist
# 						."\t".$polypylen."\t".$tpUTRlen."\t".$cfeat[$i]->gff3_string."\n";
# 					}
# 				    }
				    
# 				    push @usedpappy, $lastpolypyfeat;
# 				} else {
# 				    $WARNING && warning("Warning: edge effect detected (direction $direction) in polya model!");
# 				}
				
# 				# next FEAT; # next pPyAg - this site will only be marked as a polyasite for one cds.
# 			    } else {
# 				next FEAT; # next pPyAg - we should now be out of reach for any possible matches -- major speedup for long sequences!
# 			    }
# 			}
# 		    }
# 		} elsif($feat->primary_tag eq 'CDS') {
# 		    if(! $feat->has_tag('polya') ) {
# #			$DEBUG && debug("The following CDS is NOT in agreement with polya-model (lacks close enough polypy-ag; may have a huge 3'-UTR?)!\n".$feat->gff3_string);

# 			push @dpafeat, \$feat;
# 		    }
# 		}
# 	    }
# 	}
#     }
# }

sub get_closest_a_around_center {

    if( $no_aaa_check == 1) {
	return (-1, 0);
    }

    # note that the first or last ones may not be ok.. What about strandswitches?
    
    my $approximate_polya_site = shift;
    my $strand = shift;
    my $cds_end = shift;
    my $ppystart = shift;

    my $seqlen = length($seq->seq());
    
    # safety checks on the start and end of subseq window

    my $first_winstart = $approximate_polya_site - $polya_aaa_scanning_window >= 0 ? $approximate_polya_site - $polya_aaa_scanning_window : 0;
    my $first_winend = $approximate_polya_site;
    my $second_winstart = $approximate_polya_site+1;
    my $second_winend = $approximate_polya_site + $polya_aaa_scanning_window  < $seqlen ? $approximate_polya_site + $polya_aaa_scanning_window : $seqlen; 
    
    if ( $second_winend - $first_winstart +1 < 2 * $polya_aaa_scanning_window +1) {
	$WARNING && debug("WARNING: polya site scanning window smaller (".( $second_winend - $first_winstart +1 )
			  .") than given winsize (".(2 * $polya_aaa_scanning_window +1).") presumably due to contig break.");
    }

    # the totalwin center base is the last base of first win SO FAR

    my $first_window = lc $seq->subseq($first_winstart, $first_winend);
    my $second_window = lc $seq->subseq($second_winstart, $second_winend);

    # 1-->|2-->
    # <--2|<--1

    # swap and revcomp windows if reverse strand
    if( $strand < 0 ) {
	my $tmp = revcomp( $first_window );
	$first_window = revcomp( $second_window );
	$second_window = $tmp;
	# NOW in reverse, the totalwin center base is the first base of second win 

	# the start coordinates remain the lower coords for both windows, but swapped along with the sequence
	$tmp = $first_winstart;
	$first_winstart = $second_winstart;
	$second_winstart = $tmp;

	$tmp = $first_winend;
	$first_winend = $second_winend;
	$second_winend = $tmp;
    }

    my $polyasite = $approximate_polya_site;
    
    my $na = 0;
    my $dist = 0;
    my $first_polyasite = 0;

    # fist win 
    my @win = split / */, $first_window;

#    $DEBUG && debug("first win length is ".scalar(@win)." and end - start is ".($first_winend - $first_winstart));

    # $i is a coordinate (1 based, inclusive)
    for (my $i = @win; $i > 0; $i--) {

	my $base = $win[$i-1];
      
	if($base ne 'a') {
	    next;
	}

	# 1-->|2-->
	# <--2|<--1
	# now examining the "1" window, not the LEFTMOST window
	
	$na = 1;
	if( $strand > 0 ) {
	    $first_polyasite = $first_winstart - 1 + $i ; # winstart is the first base, and is also included on the interval. count $i bases from the base before the win.
	} else {
	    $first_polyasite = $first_winend + 1 - $i;
	}
	$dist = @win - $i; # match to the window center (first base checked in this win, which is checked in the rev direction) has $dist == 0
	($strand > 0) && $dist++;  # unless this is the reverse strand, in which case the first base checked is at a distance (rev dir of win, remember) 1 

	last;
    }
    
    if( $na == 0 ) {
	$first_polyasite = $approximate_polya_site;
	$dist = @win;
    }

    # second win
    @win = split / */, $second_window;
    my $second_dist = 0;
    my $second_na = 0;
    my $second_polyasite = 0;
    
#    print STDERR "second win length is ".scalar(@win)." and end - start is ".($second_winend - $second_winstart)."\n";
    for (my $i = 0; $i < @win ; $i++) {
      
	my $base = $win[$i];

	if($base ne 'a') {
	    next;
	}

	# 1-->|2-->
	# <--2|<--1
	# now examining the "2" subwin, not the RIGHTMOST
      
	$second_na = 1;
	if ( $strand > 0 ) {
	    $second_polyasite = $second_winstart + $i;
	} else {
	    $second_polyasite = $second_winend - $i;
	}
	$second_dist = $i; # match to the first win base has dist 1 from win-center (which is found in subwin 1 for forward bases)
	($strand < 0) && $dist++; # unless this is on the reverse, where the wincenter is first base in subwin 2..
	last;
    }

    if( $second_na == 0 ) {
	$second_polyasite = $approximate_polya_site;
        $second_dist = @win + 1;           # not found in this win
    }

    # closest of first and second?
    if($na == $second_na) { # 0 or 1
#	print STDERR "na is $na and $second_na\n";

	if ( $second_na == 1 && $second_dist < $dist ) { 
	    $polyasite = $second_polyasite;
	} elsif ($na == 1 && $dist < $second_dist) {
	    $polyasite = $first_polyasite;
	} else {
	    # pick first by default
	    $polyasite = $first_polyasite;
	}
	$DEBUG && debug("na is $na, second na $second_na dist $dist and second_dist $second_dist so setting polya $polyasite (first $first_polyasite, second $second_polyasite)");

    }

#    print STDERR "$polyasite approximate_polya_site $approximate_polya_site strand $strand win $first_winstart-$first_winend and $second_winstart-$second_winend\n";

    return ($polyasite, $na);
}

# sub get_closst_a_around_center  {
    
#     if( $no_aaa_check == 1) {
# 	return (-1, 0);
#     }
    
#     my $approximate_polya_site = shift;
#     my $strand = shift;
#     my $cds_end = shift;
#     my $ppystart = shift;

#     my $seqlen = length($seq->seq());
    
#     # safety checks on the start and end of subseq window

#     my $first_winstart = $approximate_polya_site - $polya_aaa_scanning_window >= 0 ? $approximate_polya_site - $polya_aaa_scanning_window : 0;
#     my $first_winend = $approximate_polya_site;
#     my $second_winstart = $approximate_polya_site;
#     my $second_winend = $approximate_polya_site + $polya_aaa_scanning_window + 1 < $seqlen ? $approximate_polya_site + $polya_aaa_scanning_window + 1 : $seqlen; 
    
#     # usefulnesscheck; warn when window smaller than size due to edge effect..
#     if ( $second_winend - $first_winstart +1 < 2 * $polya_aaa_scanning_window +1) {
# 	$WARNING && debug("WARNING: polya site scanning window smaller (".( $second_winend - $first_winstart +1 )
# 			  .") than given winsize (".(2 * $polya_aaa_scanning_window +1).") presumably due to contig break.");
#     }

#     my $first_aaa_window = lc substr( $seq->seq(), $first_winstart, $first_winend - $first_winstart + 1);
#     my $second_aaa_window = lc substr( $seq->seq(), $second_winstart, $second_winend - $second_winstart + 1);

#     # swap and revcomp windows if reverse strand
#     if( $strand < 0 ) {
# 	my $tmp = revcomp( $first_aaa_window );
# 	$first_aaa_window = revcomp( $second_aaa_window );
# 	$second_aaa_window = $tmp;

# 	$tmp = $first_winstart;
# 	$first_winstart = $second_winstart;
# 	$second_winstart = $tmp;

# 	$tmp = $first_winend;
# 	$first_winend = $second_winend;
# 	$second_winend = $tmp;
#     }

#     my $polyasite;
    
#     my $na = 0;
#     my $dist = 0;
#     my $first_polyasite = 0;

#     # fist win 
#     my @win = split / */, $first_aaa_window;


# #    print STDERR "first win length is ".scalar(@win)." and end - start is ".($first_winend - $first_winstart)."\n";

#   BASE:    for (my $i = @win-1; $i >= 0; $i--) {

#       $base = shift @win;
      
#       if($base ne 'a') {
# 	  next BASE;
#       }

#       $na = 1;
#       if( $strand > 0 ) {
# 	  $first_polyasite = $first_winstart - $i;
#       } else {
# 	  $first_polyasite = $first_winstart + $i;
#       }
#       $dist = $i;
#   }
    
#     if( $na == 0 ) {
# 	$first_polyasite = $approximate_polya_site;
# 	# dist already == 0
#     }

#     # second win
#     @win = split / */, $second_aaa_window;
#     my $second_dist = 0;
#     my $second_na = 0;
#     my $second_polyasite = 0;
    
# #    print STDERR "second win length is ".scalar(@win)." and end - start is ".($second_winend - $second_winstart)."\n";
#   BASEII:    for (my $i = 0; $i < @win ; $i++) {
      
#       $base = shift @win;

#       if($base ne 'a') {
# 	  next BASEII;
#       }
      
#       $second_na = 1;
#       if ( $strand > 0 ) {
# 	  $second_polyasite = $second_winstart + $i;
#       } else {
# 	  $second_polyasite = $second_winend - $i;
#       }
#       $second_dist = $i;
#   }

#     if( $second_na == 0 ) {
# 	$second_polyasite = $approximate_polya_site;
# 	# second_dist already == 0
#     }

#     # closest of first and second?
#     if($na == $second_na) { # 0 or 1
# #	print STDERR "na is $na and $second_na\n";

# 	if ( $second_dist < $dist ) { 
# 	    $polyasite = $second_polyasite;
# 	} else {
# 	    # pick first by default
# 	    $polyasite = $first_polyasite;
# 	}
#     } else {
# #	print STDERR "na is $na and $second_na so setting polya according to the non-zero\n";
# 	if ($na == 0) {
# 	    $polyasite = $second_polyasite;
# 	} else {
# 	    $polyasite = $first_polyasite;
# 	}
#     }

# #    print STDERR "$polyasite approximate_polya_site $approximate_polya_site strand $strand win $first_winstart-$first_winend and $second_winstart-$second_winend\n";

#     return ($polyasite, $na);
# }

sub get_aaa_around  {
    
    if( $no_aaa_check == 1) {
	return (-1, 0);
    }
    
    my $approximate_polya_site = shift;
    my $strand = shift;
    my $cds_end = shift;
    my $ppystart = shift;

    my $seqlen = length($seq->seq());
    
    # safety checks on the start and end of subseq window
    my $winstart;
    my $winend; 

    $winstart = $approximate_polya_site - $polya_aaa_scanning_window;
    ($winstart < 1) && ($winstart = 1);

    if( $strand > 0 ) {
	if( $winstart <= $cds_end ) {
	    $winstart = $cds_end + 1;
	}

	$winend = $winstart + 2 * $polya_aaa_scanning_window + 1;

	if($winend > $ppystart) {
	    $winend = $ppystart;
	}

	if($winend > $seqlen) {
	    $winend = $seqlen;
	    $winstart = $winend - 2 * $polya_aaa_scanning_window;
	    ( $winstart < $cds_end ) && ($winstart = $cds_end);
	    ($winstart < 1) && ($winstart = 1);
	}

    } else {

	$winend = $winstart + 2 * $polya_aaa_scanning_window + 1;

	if($winend > $seqlen) {
	    $winend = $seqlen;
	}

	if( $winend >= $cds_end ) {
	    $winend = $cds_end - 1;
	}

	$winstart = $winend - 2 * $polya_aaa_scanning_window;

	if($winstart < $ppystart) {
	    $winstart = $ppystart;
	}

	($winstart < 1) && ($winstart = 1);
    }
   
    my $aaa_window = lc substr( $seq->seq(), $winstart, $winend-$winstart+1);

    if( $strand < 0 ) {
	$aaa_window = revcomp( $aaa_window );
    }
    my $polyasite;

    my $na = 0;
    if ($aaa_window =~ m/(.*?)(a{3,})(.*)/) {
	if( $strand < 0 ) {
	    $polyasite = $winend - length($3);
	} else {
	    $polyasite = $winstart + length($1);
	}
	$na = length($2);
    } elsif ($aaa_window =~ m/(.*)a{2}(.*)/) {
	if( $strand < 0 ) {
	    $polyasite = $winend - length($2);
	} else {
	    $polyasite = $winstart + length($1);
	}
	$na = 2;
    } elsif ($aaa_window =~ m/(.*)a+(.*)/) {
	if( $strand < 0 ) {
	    $polyasite = $winend - length($2);
	} else {
	    $polyasite = $winstart + length($1);
	}
	$na = 1;
    } else {
	$polyasite = $approximate_polya_site;
    }

#    print STDERR "$polyasite approximate_polya_site $approximate_polya_site strand $strand win $winstart-$winend\n";

    return ($polyasite, $na);
}

sub polypyagscan {
    my $seq = shift;
    
    # forward
    my $thisseq = lc $$seq->seq(); 
    my $seqlen = length($thisseq);

    foreach $direction ( (1, -1) ) {
	($direction == 1) || $revcomp || return;
	($direction == -1) && ($thisseq = revcomp($thisseq)) && $DEBUG && debug("Reverse agsearch");
	# my @start;
	
	my @stop;
	my @start;
	# so, for a possible alternative to the broken regexp (without the lookahead assertion, too short motifs shadow acutal one
	$DEBUG && debug("Direction $direction; finding all py-ag.");
	while($thisseq =~ /([tcyn]+)([agrnx]?[tcynx]*){0,$max_poly_py_mismatches}([tcynx]+)(?=.{0,$exp_polypy_ag_distance}?[amrwvhdxn][grskvdbxn])/g) { # iupac enabled
	    my $start = ($direction == 1) ? $-[0]+1 : $seqlen - $+[0] + 1;
	    my $stop = ($direction == 1) ? $+[0] : $seqlen - $-[0];
	    
	    $DEBUG && print "p(Y) start $start stop $stop, so dist ".($stop - $start +1).".\n";

	    $start <= $stop or fatal("start $start > stop $stop for assigned polypyag on direction $direction."); # assert..	    

            # remember to correct polypycheck for the extra distances..	    	   
	    ($stop - $start + 1 >= $min_poly_py_length)
		|| (length($1) >= $min_qualified_pyrun)
		|| (length($2) >= $min_qualified_pyrun)
		|| (length($3) >= $min_qualified_pyrun)
		|| (length($2) == 0 && (length($1)+length($3) >= $min_qualified_pyrun))
		|| next;

	    # find closest AG dinucleotide, if the stretch survived the sizecheck
            # $DEBUG && debug("Investigating polypy-ag at $start, $stop, $direction...");

	    push @start, $start;
	    push @stop, $stop;
	    $nr_of_ppy_ag_found++;
	}

	my %dupecheck = (); # one for each direction

	$DEBUG && debug("Direction $direction; computing stats for each ppy-ag and adding them as features.");
PPY:	for(my $i = 0; $i < @stop; $i++) {
	    my $last_start = 1;
	    my $last_stop = $seqlen;

	    my $subseq = ($direction == 1) ? substr($thisseq, $stop[$i], $exp_polypy_ag_distance) : substr($thisseq, $seqlen - $start[$i] + 1, $exp_polypy_ag_distance);
	    if($subseq =~ m/(.*?)[amrwvhdxn][grskvdbxn]/) { #iupac enabled (working??) 
		# doesn't this have a problem with n-mismatch >1? eg that the internal mismatch could be "AG", and thus give rise to a very short site?

                # $DEBUG && debug("Found polypy-ag at $start[$i], $stop[$i], $direction (py-ag distance".length($1)." ).");
		my $start = ($direction == 1) ? $start[$i] : $start[$i] - length($1) - 2;
		my $stop = ($direction == 1) ? $stop[$i] + length($1) + 2 : $stop[$i];
                # $DEBUG && $direction == -1 && debug("Setting start at $start, stop at $stop: py-ag distance ".length($1)." polypy length "
		# .($stop[$i]-$start[$i]+1)." [max_poly_py_mismatches = $max_poly_py_mismatches, exp_polypy_ag_distance = $exp_polypy_ag_distance]");		
		

		# Check for duplicates
		# This could be controversial with respect to statistical properties of ppy stretches etcs
		# but for the purpouse of finding UTRs and testing if genes fit the model this is sufficient
		# and cpucycle saving..

		if($disable_dupe_checking == 1) {
		    # $DEBUG && $WARNING && print STDERR "DEBUG WARNING duplicate ppy-ag site checking disabled by user.\n";
		} else {
		   if($direction == 1) {
		       if (defined($dupecheck{$stop})) {
			   # $DEBUG && $WARNING && print STDERR "WARNING dupe found and discarded at $start-$stop ($direction).\n";   
			   next PPY; 
		       } else {
			   $dupecheck{$stop} = 1;
		       }
		   } else {
		       if (defined($dupecheck{$start})) {
			   # $DEBUG && $WARNING && print STDERR "WARNING dupe found and discarded at $start-$stop ($direction).\n";   
			   next PPY;
		       } else {
			   $dupecheck{$start} = 1;
		       }		       
		   }
	       }
		
		$$seq->add_SeqFeature(new Bio::SeqFeature::Generic ( -start => $start, -end => $stop, -strand => $direction, 
								     -primary=>'pPyAg', -source=>"papya_$version",
								     -tag => { 
									 note => $max_poly_py_mismatches." max mismatches, ".length($1)." py-ag distance, ".($stop[$i]-$start[$i]+1)." polypy length.", 
									 pyagdist => length($1),
									 polypylen => ($stop[$i]-$start[$i]+1)
									 }));
		
	    } else {
		# This would be very odd indeed given the lookahead assertion in the regexp, but better safe..
		$DEBUG && debug("No ag within $exp_polypy_ag_distance for polypy at $start[$i], $stop[$i], $direction");
	    }
	}
    }
}

sub gff3_string {
    my $feat = shift;
   
    my $strand = ($feat->strand == 1)?'+':'-'; 
    my $phase = '.';
    
    my $attribute_string ="";
   
    foreach my $tag ($feat->get_all_tags) {
	if($attribute_string ne "") {
	    $attribute_string .= ";";
	}
	#nb if you have multiple freetext tags they need doublequotes, right?
        $values=join(",",$feat->get_tag_values($tag));
#if($values =~ /[^0-9eE.-]+/
	$attribute_string.="$tag=$values";
    }

    if(!$feat->has_tag('score') ) {
	$feat->add_tag_value('score','.');
    }
    #seq is not an argument.. =)
    $gff3_string = join("\t", $seq->display_id, $feat->source_tag, $feat->primary_tag,$feat->start,$feat->end,$feat->score,$strand,$phase,$attribute_string);
    return $gff3_string;
}

sub usage {
    my $msg = shift;
    
    !defined($msg) && ($msg = "");
    $msg ne ""  && ($msg .= "\n");
    
    $msg .= "Usage: splicemodel.pl\n\t--infile <genbank filename>\n";
    $msg .= "\t[--report]\t\t\t\tModel results summary & gff partitioned into different categories\n";
    $msg .= "\t[--reportfile <filename>]\t\tImplies --report.\n";
    $msg .= "\t[--gff <filename>]\t\t\tExport conforming CDSs and signals in GFF\n";
    $msg .= "\t[--fasta <filename>]\t\t\tExport UTRs in fasta format\n";
    $msg .= "\t[--gff-utrs <filename>]\t\t\tExport UTRs in GFF\n";
    $msg .= "\t[--statistics <statistics-filename>]\n";
    $msg .= "\t[--rout]\t\t\t\tStatistics in R format\n";
    $msg .= "\t[--splitstats]\t\t\t\tSplit statisticsfile\n";
    $msg .= "\t[--uorf]\n";
    $msg .= "\t[--uorf-file <filename>]\t\tuORFs in GFF, implies --uorf.\n";
    $msg .= "\t[--print-int-alt-sites-file <filename>]\tDisregarded internal alternative splice sites\n";

    $msg .= "\t[--long]\t\t\t\tLongest possible UTRs only\n";
    $msg .= "\t[--short]\t\t\t\tShortest possible UTRs only\n";
    $msg .= "\t[--no-aaa-check]\t\t\tDisable AAA check around polyA-site\n";
    $msg .= "\t[--polya-model <fix|aaa|center>]\tPolya placement algorithm\n";
    $msg .= "\t[--coupled-polya]\t\tOnly SLadd sites are considered as polya signals - scan stops at next ds gene\n";

    $msg .= "\t[--all-ppyag]\t\t\t\tDisable removal of duplicate ppy for each ag\n";
    $msg .= "\t[--polypy-ag <distance>]\n";
    $msg .= "\t[--5pag-atg <distance>]\n";
    $msg .= "\t[--min-polypy-length <int>]\n";
    $msg .= "\t[--min-qual-polypy-run <int>]\t\tShortest stretch of contigous pY to direct qualify pY-site.\n";
    $msg .= "\t[--max-polypy-mismatches <int>]\n";
    $msg .= "\t[--max-stop-ag-distance <int>]\n";
    $msg .= "\t[--polya-pyag <int>]\n";
    $msg .= "\t[--half-polya-win <int>]\t\tWindow to check for AAA/AA/A around polyA-site\n";
    $msg .= "\t[--no-revcomp]\t\t\t\tDo not examine reverse strand\n";
    $msg .= "\t[--debug]\n";
    $msg .= "\t[--no-warn]\n";

    print STDERR $msg;
    exit;
}

