#!/usr/bin/perl -w
#
# Daniel Nilsson, 2009-04-23
#
# parse count tab files
# normalise counts
# determine relative expression, careful to weigh in confidence in count levels as well..
# perhaps simply do a test?
# especially catch large diffs in alt splice levels,
# print new tab with comparison scores.

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

use Bio::SAGE::Comparison;
use Text::NSP::Measures::2D::Fisher2::twotailed;

my $DEBUG = 0;
my $WARNING = 1;

my $tab_file_name_1 ="";
my $tab_file_name_2 ="";

my $tab_file_name_1_given =0;
my $tab_file_name_2_given =0;

my $lib_1_total_tags=0;
my $lib_2_total_tags=0;
my $summed_lib_1_total_tags = 0;
my $summed_lib_2_total_tags = 0;

my $output_file_name = "";

my $major_f_threshold = 0.6;

# my $show_top_N = 500;

while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-o') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-o requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $output_file_name = $next_arg;
            }
        }

	if ($arg eq '-M') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-M requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $lib_1_total_tags = $next_arg;
            }
	}
	
	if ($arg eq '-f') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-f requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $major_f_threshold = $next_arg;
            }
	}

	if ($arg eq '-N') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-N requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $lib_2_total_tags = $next_arg;
            }
	}
    } else {
	if (!$tab_file_name_1_given) {
	    $tab_file_name_1 = $arg;
	    $tab_file_name_1_given = 1;
	} elsif(!$tab_file_name_2_given) {
	    $tab_file_name_2 = $arg;
	    $tab_file_name_2_given = 1;
	} else {
	    print "Unknown/extra argument: $arg. Sorry, I'm afraid something went wrong on that command line. Bailing out.\n";
	}
    }
}

if ($tab_file_name_1_given == 0or $tab_file_name_2_given ==0) {
    print STDERR "Two tab file names must be given for a comparison. Bailing out.\n";
    exit 1;
}

my $outputfh = *STDOUT;
if($output_file_name ne "") {
    local *OUTPUT;
    open OUTPUT, ">$output_file_name";
    $outputfh = *OUTPUT;
} else {
    $DEBUG && print "Writing gff out to stdout.\n";
}

$WARNING && print STDERR "Read tab 1: ";
open(TABFILE,$tab_file_name_1) or die "Could not open $tab_file_name_1\n";
while(my $r=<TABFILE>) {
    chomp $r;
    
    my ($name, $start, $strand, $sladd_total_count, $sladd_total_normalised, $sladd_individual_count, $sladd_individual_pos, $sladd_individual_fputrs, $feature_desc)=parse_tab_row($r);
    if ($name eq "comment line") {
        #strange row..
	next;
    }

    $liba_tc{$name} = $sladd_total_count;
    $liba_tc_normalised{$name} = $sladd_total_normalised;
    
    $liba_start{$name} = $start;
    $liba_strand{$name} = $strand;
    $liba_sladd_individual_count{$name} = $sladd_individual_count;
    $liba_sladd_individual_pos{$name} = $sladd_individual_pos;
    $liba_sladd_individual_fputrs{$name} = $sladd_individual_fputrs;

    $liba_feature_desc{$name} = $feature_desc;
    
    $summed_lib_1_total_tags += $sladd_total_count;
}
close TABFILE;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Read tab 2: ";
open(TABFILE,$tab_file_name_2) or die "Could not open $tab_file_name_2\n";
while(my $r=<TABFILE>) {
    chomp $r;

    my ($name, $start, $strand, $sladd_total_count, $sladd_total_normalised, $sladd_individual_count,$sladd_individual_pos,$sladd_individual_fputrs,$feature_desc)=parse_tab_row($r);
    if ($name eq "comment line") {
        #strange row..
	next;
    }

    $libb_tc{$name} = $sladd_total_count;
    $libb_tc_normalised{$name} = $sladd_total_normalised;

    $libb_sladd_individual_count{$name} = $sladd_individual_count;
    $libb_sladd_individual_pos{$name} = $sladd_individual_pos;
    $libb_sladd_individual_fputrs{$name} = $sladd_individual_fputrs;

    $summed_lib_2_total_tags += $sladd_total_count;

}
close TABFILE;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Compute test statistics: ";
my $M = $lib_1_total_tags;
my $N = $lib_2_total_tags;

if($lib_1_total_tags == 0) {
    $M = $summed_lib_1_total_tags;
}

if ($lib_2_total_tags == 0) {
    $N = $summed_lib_2_total_tags;
}

my @names = keys %liba_tc;

foreach my $name (@names) {

    my ( $p, $sign ) = Bio::SAGE::Comparison::calculate_significance( $liba_tc{$name},$libb_tc{$name}, $M, $N, 1 );

    if( $sign == -1 ) { 
	$liba_tc_p{$name} = $p;
	$libb_tc_p{$name} = 1;
    } elsif ($sign == +1) {
	$liba_tc_p{$name} = 1;
	$libb_tc_p{$name} = $p;
    } elsif( $sign == 0 ) {
	if( $p <= 0.01 ) {
	    die( "Same expression should never be significant!" ); 
	}

	# nevermind, not significant anyway..
	$liba_tc_p{$name} = $p;
	$libb_tc_p{$name} = $p;
    }

    # split splice site list
    # (use two hashes with SL position as key to be able to easily compare, also with missing values) 
    my @libacounts = split (/,+/,$liba_sladd_individual_count{$name});
    my @libapos = split (/,+/,$liba_sladd_individual_pos{$name});
    my %libasite;

    for (my $i=0; $i < @libacounts; $i++) {
	$libasite{$libapos[$i]}=$libacounts[$i];
    }

    my @liba_sites = keys %libasite;
    @liba_sites = sort {$libasite{$b}<=>$libasite{$a}} @liba_sites ;

    my @libbcounts = split(/,+/,$libb_sladd_individual_count{$name});
    my @libbpos = split(/,+/,$libb_sladd_individual_pos{$name});
    my %libbsite;

    for (my $i=0; $i < @libbcounts; $i++) {
	$libbsite{$libbpos[$i]}=$libbcounts[$i];	
    }

    my @libb_sites = keys %libbsite;
    @libb_sites = sort {$libbsite{$b}<=>$libbsite{$a}} @libb_sites ;

    # make sure each top hit is represented in the other set

    for (my $i=0; $i< @libb_sites; $i++) {
	if( !defined($libasite{$libb_sites[$i]}) ) {
	    $libasite{$libb_sites[$i]} = 0;
	}
    }

    for (my $i=0; $i< @liba_sites; $i++) {
	
	if( !defined($libbsite{$liba_sites[$i]}) ) {
	    $libbsite{$liba_sites[$i]} = 0;
	}
	
    }

    # also tally the fractions of each total
    # actually, nevermind, do the major one
    if( $liba_tc{$name} >0 ) {
	$libasite_f{$name} = $libasite{$liba_sites[0]} / $liba_tc{$name};
    } else {
	$libasite_f{$name} = 1;
    }

    if ($libb_tc{$name}>0) {
	$libbsite_f{$name} = $libbsite{$libb_sites[0]} / $libb_tc{$name};
    } else {
	$libbsite_f{$name} = 1;
    }

    # compare sites between libraries
#    if (@liba_sites > 0) {
    $liba_tc_altsp{$name}="";
    $liba_tc_altsp_min{$name}=1;

#    $libb_tc_altsp{$name}="";
#    $libb_tc_altsp_min{$name}=1;

#	for (my $i=0; $i< @liba_sites; $i++) {

    
    if( ($liba_tc{$name} > 0 && $libb_tc{$name} > 0) ) {	    # otherwise its not alternative splicing, only diff control. COULD draw this further and disregard any where total count is significantly up or down in either stage, leaving only the highly similar total count ones.
	my $i = 0; #only the major sites
	if($liba_sites[$i] != $libb_sites[$i]) {

	    my $n_sa_la = $libasite{$liba_sites[$i]};
	    my $n_sa_lb = $libbsite{$liba_sites[$i]};
	    my $n_sb_la = $libasite{$libb_sites[$i]};
	    my $n_sb_lb = $libbsite{$libb_sites[$i]};

	    my $liba_tot = $n_sa_la + $n_sb_la;
	    my $sa_tot = $n_sa_la + $n_sa_lb;
	    my $table_tot = $n_sa_la + $n_sa_lb + $n_sb_la + $n_sb_lb;
	    
	    # Fischer test, two tailed
	    my $p = calculateStatistic( n11=>$n_sa_la,n1p=>$sa_tot ,np1=> $sa_tot, npp=>$table_tot);
	    
#	    my ($p,$sign) = Bio::SAGE::Comparison::calculate_significance( $libasite{$liba_sites[$i]},$libbsite{$liba_sites[$i]}, $M, $N,1 );
	    if($i > 0) {
		$liba_tc_altsp{$name}.=",";
	    }
	    $liba_tc_altsp{$name}.=$p;
	    
	    if($liba_tc_altsp_min{$name} > $p) {
		$liba_tc_altsp_min{$name} = $p;
	    }

#	    ($p,$sign) = Bio::SAGE::Comparison::calculate_significance( $libasite{$libb_sites[$i]},$libbsite{$libb_sites[$i]}, $M, $N,1 );
#	    if($i > 0) {
#		$liba_tc_altsp{$name}.=",";
#	    }
#	    $libb_tc_altsp{$name}.=$p;
#	    
#	    if($libb_tc_altsp_min{$name} > $p) {
#		$libb_tc_altsp_min{$name} = $p;
#	    }
	}
    }
    
    # compute margins..
#    $twotailed_value = calculateStatistic( n11=>$liba_sitea,n1p=>$liba_total,np1=>$sitea_total,npp=>$table_total);


    # perhaps also test within each lib, if A-C statistic non-significant between the top and runner up§
    
    
}
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Output: ";

print $outputfh "\n*** significantly different (P < 1e-5) major alternative splice sites\n\n";

print $outputfh "###Name\tStart\tStrand\tLibA_altsplice_P\tLibA_total_Raw\tLibA_total_Normalised\tLibA_sladd_Individual_count\tLibA_sladd_Individual_fpUTRlen\tLibA_sladd_Individual_pos\tLibB_total_Raw\tLibB_total_Normalised\tLibB_sladd_Individual_count\tLibB_sladd_Individual_fpUTRlen\tLibB_sladd_Individual_pos\tDesc\n";

foreach my $name (@names) { 
    if(defined($liba_tc_altsp{$name}) and ( ($liba_tc_altsp_min{$name} <0.00001 ) ) ) {
	print $outputfh $name,"\t",$liba_start{$name},"\t",$liba_strand{$name};
	print $outputfh "\t",sprintf("%.3e",$liba_tc_altsp{$name});
	print $outputfh "\t",$liba_tc{$name},"\t",$liba_tc_normalised{$name},"\t\"",$liba_sladd_individual_count{$name},"\"\t\"",$liba_sladd_individual_fputrs{$name},"\"\t\"",$liba_sladd_individual_pos{$name};
	print $outputfh "\"\t",$libb_tc{$name},"\t",$libb_tc_normalised{$name},"\t\"",$libb_sladd_individual_count{$name},"\"\t\"",$libb_sladd_individual_fputrs{$name},"\"\t\"",$libb_sladd_individual_pos{$name};
#	print $outputfh "\t",sprintf("%.3e",$libb_tc_altsp{$name}); #keep?
	print $outputfh "\"\t\"",$liba_feature_desc{$name},"\"\n";
    }
}

print $outputfh "\n*** genes with major splice site with <60% of the total counts for that gene\n\n";

print $outputfh "###Name\tStart\tStrand\tLibA_total_Raw\tLibA_total_Normalised\tLibA_tc_overrep_acP\tLibA_sladd_Individual_count\tLibA_sladd_Individual_fpUTRlen\tLibA_sladd_Individual_pos\tLibA_altsplice_P\tLibB_total_Raw\tLibB_total_Normalised\tLibB_tc_overrep_acP\tLibB_sladd_Individual_count\tLibB_sladd_Individual_fpUTRlen\tLibB_sladd_Individual_pos\tDesc\n";

my @names_in_count_order = sort { $liba_tc{$b}+$libb_tc{$b} <=> $liba_tc{$a}+$libb_tc{$a} } @names;

my $liba_less_than_f_in_major=0;
my $libb_less_than_f_in_major=0;

foreach my $name (@names_in_count_order) { 
    if(  $libasite_f{$name}<$major_f_threshold or $libbsite_f{$name}<$major_f_threshold) {
	print $outputfh $name,"\t",$liba_start{$name},"\t",$liba_strand{$name};
	print $outputfh "\t",$liba_tc{$name},"\t",$liba_tc_normalised{$name},"\t",sprintf("%.3e",$liba_tc_p{$name}),"\t\"",$liba_sladd_individual_count{$name},"\"\t\"",$liba_sladd_individual_fputrs{$name},"\"\t\"",$liba_sladd_individual_pos{$name};
	print $outputfh "\"\t",$liba_tc_altsp{$name};
	print $outputfh "\t",$libb_tc{$name},"\t",$libb_tc_normalised{$name},"\t",sprintf("%.3e",$libb_tc_p{$name}),"\t\"",$libb_sladd_individual_count{$name},"\"\t\"",$libb_sladd_individual_fputrs{$name},"\"\t\"", $libb_sladd_individual_pos{$name};
#	print $outputfh "\t",$libb_tc_altsp{$name};
	print $outputfh "\"\t\"",$liba_feature_desc{$name},"\"\n";

	if($libasite_f{$name}<$major_f_threshold) {
	    $liba_less_than_f_in_major++;
	}

	if($libbsite_f{$name}<$major_f_threshold) {
	    $libb_less_than_f_in_major++;
	}
    }
}

print $outputfh "=== Total genes with below $major_f_threshold counts in major site for lib :$tab_file_name_1:$liba_less_than_f_in_major:$tab_file_name_2:$libb_less_than_f_in_major\n";

# alternative sign up tab 

foreach my $name (@names) { 
    if(defined($liba_tc_p{$name}) and ( ($liba_tc_p{$name} <0.00001 ) ) ) {
	push @significanta, $name;
    }

    if(defined($libb_tc_p{$name}) and ( ($libb_tc_p{$name} <0.00001 ) ) ) {
	push @significantb, $name;
    }
}

$WARNING && print STDERR "Sort for highest fold change: ";
my @topa = sort { $liba_tc{$b}/($libb_tc{$b}==0?0.9:$libb_tc{$b}) <=> $liba_tc{$a}/($libb_tc{$a}==0?0.9:$libb_tc{$a}) } @significanta;
my @topb = sort { $libb_tc{$b}/($liba_tc{$b}==0?0.9:$liba_tc{$b}) <=> $libb_tc{$a}/($liba_tc{$a}==0?0.9:$liba_tc{$a}) } @significantb;
$WARNING && print STDERR "OK\n";

#$WARNING && print STDERR "Sort for most significant: ";
#@topa = sort { $liba_tc_p{$a} <=> $liba_tc_p{$b} } @names;
#@topb = sort { $libb_tc_p{$a} <=> $libb_tc_p{$b} } @names;
#$WARNING && print STDERR "OK\n";

#then print most significant ones.
#my $show_top_a = (scalar(@topa)<$show_top_N)?scalar(@topa):$show_top_N;
# print $outputfh "\n*** $show_top_a most A-C significantly upregulated in $tab_file_name_1.\n\n";

print $outputfh "\n*** significantly (P < 1e-5) upregulated genes (up in $tab_file_name_1) ordered by fold change\n\n";

print $outputfh "###Name\tStart\tStrand\tLibA_total_Raw\tLibA_total_Normalised\tLibA_tc_overrep_acP\tLibA_sladd_Individual_count\tLibA_sladd_Individual_pos\tLibB_total_Raw\tLibB_total_Normalised\tLibB_tc_overrep_acP\tLibB_sladd_Individual_count\tLibB_sladd_Individual_pos\tDesc\n";

#liba..

my $show_top_a = scalar(@topa);

for (my $i=0; $i < $show_top_a ; $i++) {

    my $name = $topa[$i];

    print $outputfh $name,"\t",$liba_start{$name},"\t",$liba_strand{$name};
    print $outputfh "\t",$liba_tc{$name},"\t",$liba_tc_normalised{$name},"\t",sprintf("%.3e",$liba_tc_p{$name}),"\t\"",$liba_sladd_individual_count{$name},"\"\t\"",$liba_sladd_individual_pos{$name};
    print $outputfh "\"\t",$libb_tc{$name},"\t",$libb_tc_normalised{$name},"\t",sprintf("%.3e",$libb_tc_p{$name}),"\t\"",$libb_sladd_individual_count{$name},"\"\t\"",$libb_sladd_individual_pos{$name};
    print $outputfh "\"\t\"",$liba_feature_desc{$name},"\"\n";
}

#my $show_top_b = (scalar(@topb)<$show_top_N)?scalar(@topb):$show_top_N;

my $show_top_b = scalar(@topb);

#print $outputfh "\n*** $show_top_b most A-C significantly upregulated from $tab_file_name_2.\n\n";

print $outputfh "\n*** significantly (P < 1e-5) downregulated genes (up in $tab_file_name_2) ordered by fold change\n\n";

print $outputfh "###Name\tStart\tStrand\tLibA_total_Raw\tLibA_total_Normalised\tLibA_tc_overrep_acP\tLibA_sladd_Individual_count\tLibA_sladd_Individual_pos\tLibB_total_Raw\tLibB_total_Normalised\tLibB_tc_overrep_acP\tLibB_sladd_Individual_count\tLibB_sladd_Individual_pos\tDesc\n";

for (my $i=0; $i < $show_top_b ; $i++) {

    my $name = $topb[$i];

    print $outputfh $name,"\t",$liba_start{$name},"\t",$liba_strand{$name};
    print $outputfh "\t",$liba_tc{$name},"\t",$liba_tc_normalised{$name},"\t",$liba_tc_p{$name},"\t\"",$liba_sladd_individual_count{$name},"\"\t\"",$liba_sladd_individual_pos{$name};
    print $outputfh "\"\t",$libb_tc{$name},"\t",$libb_tc_normalised{$name},"\t",$libb_tc_p{$name},"\t\"",$libb_sladd_individual_count{$name},"\"\t\"",$libb_sladd_individual_pos{$name};
    print $outputfh "\"\t",$liba_feature_desc{$name},"\n";
}
$WARNING && print STDERR "OK\n";

exit 0;

sub parse_tab_row {
    my $r = shift;

    if($r=~/^#+.+/) {
	#meta line..
	return ("comment line"); 
    }

    my @r=split(/\t+/, $r);

    my $name = $r[0];
    my $start = $r[1]; # GFF 1-based inclusive
    my $strand = $r[2];
    my $sladd_total_count = $r[3];
    my $sladd_normalised_count = $r[4];
    
    my $sladd_individual_count = $r[5];
    my $sladd_individual_pos = $r[6];

    my $sladd_individual_fputrs = $r[7];
    my $major_fputr = $r[8];
    my $major_fputr_count = $r[9];
    my $shortest_fputr = $r[10];
    my $longest_fputr = $r[11];

    my $feature_desc = $r[12];

    my $direction = 0;
    if( defined($strand) && $strand eq '+' ) {
	$direction = 1;
    } elsif ( defined($strand) && $strand eq '-')  {
	$direction = -1;
    } else {
	$WARNING && print STDERR "WARNING: weird strand assignment! Column error? $r\n";
	return ("comment line");
    }

    # populate a hash/object directly? nah, better done after return
    return ($name, $start, $strand, $sladd_total_count, $sladd_normalised_count,$sladd_individual_count,$sladd_individual_pos,$sladd_individual_fputrs,$feature_desc);

}
