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

my $DEBUG = 0;
my $WARNING = 1;

my $tab_file_name ="";

my $tab_file_name_given =0;
my $lib_total_tags=0;
my $summed_lib_total_tags = 0;
my $output_file_name = "";

my $ss_count_threshold = 5;

my $internal_early = -900;
my $external_max_utr_len = 2000;

my $major_f_threshold = 0.6;

my $print_all = 0;

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
                $lib_total_tags = $next_arg;
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

	if($arg eq '-A') {
	    $print_all = 1;
	}

	if ($arg eq '-i') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-i requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $internal_early = -$next_arg;
            }
	}

	if ($arg eq '-u') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-u requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $external_max_utr_len = $next_arg;
            }
	}

	if ($arg eq '-c') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-c requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $ss_count_threshold = $next_arg;
            }
	}

    } else {
	if (!$tab_file_name_given) {
	    $tab_file_name = $arg;
	    $tab_file_name_given = 1;
      	} else {
	    print "Unknown/extra argument: $arg. Sorry, I'm afraid something went wrong on that command line. Bailing out.\n";
	}
    }
}

if ($tab_file_name_given == 0) {
    print STDERR "A tab file name must be given. Bailing out.\n";
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

$WARNING && print STDERR "Read tab: ";
open(TABFILE,$tab_file_name) or die "Could not open $tab_file_name\n";
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
    
    $summed_lib_total_tags += $sladd_total_count;
}
close TABFILE;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Compute test statistics: ";
my $M = $lib_total_tags;

if($lib_total_tags == 0) {
    $M = $summed_lib_total_tags;
}

my @names = keys %liba_tc;

foreach my $name (@names) {

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

    # also tally the fractions of each total

    # actually, nevermind, do the major one

    if( $liba_tc{$name} >0 ) {
	$libasite_f{$name} = $libasite{$liba_sites[0]} / $liba_tc{$name};
    } else {
	$libasite_f{$name} = 1;
    }
}
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Output: ";

print $outputfh "###Name\tStart\tStrand\tLibA_total_Raw\tLibA_total_Normalised\tLibA_sladd_ind_norm_count\tLibA_sladd_ind_fpUTRlen\tLibA_sladd_ind_pos\tDesc\n";

my $liba_less_than_f_in_major=0;
my $internal_only=0;
my $external_only=0;
my $ext_and_int=0;
my $weirdo=0;

foreach my $name (@names) { 
    my @individual_ss_counts=();
    my @individual_fputrs=();
    my @individual_pos=();

    if(  $libasite_f{$name}<$major_f_threshold) {
	
	if ($liba_sladd_individual_count{$name} =~ m/,/) {
	    @individual_ss_counts=split(/,+/,$liba_sladd_individual_count{$name});
	}

	if ($liba_sladd_individual_fputrs{$name} =~ m/,/) {
	    @individual_fputrs=split(/,+/,$liba_sladd_individual_fputrs{$name});
	}

	if ($liba_sladd_individual_pos{$name} =~ m/,/) {
	    @individual_pos=split(/,+/,$liba_sladd_individual_pos{$name});
	} else {
	    push @individual_pos, $liba_sladd_individual_pos{$name};
	}

	my $has_internal=0;
	my $has_external=0;
               
	for($i=0; $i<@individual_ss_counts; $i++ ) {
	    
	    # normalise the tag count while we're at it
	    $individual_ss_counts[$i] = sprintf("%.0f", $individual_ss_counts[$i] / $lib_total_tags * 1000000 );

	    # internal site, interesting for signal trimming
	    if( $individual_fputrs[$i] < 0 && $individual_fputrs[$i] > $internal_early && $individual_ss_counts[$i] >=$ss_count_threshold) {
		$has_internal++;
	    }
	    
	    # require one external site, <2kb
	    # actually, also all-internal can be interesting, but perhaps we dont bother with them right now?
	    if( $individual_fputrs[$i] > 0 && $individual_fputrs[$i] <=  $external_max_utr_len && $individual_ss_counts[$i] >=$ss_count_threshold ) {
		$has_external++;
	    }
	}

	if(($print_all==0 && $has_internal >0 && $has_external>0) or ($print_all == 1 && ($has_internal >1 || $has_external > 1))) {
	    print $outputfh $name,"\t",$liba_start{$name},"\t",$liba_strand{$name};
	    print $outputfh "\t",$liba_tc{$name},"\t",$liba_tc_normalised{$name},"\t\"",join(",", @individual_ss_counts),"\"\t\"",$liba_sladd_individual_fputrs{$name},"\"\t\"",$liba_sladd_individual_pos{$name};
	    print $outputfh "\"\t\"",$liba_feature_desc{$name},"\"\n";
	}
	 
	if($has_internal>0 && $has_external>0) {
	    $ext_and_int++;
	} elsif($has_internal >0) {
	    $internal_only++;
	} elsif($has_external >0) {
	    $external_only++;
	} else {
	    $weirdo++;
	}

	if($libasite_f{$name}<$major_f_threshold) {
	    $liba_less_than_f_in_major++;
	}
    }
}

print $outputfh "=== Total genes with below $major_f_threshold counts in major site (max UTR len $external_max_utr_len, min sas count $ss_count_threshold, require internal-early before bp $internal_early for lib :$tab_file_name:$liba_less_than_f_in_major:$external_only:$ext_and_int:$internal_only:$weirdo (accepted ",($external_only+$internal_only+$ext_and_int),")\n";

if($print_all == 1) {
    print $outputfh "=== Printing all sites with counts in major site below threshold ($major_f_threshold), and more than one site within $external_max_utr_len bp or within the early threshold ($internal_early) internally over the count threshold of $ss_count_threshold.\n";
} else {
    print $outputfh "=== Printing sites with counts in major site below threshold ($major_f_threshold), and at least one site within $external_max_utr_len bp as well as at least one internally within the early threshold ($internal_early), over the count threshold of $ss_count_threshold.\n";
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
