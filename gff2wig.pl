#!/usr/bin/perl -w
#
# Daniel Nilsson, 2010-02-26
#

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $DEBUG = 0;
my $WARNING = 1;

my $tag_map_feature_name = "nucleotide_match";
my $tag_map_file_name="";
my $sequence_file_name = "";
my $output_file_name = "";

if (@ARGV == 0 ) {
    print STDERR "USAGE: gff2wig.pl -s ref_sequence_file -t tag_map_gff [-f tag_feature_name] [-o output_file]\n";
    exit;
}

# preprocess command line arguments
while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
	
        if ($arg eq '-t') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-t requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $tag_map_file_name = $next_arg;
            }
        }

        if ($arg eq '-f') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-f requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $tag_map_feature_name = $next_arg;
            }
        }

        if ($arg eq '-s') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-s requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $sequence_file_name = $next_arg;
            }
        }
        if ($arg eq '-o') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-o requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $output_file_name = $next_arg;
            }
        }

    }
}

# check input
if($tag_map_file_name eq "") {
    print STDERR "No tag map gff file name given. Bailing out.\n";
    exit 1;
}

if($sequence_file_name eq "") {
    print STDERR "No sequence file name given. Bailing out.\n";
    exit 1;
}

# open files
my $outputfh = *STDOUT;
if($output_file_name ne "") {
    local *OUTPUT;
    open OUTPUT, ">$output_file_name";
    $outputfh = *OUTPUT;
} else {
    $DEBUG && print "Writing gff out to stdout.\n";
}

$WARNING && print STDERR "Read sequence file: ";
my $fastain = Bio::SeqIO->new('-format' => 'fasta', -file => $sequence_file_name);
my $seq = $fastain->next_seq();
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Processing ", $seq->display_id,".\n";

my @coverage = (0)x$seq->length;
my $tag_maps=0;

$WARNING && print STDERR "Read tag GFF file: ";
open(TAGGFF,$tag_map_file_name) or die "Could not open $tag_map_file_name\n";
while(my $r=<TAGGFF>) {
    chomp $r;
    parse_gff_row($r, $tag_map_feature_name);
}
close TAGGFF;
$WARNING && print STDERR "OK\n";

my $ref_name = $seq->display_id;
my $ref_size = $seq->length;

print $outputfh "# Coverage per pos, unbinned uncomressed and unsmoothed..\n";
print $outputfh "# A total of $tag_maps tags mapped to this sequence.\n";
print $outputfh "track type=wiggle_0 name=\"$tag_map_feature_name\" desc=\"$tag_map_feature_name for $ref_name from $tag_map_file_name\" visibility=full transform=logtransform\n";
print $outputfh "fixedStep chrom=$ref_name start=1 step=1 span=1\n";
for ( my $j = 0; $j < @coverage ; $j++ ) {
    print $outputfh $coverage[$j], "\n";
}

exit 0;

sub parse_gff_row {
    my $r = shift;
    my $filter_feature_name = shift;

    if($r=~/^#+.+/) {
	#meta line..
	return;
    }

    my @r=split(/\t+/, $r);

    my $ref_name = $r[0];
    my $feature_source = $r[1];
    my $feature_name = $r[2];    
    my $ref_start = $r[3]; # GFF 1-based inclusive
    my $ref_end = $r[4];
    my $score = $r[5];
    my $strand = $r[6];
    my $frame = $r[7];
    my $feature_desc = $r[8];    

#    if( $feature_name eq "contig" ) {
#	$DEBUG && print "landmark $feature_name - $ref_name - ignored.\n";
#	return;
#    }

    my $direction = 0;
    if( defined($strand) && $strand eq '+' ) {
	$direction = 1;
    } elsif ( defined($strand) && $strand eq '-')  {
	$direction = -1;
    } else {
	if ($r =~ /^>/) {
	    #included fasta seq name.. Ignore for now..
	    $WARNING && print STDERR "WARNING: found fasta seq in gff for $r. Ignoring.\n";
	} else {
	    $WARNING && print STDERR "WARNING: weird strand assignment! Column error? $r\n";
	}
	return;
    }
    
    if ($feature_name eq $filter_feature_name) {
		
	# callback
	for( my $j = $ref_start - 1; $j < $ref_end ; $j++ ) {
	    $coverage[$j]+= $direction;
	}
	$tag_maps++;
    }
}



