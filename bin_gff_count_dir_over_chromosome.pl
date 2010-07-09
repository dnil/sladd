#!/usr/bin/perl -w

# Daniel Nilsson, 2009-05-01
#

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $DEBUG = 0;
my $WARNING = 1;

my $tagged_gff_file_name ="";
my $tag_map_feature_name = "nucleotide_match";
my $tag_map_file_name="";
my $sequence_file_name = "";
my $output_file_name = "";

while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-g') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-g requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $tagged_gff_file_name = $next_arg;
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

if($tagged_gff_file_name eq "") {
    print STDERR "No tagged gff file name given. Bailing out.\n";
    exit 1;
}

if($sequence_file_name eq "") {
    print STDERR "No sequence file name given. Bailing out.\n";
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

$WARNING && print STDERR "Read sequence file: ";
my $fastain = Bio::SeqIO->new('-format' => 'fasta', -file => $sequence_file_name);
my $seq = $fastain->next_seq();
$WARNING && print STDERR "OK\n";

# presumably the crunched version..
$WARNING && print STDERR "Parse tagged GFF file for tags: ";
open(TAGGEDGFF,$tagged_gff_file_name) or die "Could not open $tagged_gff_file_name\n";
while(my $r=<TAGGEDGFF>) {
    chomp $r;
    parse_gff_row($r, $tag_map_feature_name);
}
close TAGGEDGFF;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Retrieve feature references: ";
my @tag_maps;

foreach my $feat ($seq->all_SeqFeatures()) {
    if( $feat->primary_tag eq $tag_map_feature_name) {
	push @tag_maps, \$feat;
    } 
}
print STDERR "Found ".scalar(@tag_maps)." unique tags.\n";
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Sort tags: ";
@tag_maps = sort {($$a->strand == 1 && $$a->start <=> $$b->start) || ($$a->strand == -1 && $$a->end <=> $$b->end)} @tag_maps;
$WARNING && print STDERR "OK\n";

#print $outputfh "##gff-version 3\n";

my $ref_name = $seq->display_id;
my $ref_size = $seq->length;

print $outputfh "###seq_name\tstart\tend\n";
print $outputfh "###$ref_name\t1\t$ref_size\n";
print $outputfh "###bin_start\tbin_end\tbin_count(TPM)\n";

#my @bins;

my $bin_size=2000;

my $current_bin_start=1; # 1based inclusive coords
my $current_bin_end= $current_bin_start + $bin_size -1;
if( $current_bin_end > $ref_size) {
    $current_bin_end = $ref_size;
}
my $bin_count=0;
my $bin_dir =0;
my $total_tag_count=0;
foreach my $tag (@tag_maps) {
    $total_tag_count += $$tag->score;

#    print STDERR "tag at ".$$tag->start." with total count now at $total_tag_count.\n";

    if( ($$tag->strand == 1) && ($$tag->start <= $current_bin_end ) or (($$tag->strand == -1) && ($$tag->end <= $current_bin_end)) ) {

	$bin_count += $$tag->score;

	$bin_dir += $$tag->score * $$tag->strand;

#	print STDERR "added to bin count for $current_bin_start $current_bin_end: $bin_count.\n";
    } else {
	
	while( ($$tag->strand == 1) && ($$tag->start > $current_bin_end ) or (($$tag->strand == -1) && ($$tag->end > $current_bin_end)) ) {
	    # init next bin
	    # first print old
	    print $outputfh "$current_bin_start\t$current_bin_end\t$bin_count\t$bin_dir\n";
	    
	    # move start
	    $current_bin_start = $current_bin_start + $bin_size;
	    if($current_bin_start >= $ref_size ) {
		print STDERR "Oops, didn't expect that; $current_bin_start>=$ref_size $current_bin_end $bin_count. You sure you've got the right file names lined up together on that commandline, pal?\n";
		exit;
	    }
	    #move end
	    $current_bin_end = $current_bin_start + $bin_size -1 ;
	    if( $current_bin_end > $ref_size) {
		$current_bin_end = $ref_size;
	    }
	    # zero count
	    $bin_count=0;
	    $bin_dir=0;
	}

	if( ($$tag->strand == 1) && ($$tag->start <= $current_bin_end ) or (($$tag->strand == -1) && ($$tag->end <= $current_bin_end)) ) {
	    $bin_count += $$tag->score;
	    $bin_dir += $$tag->score * $$tag->strand;
	} else {
	    print STDERR "Now that wasn't expected either: did we run out of chromosome??\n";
	}
    }
}

# then see if there is chromosome left at the end of the tags
while($current_bin_end < $ref_size ) {	
    print $outputfh "$current_bin_start\t$current_bin_end\t$bin_count\t$bin_dir\n";
    $current_bin_start = $current_bin_start + $bin_size;
    $current_bin_end = $current_bin_start + $bin_size -1;
    if( $current_bin_end > $ref_size) {
	$current_bin_end = $ref_size;
    }
    $bin_count =0;
    $bin_dir=0;
}

# and print the last one 
print $outputfh "$current_bin_start\t$current_bin_end\t$bin_count\t$bin_dir\n";

print $outputfh "###tag_tagfile\tn_unique_tags\ttotal_counts\n";
print $outputfh "###$tagged_gff_file_name\t".scalar(@tag_maps)."\t$total_tag_count\n";

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
	$WARNING && print STDERR "WARNING: weird strand assignment! Column error? $r\n";	
	return;
    }
    
    if ($feature_name eq $filter_feature_name) {
	my $feat = new Bio::SeqFeature::Generic ( -start => $ref_start, -end => $ref_end, -strand => $direction,
						  -primary=>$feature_name,
						  -source=>$feature_source );
	if(defined($score)) {
	    $feat->add_tag_value('score',$score);
	}

	if (defined($feature_desc) && $feature_desc ne "") {
	    my @feature_strings;
	    if( $feature_desc =~ /[^=";\s]+\s*=\s*\"[^"\s]+\"\s*(?:;|$ )/ ) {
		(@feature_strings) = ($feature_desc =~ /((?:[^=";\s]+\s*=\s*\"[^"]+\"\s*)|(?:[^=";\s]+\s*=\s*[^"=\s;]+\s*)|(?:[^=";]+\s*))(?:;|$)/g );
		$DEBUG && print "feature desc has \"-enclosed value.\n";
	    } else {
		@feature_strings = split(/[;]+/, $feature_desc); # q&d, could have problems with ; in ""..
	    }
	    $DEBUG && print $feature_desc,"\n";
	    foreach my $feature_string (@feature_strings) {
		my ($key, $value);
		
		if( $feature_string =~ /^\s*([^=\s]+)\s*$/ ) { # boolean tags.. like "pseudo;"
		    $key = $1; 
		    $value = 1;
		    $DEBUG && print "boolean tag.\n";
		} else {
		    ($key, $value) = ($feature_string =~ m/([^=]+)\s*=\s*(.+)/);
		}

		$value =~ s/"//g; 
		$value =~ s/;//g;

		$DEBUG && print "key: $key is value: $value\n";
		if($value=~m/,/) {
		    my @values = split (/,+/, $value);
		    foreach my $tmp_value (@values) {
			$feat->add_tag_value($key,$tmp_value);
		    }
		} else {
		    $feat->add_tag_value($key,$value);
		}
	    }
	}

	$seq->add_SeqFeature($feat);
    }
}
