#!/usr/bin/perl -w
#
# Daniel Nilsson, 2009-04-23
#

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $DEBUG = 0;
my $WARNING = 1;

my $tagged_gff_file_name ="";
my $gene_feature_name = "gene";
my $sequence_file_name = "";
my $gff_out = "";
my $output_file_name = "";

my $library_size = 1000000;

if (@ARGV == 0 ) {
    print STDERR "USAGE: tagged_rnaseq_gff_to_counts_tab.pl -g tagged_gff [-c gene_feature_name <$gene_feature_name>] [-N library_size <$library_size>] [-o output_file <STDOUT>]\n";
    exit;
}

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

        if ($arg eq '-c') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-c requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $gene_feature_name = $next_arg;
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

        if ($arg eq '-N') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-N requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $library_size = $next_arg;
            }
        }
    }  
}

if($tagged_gff_file_name eq "") {
    print STDERR "No tagged gff file name given. Bailing out.\n";
    exit 1;
}

my $outputfh = *STDOUT;
if($output_file_name ne "") {
    local *OUTPUT;
    open OUTPUT, ">$output_file_name";
    $outputfh = *OUTPUT;
} else {
    $DEBUG && print "Writing output to stdout.\n";
}

my $total_mapped_count =0;
my $total_nu_over_gene_length = 0;

$WARNING && print STDERR "Read tagged GFF file: ";
open(TAGGEDGFF,$tagged_gff_file_name) or die "Could not open $tagged_gff_file_name\n";
while(my $r=<TAGGEDGFF>) {
    chomp $r;
    parse_gff_row($r, $gene_feature_name);
}
close TAGGEDGFF;
$WARNING && print STDERR "OK\n";

# adjust from previous library size estimate to current total mapped count, ie mapped in this particular feature, which is a subset of all mapped or all in library..
if ($total_mapped_count > 0 ) {
    $total_nu_over_gene_length = $total_nu_over_gene_length * $library_size / $total_mapped_count;
} else {
    print STDERR "WARNING: no tag counts found for gff file $tagged_gff_file_name. Was this really the expected result?\n";
}
print $outputfh "Total mapped count:$total_mapped_count\nTotal nu over gene length:$total_nu_over_gene_length\n";

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
	
	# specific accounting
	if ($feat->has_tag('tagcount')) {
	    my ($gene_tag_count) = ($feat->get_tag_values('tagcount'));
	    $total_mapped_count += $gene_tag_count;
	}
	
	if ($feat->has_tag('tagcountnu')) {
	    my ($gene_tag_count_nu) = ($feat->get_tag_values('tagcountnu'));
	    my $gene_length = $feat->length;
	    if ( $gene_length > 0 ) {
		$total_nu_over_gene_length += $gene_tag_count_nu / $gene_length;
	    } else {
		print STDERR "WARNING: found 0 length feature ", $feature_name, " at ", $ref_start, " in ", $ref_name, ". No nu/len update added for this feature.\n";
	    }
	}
    }
}
