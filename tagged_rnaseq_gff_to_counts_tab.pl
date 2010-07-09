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
my $total_tag_count = $library_size;

my $total_nu_over_gene_length = 0.0001;

if (@ARGV == 0 ) {
    print STDERR "USAGE: tagged_rnaseq_gff_to_counts_tab.pl -g tagged_gff [-c gene_feature_name] -s sequence_file [-o output_file] [-T total_mapped_tags] [-N library_size] [-n nu_tau_length_fraction_correction]\n";
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
        if ($arg eq '-s') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-s requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $sequence_file_name = $next_arg;
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

        if ($arg eq '-n') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-n requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $total_nu_over_gene_length = $next_arg;
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

        if ($arg eq '-T') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-T requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $total_tag_count = $next_arg;
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

$WARNING && print STDERR "Read tagged GFF file: ";
open(TAGGEDGFF,$tagged_gff_file_name) or die "Could not open $tagged_gff_file_name\n";
while(my $r=<TAGGEDGFF>) {
    chomp $r;
    parse_gff_row($r, $gene_feature_name);
}
close TAGGEDGFF;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Retrieve feature references: ";
my @genes;

foreach my $feat ($seq->all_SeqFeatures()) {
    if( $feat->primary_tag eq $gene_feature_name) {
	push @genes, \$feat;
    }
}
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Agglomerate counts: ";

# use total count of tags mapped to genes (from input, need to go through all gffs to calculate this) rather than library size as total tags in transcriptome (we are ignorant of UTR/actual transcripts at this stage..)
foreach my $gene (@genes) {

    my $gene_tag_count = 0;

    if($$gene->has_tag('tagcount')) {
	($gene_tag_count) = ($$gene->get_tag_values('tagcount'));
    }

    my $gene_tag_count_nu = $gene_tag_count/$total_tag_count; #read latter from combined raw gene gff. one round for reading tagcounts & producing a gene_nu / gene_len
    
    if($$gene->has_tag('tagcountnu')) {
	$$gene->remove_tag('tagcountnu');
    }
    $$gene->add_tag_value('tagcountnu',$gene_tag_count_nu);
    
    my $gene_length = $$gene->length ;
    my $gene_fraction_of_transcripts = $gene_tag_count_nu / $gene_length / $total_nu_over_gene_length ; #read latter from combined raw gene gff
    my $tpm = $gene_fraction_of_transcripts * 1000000;
    my $npm = $gene_tag_count_nu * 1000000;

    if($$gene->has_tag('TPM')) {
	$$gene->remove_tag('TPM');
    }
    $$gene->add_tag_value('TPM',sprintf("%.0f", $tpm));

    if($$gene->has_tag('NPM')) {
	$$gene->remove_tag('NPM');
    }
    $$gene->add_tag_value('NPM',sprintf("%.0f", $npm));

    my $normalised_coverage = 0;
    if ($$gene->has_tag('normcoverage')) {
	($normalised_coverage) = ($$gene->get_tag_values('normcoverage'));
    }

    my $coverage = 0;
    if ($$gene->has_tag('coverage')) {
	($coverage) = ($$gene->get_tag_values('coverage'));
    }
    
    my $name=$$gene->start." ".$$gene->end;
    if($$gene->has_tag('Name') and ($$gene->get_tag_values('Name') ne "cds")) {
	@namearr = $$gene->get_tag_values('Name');
	$name= $namearr[0];
    } elsif($$gene->has_tag('ID')) {
	@namearr = $$gene->get_tag_values('ID');
	$name= $namearr[0];
    }
	
    my $description = ".";
    if($$gene->has_tag('description')) {
	my @descarr = $$gene->get_tag_values('description');
	$description= $descarr[0];
    } elsif($$gene->has_tag('product')) {
	my @descarr = $$gene->get_tag_values('product');
	$description= $descarr[0];
    }
    if($$gene->has_tag('note')) {
	my @descarr = $$gene->get_tag_values('note');
	$description.= " ".join (" ", @descarr);
    }

    my $strand = ($$gene->strand == 1)?'+':'-';
    my $atg_pos =  ($$gene->strand == 1)?$$gene->start:$$gene->end;

    print $outputfh $name,"\t",$atg_pos,"\t",$strand,"\t",sprintf("%.0f", $npm),"\t",sprintf("%.0f", $tpm);
    print $outputfh "\t.\t.\t.\t.\t.\t.\t.";
    print $outputfh "\t",$description,"\n";
}
$WARNING && print STDERR "OK\n";

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


# sub parse_gff_row {
#     my $r = shift;
#     my $filter_feature_name = shift;

#     if($r=~/^#+.+/) {
# 	#meta line..
# 	return; 
#     }

#     my @r=split(/\t+/, $r);

#     my $ref_name = $r[0];
#     my $feature_source = $r[1];
#     my $feature_name = $r[2];    
#     my $ref_start = $r[3]; # GFF 1-based inclusive
#     my $ref_end = $r[4];
#     my $score = $r[5];
#     my $strand = $r[6];
#     my $frame = $r[7];
#     my $feature_desc = $r[8];    

#     my $direction = 0;
#     if( defined($strand) && $strand eq '+' ) {
# 	$direction = 1;
#     } elsif ( defined($strand) && $strand eq '-')  {
# 	$direction = -1;
#     } else {
# 	$WARNING && print STDERR "WARNING: weird strand assignment! Column error? $r\n";	
# 	return;
#     }
    
#     if ($feature_name eq $filter_feature_name) {
# 	my $feat = new Bio::SeqFeature::Generic ( -start => $ref_start, -end => $ref_end, -strand => $direction,
# 						  -primary=>$feature_name,
# 						  -source=>$feature_source );

# 	if (defined($feature_desc) && $feature_desc ne "") {
# 	    my @feature_strings = split(/[;]+/, $feature_desc); # q&d, could have problems with ; in ""..
# 	    $DEBUG && print $feature_desc,"\n";
# 	    foreach my $feature_string (@feature_strings) {
# 		my ($key, $value) = ($feature_string =~ m/^\s*([^=\s]+)\s*=\s*(.+)\s*$/);
		
#                 #print "key: $key is value: $value\n"
# 		if($value=~m/,+/) {
# 		    my @values = split (/,+/, $value);
# #		    print "..so add ",join(" ", @values),"\n";
# 		    foreach my $tmp_value (@values) {
# 			$feat->add_tag_value($key,$tmp_value);
# 		    }
# 		} else {
# 		    $feat->add_tag_value($key,$value);
# 		}

# 	    }
# 	}

# 	$seq->add_SeqFeature($feat);       	
#     }
# }
