#!/usr/bin/perl -w
#
# Daniel Nilsson, 2009-04-23
#

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $DEBUG = 0;
my $WARNING = 1;

my $gene_gff_file_name ="";
my $gene_feature_name = "gene";
my $tag_map_feature_name = "nucleotide_match";
my $tag_map_file_name="";
my $sequence_file_name = "";
my $gff_out = "";
my $output_file_name = "";
my $library_size = 1000000;

if (@ARGV == 0 ) {
    print STDERR "USAGE: per_gene_coverage_info.pl -s ref_sequence_file -g reference_gff [-c gene_feature_name] -t tag_map_gff [-f tag_feature_name] [-o output_file] [-N library_size]\n";
    exit;
}

# preprocess command line arguments
while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-g') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-g requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $gene_gff_file_name = $next_arg;
            }
        }
	
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

        if ($arg eq '-c') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-c requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $gene_feature_name = $next_arg;
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

# check input
if($gene_gff_file_name eq "") {
    print STDERR "No gene gff file name given. Bailing out.\n";
    exit 1;
}

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

$WARNING && print STDERR "Read gene GFF file: ";
open(GENEGFF,$gene_gff_file_name) or die "Could not open $gene_gff_file_name\n";
while(my $r=<GENEGFF>) {
    chomp $r;
    parse_gff_row($r, $gene_feature_name);
}
close GENEGFF;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Read tag GFF file: ";
open(TAGGFF,$tag_map_file_name) or die "Could not open $tag_map_file_name\n";
while(my $r=<TAGGFF>) {
    chomp $r;
    parse_gff_row($r, $tag_map_feature_name);
}
close TAGGFF;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Retrieve feature references: ";
my @genes;
my @tag_maps;

foreach my $feat ($seq->all_SeqFeatures()) {
    if( $feat->primary_tag eq $tag_map_feature_name) {
	push @tag_maps, \$feat;
    } 
    elsif( $feat->primary_tag eq $gene_feature_name) {
	push @genes, \$feat;
    }
}
$WARNING && print STDERR "OK\n";

print $outputfh "##gff-version 3\n";
my $ref_name = $seq->display_id;
my $ref_size = $seq->length;    
print $outputfh "##sequence-region $ref_name 1 $ref_size\n";
print $outputfh "$ref_name\tTriTrypDB\tcontig\t1\t$ref_size\t.\t+\t.\tName=$ref_name\n";

if (@tag_maps == 0) {
    print STDERR "WARNING: Number of tags is 0. Error in input? Continuing, but be advised.\n";
}

if (@genes == 0) {
    print STDERR "Number of genes is 0. Error in input? Bailing out.\n";
    print "Number of genes is 0 .\n";   
    exit 1;
}

$WARNING && print STDERR "Sort genes: ";
@genes = sort { $$a->strand <=> $$b->strand || ($$a->strand == 1 && $$a->start <=> $$b->start) || ($$a->strand == -1 && $$b->end <=> $$a->end)} @genes;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Sort tags: ";
@tag_maps = sort { $$a->strand <=> $$b->strand || ($$a->strand == 1 && $$a->start <=> $$b->start) || ($$a->strand == -1 && $$b->end <=> $$a->end)} @tag_maps;
$WARNING && print STDERR "OK\n";

my $current_tag = 0;
my $mappedtagcount = 0;

my $ngenes = scalar(@genes);
my $ntags = scalar(@tag_maps);

$WARNING && print STDERR "Tag features: ";
#$DEBUG && print "N genes = $ngenes, N tags = $ntags.\n";

for (my $i = 0; $i < $ngenes ; $i++) {

    my $current_gene = $genes[$i];
    my $current_tag_map = $tag_maps[$current_tag];

    my $current_strand = $$current_gene->strand;    

    $DEBUG && print "DEBUG: Current gene (nr $i) ".($$current_gene->gff_string)."\n";
    
    if ( $ntags == 0 ) {
	# there were no tags to map
	last;
    }

    if ( $current_tag == $ntags-1 ) {
	if ($$current_tag_map->has_tag('mapped') or $$current_tag_map->has_tag('mapped_intergenic') ) {
	    # if the last tag has been mapped already, there is noting to assign to this gene, so abort?
	    last;
	}
    }

    if($current_strand == -1) {

	my $keep_checking_tags = 1;

	if ($i == $ngenes-1 or $$current_gene->strand != ${$genes[$i+1]}->strand) {
	    # last gene on this strand
	    # keep checking all tags on this strand and mark them as potential extras??
	    # no, polyA-only splicesites will occurr eleswhere (any end of TU, also when many TUs present on same chr)
	    $DEBUG && print STDERR "last round on $current_strand..\n";
	}

	while ( $keep_checking_tags ) {

	    $DEBUG && print "DEBUG: current_tag (nr $current_tag): ".($$current_tag_map->gff_string)."\n";

	    if ($$current_tag_map->strand == $current_strand) {

		if($$current_tag_map->end >= $$current_gene->end) {
		    
		    if( $$current_tag_map->start <= $$current_gene->end ) {			
			# tag overlaps rev strand gene end (ie start codon)

			# good tag, attach to gene.
			if(! $$current_gene->has_tag('tagged_part')) {
			    $$current_gene->add_tag_value('tagged_part', 1);
			}

			# update gene score: add the covered number of bases to the running count
			my $coverage=0;
			if($$current_gene->has_tag('coverage')) {
			    ($coverage) = ($$current_gene->get_tag_values('coverage'));
			    $$current_gene->remove_tag('coverage');
			}
			# what if the tag covers the entire gene? This would not be common with present day tag lengths, but still, let's deal with it.
			my $start = ($$current_tag_map->start > $$current_gene->start) ? $$current_tag_map->start : $current_gene->start;
			$coverage += ( $$current_gene->end - $start +1 );
			$$current_gene->add_tag_value('coverage',sprintf("%.0f", $coverage));

			my $tagcount=0;
			if($$current_gene->has_tag('tagcount')) {
			    ($tagcount) = ($$current_gene->get_tag_values('tagcount'));
			    $$current_gene->remove_tag('tagcount');
			}
			$tagcount++;
			$$current_gene->add_tag_value('tagcount',sprintf("%.0f", $tagcount));
			$mappedtagcount++;

			# track tag as used - store atg coord for want of a better gene description
			$$current_tag_map->add_tag_value('mapped',$$current_gene->end);
			
#			if($$current_gene->has_tag('Name') and ($$current_gene->get_tag_values('Name') ne "cds")) {
#			    $$current_tag_map->add_tag_value('mapped_name',$$current_gene->get_tag_values('Name'));
#			} elsif($$current_gene->has_tag('ID')) {
#			    $$current_tag_map->add_tag_value('mapped_name',$$current_gene->get_tag_values('ID'));
#			}
		    } else {
			$$current_tag_map->add_tag_value('mapped_intergenic',$$current_gene->end);
		    }

		    if ($current_tag < $ntags-1) {
			$current_tag++;
			$current_tag_map = $tag_maps[$current_tag];
		    } elsif ($current_tag == $ntags-1) {
			# out of tags, next gene please.
			$keep_checking_tags = 0;
			$DEBUG && print "DEBUG: Rather done with gene (nr $i) ".($$current_gene->gff_string)."\n";
		    }

		} elsif ($$current_tag_map->end < $$current_gene->end && $$current_tag_map->end >= $$current_gene->start) {
		   
		    if($$current_tag_map->start > $$current_gene->start) {
			# tag is within the gene

			# good tag, attach to gene.
			if(! $$current_gene->has_tag('tagged')) {
			    $$current_gene->add_tag_value('tagged', 1);
			}
			# update gene score: add the covered number of bases to the running count
			my $coverage=0;
			if($$current_gene->has_tag('coverage')) {
			    ($coverage) = ($$current_gene->get_tag_values('coverage'));
			    $$current_gene->remove_tag('coverage');
			}
			$coverage += $$current_tag_map->length;
			$$current_gene->add_tag_value('coverage',sprintf("%.0f", $coverage));

		    } else {
			# tag overlaps rev strand start of gene (stop codon)
			
			# good tag, attach to gene
			if(! $$current_gene->has_tag('tagged_part')) {
			    $$current_gene->add_tag_value('tagged_part', 1);
			}

			# update gene score: add the *covered* number of bases to the running count (tag start to gene end)
			my $coverage=0;
			if($$current_gene->has_tag('coverage')) {
			    ($coverage) = ($$current_gene->get_tag_values('coverage'));
			    $$current_gene->remove_tag('coverage');
			}
			$coverage += ( $$current_tag_map->end - $$current_gene->start +1 );
			$$current_gene->add_tag_value('coverage',sprintf("%.0f", $coverage));

		    }

		    my $tagcount=0;
		    if($$current_gene->has_tag('tagcount')) {
			($tagcount) = ($$current_gene->get_tag_values('tagcount'));
			$$current_gene->remove_tag('tagcount');
		    }
		    $tagcount++;
		    $$current_gene->add_tag_value('tagcount',sprintf("%.0f", $tagcount));
		    $mappedtagcount++;

		    #track tag as used  - store atg coord for want of a better gene description
		    $$current_tag_map->add_tag_value('mapped',$$current_gene->end);

#		    if($$current_gene->has_tag('Name') and ($$current_gene->get_tag_values('Name') ne "cds")) {
#			$$current_tag_map->add_tag_value('mapped_name',$$current_gene->get_tag_values('Name'));
#		    }elsif($$current_gene->has_tag('ID')) {
#			$$current_tag_map->add_tag_value('mapped_name',$$current_gene->get_tag_values('ID'));
#		    }

		    if ($current_tag < $ntags-1) {
			$current_tag++;
			$current_tag_map = $tag_maps[$current_tag];
		    } elsif($current_tag == $ntags-1) {
			# out of tags, next gene please.
			$keep_checking_tags = 0;
			$DEBUG && print "DEBUG: Rather done with gene (nr $i) ".($$current_gene->gff_string)."\n";
		    }

		} else {
		    #tag is on the right strand, but out of gene: keep tag and move to checking next gene
		    $keep_checking_tags = 0;
		    $DEBUG && print "DEBUG: Rather done with gene (nr $i) ".($$current_gene->gff_string)."\n";
		}

	    } else {
		#so, what to do when current_tag is not on the right strand any longer? simply keep counting up the genes, I guess. 
		#tags are also strand-sorted after all.
		$keep_checking_tags = 0;
		$DEBUG && print "DEBUG: Rather done with gene (nr $i) ".($$current_gene->gff_string)."\n";
	    }
	}

    } elsif ($current_strand == 1) {
	
	my $keep_checking_tags = 1;
	
	if ($i == $ngenes-1 or $current_strand != ${$genes[$i+1]}->strand) {
	    # now checking last gene on + strand.. well, nothing much to do about it?
	    $DEBUG && print STDERR "last round on $current_strand..\n";
	}

	while( $keep_checking_tags ) {
	    $DEBUG && print "DEBUG: current_tag (nr $current_tag): ".($$current_tag_map->gff_string)."\n";
	    if ($$current_tag_map->strand == $current_strand) {

		if($$current_tag_map->start < $$current_gene->start ) {
		    if( $$current_tag_map->end >= $$current_gene->start ) {
			# partially on CDS - count towards coverage within.

			# good tag, attach to gene.
			if(! $$current_gene->has_tag('tagged_part')) {
			    $$current_gene->add_tag_value('tagged_part', 1);
			}

			# update gene score: add the covered number of bases to the running count
			my $coverage=0;
			if($$current_gene->has_tag('coverage')) {
			    ($coverage) = ($$current_gene->get_tag_values('coverage'));
			    $$current_gene->remove_tag('coverage');
			}
			# what if the tag covers the entire gene? This would not be common with present day tag lengths, but still, let's deal with it.
			my $end = ($$current_tag_map->end < $$current_gene->end) ? $$current_tag_map->end : $current_gene->end;
			$coverage += ( $end - $$current_gene->start +1 );
			$$current_gene->add_tag_value('coverage',sprintf("%.0f", $coverage));
			
			my $tagcount=0;
			if($$current_gene->has_tag('tagcount')) {
			    ($tagcount) = ($$current_gene->get_tag_values('tagcount'));
			    $$current_gene->remove_tag('tagcount');
			}
			$tagcount++;
			$$current_gene->add_tag_value('tagcount',sprintf("%.0f", $tagcount));
			$mappedtagcount++;

			#track tag as used - store atg coord for want of a better gene description
			$$current_tag_map->add_tag_value('mapped',$$current_gene->start);
			
#			if($$current_gene->has_tag('Name') and ($$current_gene->get_tag_values('Name') ne "cds")) {
#			    $$current_tag_map->add_tag_value('mapped_name',$$current_gene->get_tag_values('Name'));
#			} elsif($$current_gene->has_tag('ID')) {
#			    $$current_tag_map->add_tag_value('mapped_name',$$current_gene->get_tag_values('ID'));
#			}
		    } else {
			$$current_tag_map->add_tag_value('mapped_intergenic',$$current_gene->end);
		    }
		    
		    if ($current_tag < $ntags-1) { 
			$current_tag++;
			$current_tag_map = $tag_maps[$current_tag];
		    } elsif($current_tag == $ntags-1) {
			# out of tags, next gene please.
			$keep_checking_tags = 0;
			$DEBUG && print "DEBUG: Rather done with gene (nr $i) ".($$current_gene->gff_string)."\n";
		    }

		} elsif ($$current_tag_map->start >= $$current_gene->start && $$current_tag_map->start <= $$current_gene->end ) {
		   
		    if($$current_tag_map->end <= $$current_gene->end) {
			# tag is within the gene

			# good tag, attach to gene.
			if(! $$current_gene->has_tag('tagged')) {
			    $$current_gene->add_tag_value('tagged', 1);
			}
			
			# update gene score: add the covered number of bases to the running count
			my $coverage=0;
			if($$current_gene->has_tag('coverage')) {
			    ($coverage) = ($$current_gene->get_tag_values('coverage'));
			    $$current_gene->remove_tag('coverage');
			}
			$coverage += $$current_tag_map->length;
			$$current_gene->add_tag_value('coverage',sprintf("%.0f", $coverage));
			
		    } else {
			# tag overlaps end of gene

			# good tag, attach to gene.
			if(! $$current_gene->has_tag('tagged_part')) {
	 		    $$current_gene->add_tag_value('tagged_part', 1);
			}

			# update gene score: add the *covered* number of bases to the running count (tag start to gene end)
			my $coverage=0;
			if($$current_gene->has_tag('coverage')) {
			    ($coverage) = ($$current_gene->get_tag_values('coverage'));
			    $$current_gene->remove_tag('coverage');
			}
			$coverage += ( $$current_gene->end - $$current_tag_map->start +1 );
			$$current_gene->add_tag_value('coverage',sprintf("%.0f", $coverage));
		    }

		    my $tagcount=0;
		    if($$current_gene->has_tag('tagcount')) {
			    ($tagcount) = ($$current_gene->get_tag_values('tagcount'));
			    $$current_gene->remove_tag('tagcount');
			}
		    $tagcount++;
		    $$current_gene->add_tag_value('tagcount',sprintf("%.0f", $tagcount));
		    $mappedtagcount++;
		    
		    #track tag as used  - store atg coord for want of a better gene description
		    $$current_tag_map->add_tag_value('mapped',$$current_gene->start);

#		    if($$current_gene->has_tag('Name') and ($$current_gene->get_tag_values('Name') ne "cds")) {
#			$$current_tag_map->add_tag_value('mapped_name',$$current_gene->get_tag_values('Name'));
#		    } elsif($$current_gene->has_tag('ID')) {
#			$$current_tag_map->add_tag_value('mapped_name',$$current_gene->get_tag_values('ID'));
#		    }

		    if ($current_tag < $ntags-1) {
			$current_tag++;
			$current_tag_map = $tag_maps[$current_tag];
		    } elsif($current_tag == $ntags-1) {
			# out of tags, next gene please.
			$keep_checking_tags = 0;
			$DEBUG && print "DEBUG: Rather done with gene (nr $i) ".($$current_gene->gff_string)."\n";
		    }

		} else {
		    #tag is on the right strand, but out of gene: keep tag and move to checking next gene
		    $keep_checking_tags = 0;		    
		    $DEBUG && print "DEBUG: Rather done with gene (nr $i) ".($$current_gene->gff_string)."\n";
		}
	    } else {
		#so, what to do when current_tag is not on the right strand? like for the last gene in a TU before a SS? 
		#or more correctly, for the last gene on the chr - strand before starting with 

#		$keep_checking_tags = 0;
		
		if ($current_tag < $ntags-1) { 
		    $current_tag++;
		    $current_tag_map = $tag_maps[$current_tag];
		} else {
		    # out of tags, next gene please.
		    $keep_checking_tags = 0;
		    $DEBUG && print "DEBUG: Rather done with gene (nr $i) ".($$current_gene->gff_string)."\n";
		}
	    }
	} 
    }
}

$WARNING && print STDERR "OK\n";

# how many tags remain unassigned?
# how many genes had internal assignments? 

$WARNING && print STDERR "Count tag results: ";
my $tag_maps_gene_mapped = 0;
my $tag_maps_gene_mapped_internal = 0;
my $tag_maps_mapped_intergenic = 0;

# print $outputfh "##gff-version 3\n";
# my $ref_name = $seq->display_id;
# my $ref_size = $seq->length;
# print $outputfh "##sequence-region $ref_name 1 $ref_size\n";
# print $outputfh "$ref_name\tTriTrypDB\tcontig\t1\t$ref_size\t.\t+\t.\tName=$ref_name\n";

foreach my $tag_map (@tag_maps) {
    if($$tag_map->has_tag('mapped')) {
	$tag_maps_gene_mapped++;
    } elsif( $$tag_map->has_tag('mapped_intergenic') ) {
	$tag_maps_mapped_intergenic++;
    }
#    print $outputfh $$tag_map->gff_string,"\n";
#    print $outputfh gff3_string($$tag_map),"\n";
}
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Count gene results: ";
my $gene_mapped = 0;

#my $total_nu_over_length=0;

foreach my $gene (@genes) {
    if($$gene->has_tag('tagged') or $$gene->has_tag('tagged_part') ) {
	$gene_mapped++;

	# compute per column coverage
	my $coverage = 0;

	if($$gene->has_tag('coverage')) {
	    ($coverage) = ($$gene->get_tag_values('coverage'));
	    $$gene->remove_tag('coverage');
	}

	my $gene_length = $$gene->end - $$gene->start + 1;
	my $gene_column_coverage = $coverage / $gene_length;

	# "raw" per gene column coverage
	$$gene->add_tag_value('coverage',sprintf("%.0f", $gene_column_coverage)); # or do decimal?

	# normalise per gene per column to TPM
	my $normalised_coverage = $gene_column_coverage / $library_size * 1000000; # NPM column fraction
	$$gene->add_tag_value('normcoverage',sprintf("%.0f", $normalised_coverage));

	my $gene_tag_count = 0;

	if($$gene->has_tag('tagcount')) {
	    ($gene_tag_count) = ($$gene->get_tag_values('tagcount'));
	}

	my $gene_tag_count_nu = $gene_tag_count/$library_size; # this is rather aproximate, as the numer of reads on the actual CDSs would be less than the reported library size
	
	if($$gene->has_tag('tagcountnu')) {
	    $$gene->remove_tag('tagcountnu');
	}
	$$gene->add_tag_value('tagcountnu',$gene_tag_count_nu);

#	$total_nu_over_len += $gene_tag_count_nu / $gene_length;
    }

#    print $outputfh $$gene->gff_string,"\n";

}

#loop again for final adjustments
foreach my $gene (@genes) {

    my $gene_tag_count_nu = 0;

    if($$gene->has_tag('tagcountnu')) {
	($gene_tag_count_nu) = ($$gene->get_tag_values('tagcountnu'));
    }
	
    my $gene_length = $$gene->length ;
#	my $gene_fraction_of_transcripts = $gene_tag_count_nu / $gene_length / $total_nu_over_len ; # need complete set of genes for this..
#	my $tpm = $gene_fraction_of_transcripts * 1000000;
    my $npm = $gene_tag_count_nu * 1000000;

#	if($$gene->has_tag('TPM')) {
#	    $$gene->remove_tag('TPM');
#	}
#	$$gene->add_tag_value('TPM',sprintf("%.0f", $tpm));

    if($$gene->has_tag('NPM')) {
	$$gene->remove_tag('NPM');
    }
    $$gene->add_tag_value('NPM',sprintf("%.0f", $npm));

    # set score to tags per million
    if($$gene->has_tag('score')) {
	$$gene->remove_tag('score');
    }
    $$gene->add_tag_value('score',sprintf("%.0f", $npm));
        
    print $outputfh gff3_string($$gene),"\n";
}
$WARNING && print STDERR "OK\n";

print "There were ".scalar(@tag_maps)." tag_maps for ".scalar(@genes)." genes.\n";
print "$gene_mapped genes (feature type $gene_feature_name) had RNAseq reads mapped onto them.\n";
print "$tag_maps_gene_mapped tags were mapped to genes (feature type $gene_feature_name), $tag_maps_mapped_intergenic were labeled intergenic.\n";

exit 0;

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
    
    $gff3_string = join("\t", $seq->display_id, $feat->source_tag, $feat->primary_tag,$feat->start,$feat->end,$feat->score,$strand,$phase,$attribute_string);
    return $gff3_string;
}

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
#		$DEBUG && print "feature desc has \"-enclosed value.\n";
	    } else {
		@feature_strings = split(/[;]+/, $feature_desc); # q&d, could have problems with ; in ""..
	    }
#	    $DEBUG && print $feature_desc,"\n";
	    foreach my $feature_string (@feature_strings) {
		my ($key, $value);
		if( $feature_string =~ /^\s*([^=\s]+)\s*$/ ) { # boolean tags.. like "pseudo;"
		    $key = $1; 
		    $value = 1;
#		    $DEBUG && print "boolean tag.\n";
		} else {
		    ($key, $value) = ($feature_string =~ m/([^=]+)\s*=\s*(.+)/);
		}				
#		$DEBUG && print "key: $key is value: $value\n";

		$value =~ s/"//g; 
		$value =~ s/;//g;

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



