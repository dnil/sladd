#!/usr/bin/perl -w
#
# Daniel Nilsson, 2009-05-01
#

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $DEBUG = 0;
my $WARNING = 1;

my $tagged_gff_file_name ="";
my $gene_feature_name = "gene";
my $tag_map_feature_name = "nucleotide_match";
my $tag_map_file_name="";
my $sequence_file_name = "";
my $gff_out = "";
my $output_file_name = "";

my $library_size = 1000000;

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

        if ($arg eq '-N') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-N requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $library_size = $next_arg;
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

$WARNING && print STDERR "Parse tagged GFF file for genes: ";
open(TAGGEDGFF,$tagged_gff_file_name) or die "Could not open $tagged_gff_file_name\n";
while(my $r=<TAGGEDGFF>) {
    chomp $r;
    parse_gff_row($r, $gene_feature_name);
}
close TAGGEDGFF;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Parse tagged GFF file again for tags: ";
open(TAGGEDGFF,$tagged_gff_file_name) or die "Could not open $tagged_gff_file_name\n";
while(my $r=<TAGGEDGFF>) {
    chomp $r;
    parse_gff_row($r, $tag_map_feature_name);
}
close TAGGEDGFF;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Retrieve feature references: ";
my @genes;
my @tag_maps;

foreach my $feat ($seq->all_SeqFeatures()) {
    if( $feat->primary_tag eq $tag_map_feature_name) {
	push @tag_maps, \$feat;
    } 
    if( $feat->primary_tag eq $gene_feature_name) {
	push @genes, \$feat;
    }
}
print STDERR "Found ".scalar(@tag_maps)." tags and ".scalar(@genes)." genes.\n";
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Sort genes: ";
@genes = sort { $$a->strand <=> $$b->strand || ($$a->strand == 1 && $$a->start <=> $$b->start) || ($$a->strand == -1 && $$b->end <=> $$a->end)} @genes;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Sort tags: ";
@tag_maps = sort { $$a->strand <=> $$b->strand || ($$a->strand == 1 && $$a->start <=> $$b->start) || ($$a->strand == -1 && $$b->end <=> $$a->end)} @tag_maps;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Agglomerate counts: ";

print $outputfh "##gff-version 3\n";

my $ref_name = $seq->display_id;
my $ref_size = $seq->length;

print $outputfh "##sequence-region $ref_name 1 $ref_size\n";
print $outputfh "$ref_name\tTriTrypDB\tcontig\t1\t$ref_size\t.\t+\t.\tName=$ref_name\n";

foreach my $gene (@genes) {
    my %sladds;

    if($$gene->has_tag('SLadd')) {

	my @sladds = $$gene->get_tag_values('SLadd');
#	my $sladdstring = $sladdarray[0];
#	my @sladds = split(/\s+/, $sladdstring);
	
	foreach my $sladd (@sladds) {
#	    print $sladd, "\n"; #DEBUG
	    
	    (!defined($sladds{$sladd})) && ($sladds{$sladd}=0);
	    $sladds{$sladd}++;
	    if( $$gene->strand == -1) {
		$fputrlen{$sladd} = $sladd - ($$gene->end+1) +1 ;
	    } else {
		$fputrlen{$sladd} = ($$gene->start-1) - $sladd + 1 ;
	    }
	}

	# add only one length per unique site - have the individual counts anyway
	$$gene->remove_tag('SLadd');
	foreach my $sladd (keys %sladds) {
	    $$gene->add_tag_value('fpUTRlen', $fputrlen{$sladd});
	    $$gene->add_tag_value('SLadd', $sladd);
	}
    }

    if($$gene->has_tag('SLadd_internal') ) {
	my @sladds = $$gene->get_tag_values('SLadd_internal');
#	my $sladdstring = $sladdarray[0];
#	my @sladds = split(/\s+/, $sladdstring);

	foreach my $sladd (@sladds) {
	    (!defined($sladds{$sladd})) && ($sladds{$sladd}=0);
	    $sladds{$sladd}++;
	}

	# compress the SLadd tag
	$$gene->remove_tag('SLadd_internal');
	foreach my $sladd (keys %sladds) {
#	    $$gene->add_tag_value('fpUTRlen', $fputrlen{$sladd});
	    $$gene->add_tag_value('SLadd_internal', $sladd);
	}
    }

    my $total = 0;
    if ( keys %sladds > 0 ) {
	foreach my $sladd (sort keys %sladds) { 
	    $$gene->add_tag_value('SLaddcounts',$sladds{$sladd});
	    $total+=$sladds{$sladd};
	}
    }

    my $normalised_total = $total / $library_size * 1000000;

    $$gene->add_tag_value('SLadd_totalcount',$total);
#    $$gene->set_attributes('-score'=>$total);
    $$gene->add_tag_value('SLadd_norm_totalcount',sprintf("%.0f", $normalised_total));

    if($$gene->has_tag('score')) {
	$$gene->remove_tag('score');
    }
    $$gene->add_tag_value('score',sprintf("%.0f", $normalised_total));
    
#    print $outputfh $$gene->gff_string,"\n";
    print $outputfh gff3_string($$gene),"\n";
}
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Crunch tags: ";

my $last_pos=-1;
my $tag_count = 0;

my @condensed_tags = ();

my $current_new_tag;

foreach my $tag (@tag_maps) {

    if (( ($$tag->strand == -1) and ($$tag->end != $last_pos)) or (($$tag->strand == 1) and ($$tag->start != $last_pos)) ) {

	if ($last_pos != -1){ 
	    # finalise & print.
	    if($current_new_tag->has_tag('score')) {
		my $total = $current_new_tag->score;	
		$current_new_tag->remove_tag('score');
		my $normalised_score = $total / $library_size * 1000000;
		$current_new_tag->add_tag_value('score',sprintf("%.0f", $normalised_score));
		$current_new_tag->add_tag_value('raw_count',($total));
		$current_new_tag->add_tag_value('normalised_count',sprintf("%.0f", $normalised_score));
	    }
	    print $outputfh gff3_string($current_new_tag),"\n";
	}

	# open new
	$current_new_tag = $$tag;
	if($current_new_tag->has_tag('score')) {
	    $current_new_tag->remove_tag('score');
	}
	if($current_new_tag->has_tag('Readqual')) {
	    $current_new_tag->remove_tag('Readqual');
	}
	$current_new_tag->add_tag_value('score',1);
	$last_pos = ($$tag->strand == -1)?$$tag->end:$$tag->start;

    } else {
	# condense
	
	# increment score for count
#	$current_new_tag->score();
#	$current_new_tag->set_attributes('-score'=>( ));

	if($current_new_tag->has_tag('score')) {
	    $old_score = $current_new_tag->score;
	    $current_new_tag->remove_tag('score');
	    $current_new_tag->add_tag_value('score',($old_score+1));
	}

	# join tags
	foreach my $key ($$tag->get_all_tags) {
	    if ( $key eq "Readseq" ) {
		$current_new_tag->add_tag_value($key,$$tag->get_tag_values($key))
	    }
	}
    }    
}

#and for the last one..
if ($last_pos != -1){ 
    if($current_new_tag->has_tag('score')) {
	my $total = $current_new_tag->score;	
	$current_new_tag->remove_tag('score');
	my $normalised_score = $total / $library_size * 1000000;
	$current_new_tag->add_tag_value('score',sprintf("%.0f", $normalised_score));
	$current_new_tag->add_tag_value('raw_count',($total));
	$current_new_tag->add_tag_value('normalised_count',sprintf("%.0f", $normalised_score));
    }	
    print $outputfh gff3_string($current_new_tag),"\n";
}

$WARNING && print STDERR "OK\n";

exit 0;

sub gff3_string {
    my $feat = shift;
   
    my $strand = ($feat->strand == 1)?'+':'-'; 
    my $phase = '.';
    
    my $attribute_string ="";
   
    foreach my $tag ($feat->get_all_tags) {
	if($tag eq "score") { next; } #supress printing score-tag; has an own column already..

	if($attribute_string ne "") {
	    $attribute_string .= ";";
	}
	#nb if you have multiple freetext tags they need doublequotes, right?
        $values=join(",",$feat->get_tag_values($tag));
	
#	$feat->add_tag_value($key,$tmp_value);
#if($values =~ /[^0-9eE.-]+/
	$attribute_string.="$tag=$values";
    }
    #seq is not an argument.. =)
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
# 		my ($key, $value) = ($feature_string =~ m/([^=\s]+)\s*[=]+\s*(.+)/);
# #print "key: $key is value: $value\n"
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
