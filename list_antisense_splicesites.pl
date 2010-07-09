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
    }  
}

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
@genes = sort { ($$a->start <=> $$b->start) } @genes;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Sort tags: ";
@tag_maps = sort { ($$a->strand == 1 && $$a->start <=> $$b->start) || ($$a->strand == -1 && $$b->end <=> $$a->end)} @tag_maps;
$WARNING && print STDERR "OK\n";

my $current_tag = 0;

$WARNING && print STDERR "Determine strandedness: ";

my @gene_strandedness = split(/ */,'0'x$ref_size);

my $ngenes = scalar(@genes);
my $ntags = scalar(@tag_maps);

#$DEBUG && print "N genes = $ngenes, N tags = $ntags.\n";

for (my $i = 0; $i < $ngenes ; $i++) {

    my $current_gene = $genes[$i];
    my $previous_gene;
    if ($i > 0) {
	$previous_gene = $genes[$i-1];
    }    

    $DEBUG && print "DEBUG: Current gene (nr $i) ".($$current_gene->gff_string)."\n";

    if ($i != 0) {
	# not first gene

	my $strand = 0;
	if ($$previous_gene->strand == $$current_gene->strand) {
	    $strand = $$current_gene->strand;
	}

	for (my $j = $$previous_gene->end; $j < $$current_gene->start; $j++) {
	    $gene_strandedness[$j] = $strand;
	}
    } elsif ($i == 0 ) {
	#first gene
	for (my $j = 0; $j < $$current_gene->start; $j++) {
	    $gene_strandedness[$j] = 0;
	}
    }

    if ($i == $ngenes-1) {
	#last gene
	for (my $j = $$current_gene->end; $j < $ref_size; $j++) {
	    $gene_strandedness[$j] = 0;
	}
    }

    #anyway, always mark current gene
    for (my $j = $$current_gene->start; $j < $$current_gene->end; $j++) {	
	$gene_strandedness[$j] = $$current_gene->strand;       
    }
   
}    
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Check tag strandedness: ";

my $tags_in_gray = 0;
my $tags_in_sense = 0;
my $tags_against_sense = 0;

for (my $i = 0; $i < $ntags ; $i++) {

    my $current_tag_map = $tag_maps[$i];
        
    my $tag_strand = $$current_tag_map->strand;

    # sense defaults to true
    my $tag_sense = 1;

    for (my $j = $$current_tag_map->start; $j < $$current_tag_map->end; $j++) {
	
        #is set to gray if any overlapped pos is gray
	if( $gene_strandedness[$j] == 0 ) {
	    $tag_sense = 0;
	}

	# and set to antisense if antisense anywhere, unless gray 
	if ($tag_sense != 0 && $tag_strand != $gene_strandedness[$j]) {
	    $tag_sense = -1;
	}
    }

    if($tag_sense == 1) {	
#	$$current_tag_map->add_tag_value('strandedness', 'sense');
	$tags_in_sense++;
    } else {
	if($tag_sense == 0) {
	    $$current_tag_map->add_tag_value('strandedness', 'unclear');
	    $tags_in_gray++;
	} elsif ($tag_sense == -1) {
	    $$current_tag_map->add_tag_value('strandedness', 'antisense');
	    $tags_against_sense++;
	}
	
	print $outputfh gff3_string($$current_tag_map),"\n";
    }
}
$WARNING && print STDERR "OK\n";

print "In summary, $tags_in_sense tags in sense, $tags_in_gray tags in grayzones and $tags_against_sense tags map anti-sense.\n";

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
		$DEBUG && print "key: $key is value: $value\n";

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
