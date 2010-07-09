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
    parse_gff_row($r);
}
close GENEGFF;
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Retrieve feature references: ";
my @genes;

foreach my $feat ($seq->all_SeqFeatures()) {
    if( $feat->primary_tag eq "mRNA" || $feat->primary_tag eq "rRNA" || $feat->primary_tag eq "snRNA"
	|| $feat->primary_tag eq "tRNA" || $feat->primary_tag eq "transcript") {
	push @genes, \$feat;
    } 
}
$WARNING && print STDERR "OK\n";

print $outputfh "##gff-version 3\n";
my $ref_name = $seq->display_id;
my $ref_size = $seq->length;
print $outputfh "##sequence-region $ref_name 1 $ref_size\n";
print $outputfh "$ref_name\tTriTrypDB\tcontig\t1\t$ref_size\t.\t+\t.\tName=$ref_name\n";

if (@genes == 0) {
    print STDERR "Number of genes is 0. Error in input? Bailing out.\n";
    print "Number of genes is 0 .\n";   
    exit 1;
}

$WARNING && print STDERR "Sort genes: ";
@genes = sort { $$a->start <=> $$b->start } @genes;
$WARNING && print STDERR "OK\n";

my $ngenes = scalar(@genes);

$WARNING && print STDERR "Add TUs: ";

my $tu_open=0;
my $tu_type="";

my $tus=0;
   
my $tu_start;
my $tu_end;
my $tu_direction;

my $tu_genes=0;
my @tu_genes;
    
my $flag_tu_VSG =0;
my $flag_tu_ESA =0;
my $flag_tu_rRNA =0;	
my $flag_tu_snRNA=0;
my $flag_tu_tRNA =0;
my $flag_tu_PARP =0;
my $flag_tu_slRNA = 0;	

my $current_gene;
my $current_description ="";

sub finalise_tu;
sub start_new_tu;
sub add_gene_to_tu;

for (my $i = 0; $i < $ngenes ; $i++) {

    $current_gene = $genes[$i];

    my $current_strand = $$current_gene->strand;

    if ($$current_gene->has_tag('description') ) {
	($current_description)= ($$current_gene->get_tag_values('description'));
	$DEBUG && print STDERR "DEBUG: got description for gene $i: $current_description.\n";
    } else {
	$current_description="";
    }
   

#    $DEBUG && print STDERR "DEBUG: Current gene (nr $i) ".($$current_gene->gff_string)."\n";

    if ( $tu_open == 0 ) {
	start_new_tu;  # i == 0 or after PAG
    } elsif( $current_description =~ m/PAG/ or $current_description =~ m/procyclin\-associated\+gene/i ) {
	$flag_tu_PARP=1;
	add_gene_to_tu;
	finalise_tu;
    } elsif ( $flag_tu_PARP == 0 and ($current_description =~ m/procyclin/i or $current_description =~ m/PARP/) ) {
	# always start a new unit, even if flag_PARP == 1?
	finalise_tu;
	start_new_tu;
	$flag_tu_PARP=1;
    } elsif ( $tu_type ne "mRNA" and $$current_gene->primary_tag eq "mRNA") {	
	finalise_tu;
	start_new_tu;
    } elsif ($tu_type eq "mRNA" and $$current_gene->primary_tag ne "mRNA") {       
	finalise_tu;
	start_new_tu;
    } elsif ($$current_gene->strand ne $tu_direction) {
	finalise_tu;
	start_new_tu;
    } else {
	add_gene_to_tu;
    }

    if ( $current_description =~ m/VSG/i or $current_description =~ m/variant\+surface\+glycoprotein/i) {
	$flag_tu_VSG = 1;
    }

    if ( $current_description =~ m/\bESAG/i or $current_description =~ m/expression\+site\-associated/i) {
	$flag_tu_ESA = 1;
    }

    if ($i == $ngenes-1) {
	# last gene..
	if($tu_open == 1) {
	    finalise_tu;
	}
    }
}

$WARNING && print STDERR "OK\n";

my @tu;

$WARNING && print STDERR "Print TUs: ";
foreach my $feat ($seq->all_SeqFeatures()) {
    if($feat->primary_tag eq 'TU') {
	push @tu, \$feat;
	print $outputfh gff3_string($feat),"\n";
    }

}
$WARNING && print STDERR "OK\n";

# which were the above - lists?

print "There were ".scalar(@tu)." TUs for ".scalar(@genes)." genes.\n";
exit 0;

sub add_gene_to_tu {
    
    $tu_end = $$current_gene->end;

    my $gene_name;
    if($$current_gene->has_tag('Name')) {
	($gene_name)=($$current_gene->get_tag_values('Name'));
	if ($gene_name eq "cds") {
	    if($$current_gene->has_tag('ID')) {
		($gene_name)=($$current_gene->get_tag_values('ID'));
	    }
	}
    }
    push @tu_genes,$gene_name;
    
    $tu_genes++;
}

sub finalise_tu {
    # finalise current tu
    
    $feature_source = "processTU";
    $feature_name = "TU";
    
    if($tu_open == 0) {
	print STDERR "WARNING: RULES INCOMPLETE? Attempt to finalise an already closed tu ($tus).";
	return;
    }
	
    my $tu = new Bio::SeqFeature::Generic ( -start => $tu_start, -end => $tu_end, -strand => $tu_direction,
					    -primary=>$feature_name,
					    -source=>$feature_source );
    
    if($flag_tu_VSG == 1) {
	$tu->add_tag_value('VSG',1);
    }

    if($flag_tu_rRNA == 1) {
	$tu->add_tag_value('rRNA',1);
    }

    if($flag_tu_PARP == 1) {
	$tu->add_tag_value('PARP',1);
    }

    if($flag_tu_ESA == 1) {
	$tu->add_tag_value('ESAG',1);
    }

    if($flag_tu_tRNA == 1) {
	$tu->add_tag_value('tRNA',1);
    }

    if($flag_tu_snRNA == 1) {
	$tu->add_tag_value('snRNA',1);
    }

    if($flag_tu_slRNA == 1) {
	$tu->add_tag_value('slRNA',1);
    }

    $tu->add_tag_value('genes',$tu_genes);
    $tu->add_tag_value('type',$tu_type);
    $tu->add_tag_value('TU_id',$ref_name."_".$tus);

    for my $gene_name (@tu_genes) {
	$tu->add_tag_value('gene_name',$gene_name);
    }

    $seq->add_SeqFeature($tu);
    $tu_open=0;
}

sub start_new_tu {
    $tus++;
   
    $tu_start = $$current_gene->start;
    $tu_end=$$current_gene->end;
    $tu_direction=$$current_gene->strand;

    $tu_genes=1;
    
    $flag_tu_VSG =0;
    $flag_tu_ESA =0;
    $flag_tu_rRNA =0;	
    $flag_tu_snRNA=0;
    $flag_tu_slRNA=0;
    $flag_tu_tRNA =0;
    $flag_tu_PARP =0;

    my $gene_name;
    if($$current_gene->has_tag('Name')) {
	($gene_name)=($$current_gene->get_tag_values('Name'));
	if ($gene_name eq "cds") {
	    if($$current_gene->has_tag('ID')) {
		($gene_name)=($$current_gene->get_tag_values('ID'));
	    }
	}
    }
    @tu_genes = ($gene_name);

    $tu_type = $$current_gene->primary_tag;
    
    if ( $tu_type eq "rRNA" ) {
	$flag_tu_rRNA = 1;	    
    } elsif ( $tu_type eq "tRNA" ) {
	# tRNA may have promoter and terminator elements
	$flag_tu_tRNA = 1;	    
    } elsif ( $tu_type eq "snRNA" ) {
	$flag_tu_snRNA = 1;
    } elsif ( $tu_type eq "transcript" ) {
	# spliced leader array
	$flag_tu_slRNA = 1;	
    }

    $tu_open=1;
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
    
    my $score = '.';
    if ($feat->has_tag('genes') ) {
	($score)= ($feat->get_tag_values('genes'));
    } 

    $gff3_string = join("\t", $seq->display_id, $feat->source_tag, $feat->primary_tag,$feat->start,$feat->end,$score,$strand,$phase,$attribute_string);
    return $gff3_string;
}

sub parse_gff_row {
    my $r = shift;

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
#	    $DEBUG && print "feature desc has \"-enclosed value.\n";
	} else {
	    @feature_strings = split(/[;]+/, $feature_desc); # q&d, could have problems with ; in ""..
	}
#	$DEBUG && print $feature_desc,"\n";
	foreach my $feature_string (@feature_strings) {
	    my ($key, $value);
	    if( $feature_string =~ /^\s*([^=\s]+)\s*$/ ) { # boolean tags.. like "pseudo;"
		$key = $1; 
		$value = 1;
#		$DEBUG && print "boolean tag.\n";
	    } else {
		($key, $value) = ($feature_string =~ m/([^=]+)\s*=\s*(.+)/);
	    }				
#	    $DEBUG && print "key: $key is value: $value\n";
	    
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



