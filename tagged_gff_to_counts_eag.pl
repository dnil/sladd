#!/usr/bin/perl -w -I /Users/daniel/install/splicemodel
#
# Daniel Nilsson, 2009-04-23
#

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

use seqfunk;

my $DEBUG = 0;
my $WARNING = 1;

my $tagged_gff_file_name ="";
my $gene_feature_name = "gene";
my $tag_map_feature_name = "nucleotide_match";
my $tag_map_file_name="";
my $sequence_file_name = "";
my $gff_out = "";
my $output_file_name = "";

my $utrlencutoff = 5000;

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

        if ($arg eq '-c') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-c requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $gene_feature_name = $next_arg;
            }
        }

        if ($arg eq '-N') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-s requires an argument, but non given. Bailing out.\n";
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
my @tag_maps;

foreach my $feat ($seq->all_SeqFeatures()) {
#    if( $feat->primary_tag eq $tag_map_feature_name) {
#	push @tag_maps, \$feat;
#    } 
    if( $feat->primary_tag eq $gene_feature_name) {
	push @genes, \$feat;
    }
}
$WARNING && print STDERR "OK\n";

$WARNING && print STDERR "Agglomerate counts: ";
foreach my $gene (@genes) {
    my %sladds;
    my %fputrlen;
    my %altstarts;
    my %eags; 

    my $shortest_fputr=-1;
    my $longest_fputr=-1;
    my $major_fputr=-1;
    my $major_fputr_count=-1;
    my $major_sladd=-1;

    if($$gene->has_tag('SLadd')) {
	my @sladds = $$gene->get_tag_values('SLadd');
#	if (@sladdarray > 1) {
#	    print "OOPS: sladdarray is now working ok ".scalar(@sladdarray)."\n";
#	}
#	my $sladdstring = $sladdarray[0];
#	my @sladds = split(/\s+/, $sladdstring);

	foreach my $sladd (@sladds) {
#	    print $sladd, "\n"; #DEBUG
	    
	    (!defined($sladds{$sladd})) && ($sladds{$sladd}=0);
	    $sladds{$sladd}++;
	    if( $$gene->strand == -1) {
		$fputrlen{$sladd} = $sladd - ($$gene->end+1) +1 ;
	    } else {
		$fputrlen{$sladd} = ($$gene->start - 1) - $sladd + 1 ;
	    }

	    $$gene->add_tag_value('fpUTRlen', $fputrlen{$sladd});

	    if($shortest_fputr == -1) {
		$shortest_fputr = $fputrlen{$sladd};
	    } elsif($fputrlen{$sladd} < $shortest_fputr ) {
		$shortest_fputr =$fputrlen{$sladd};
	    }

	    if($longest_fputr == -1) {
		$longest_fputr = $fputrlen{$sladd};
		$longest_sladd = $sladd;
	    } elsif($fputrlen{$sladd} > $longest_fputr && $fputrlen{$sladd} < $utrlencutoff) {
		$longest_fputr =$fputrlen{$sladd};
		$longest_sladd = $sladd;
	    }
 
	    if($major_fputr == -1) {
		$major_sladd = $sladd;
		$major_fputr = $fputrlen{$sladd};
		$major_fputr_count = $sladds{$sladd};
	    } elsif($sladds{$sladd} > $major_fputr_count && $fputrlen{$sladd} < $utrlencutoff ) {
		$major_sladd = $sladd;
		$major_fputr = $fputrlen{$sladd};
		$major_fputr_count = $sladds{$sladd};
	    }
	    
	    if($fputrlen{$sladd}<$utrlencutoff) {
		$eags{$sladd} = check_eag($sladd,$gene,$seq);
	    } else {
		$eags{$sladd} = -1;
	    }
	}
    }

    if($$gene->has_tag('SLadd_internal') ) {

	my @sladds = $$gene->get_tag_values('SLadd_internal');
#	my $sladdstring = $sladdarray[0];
#	my @sladds = split(/\s+/, $sladdstring);

	foreach my $sladd (@sladds) {
	    (!defined($sladds{$sladd})) && ($sladds{$sladd}=0);
	    $sladds{$sladd}++;

	    #add negative values for utrlen
	    if( $$gene->strand == -1) {
		$fputrlen{$sladd} = $sladd - ($$gene->end+1) +1 ;
	    } else {
		$fputrlen{$sladd} = ($$gene->start - 1) - $sladd + 1 ;
	    }

	    if($longest_fputr == -1) {
		$longest_sladd = $sladd;
		$longest_fputr = $fputrlen{$sladd};
	    } elsif($fputrlen{$sladd} > $longest_fputr && $fputrlen{$sladd} < $utrlencutoff) {
		$longest_sladd = $sladd;
		$longest_fputr =$fputrlen{$sladd};
	    }

	    #don't count shortest for internal splice sites.

	    if($major_fputr == -1) {
		$major_sladd = $sladd;
		$major_fputr = $fputrlen{$sladd};
		$major_fputr_count = $sladds{$sladd};
	    } elsif($sladds{$sladd} > $major_fputr_count) {
		$major_sladd = $sladd;
		$major_fputr = $fputrlen{$sladd};
		$major_fputr_count = $sladds{$sladd};
	    }
	    
	    $eags{$sladd}= -1;
	}
    }

    my $total = 0;
    if ( keys %sladds > 0 ) {
	foreach my $sladd (sort keys %sladds) {
	    $$gene->add_tag_value('SLaddcounts',$sladds{$sladd});
	    $total+=$sladds{$sladd};
	}
    }
    $$gene->add_tag_value('SLadd_totalcount',$total);

    my $name=$$gene->start." ".$$gene->end;
    if($$gene->has_tag('Name')) {
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
    
    my $normalised_total = $total / $library_size * 1000000;

    print $outputfh $name,"\t",$atg_pos,"\t",$strand,"\t",$total,"\t",sprintf("%.0f", $normalised_total);
    if($total > 0) {
#	my @addcountsarr =$$gene->get_tag_values('SLaddcounts');
	#sorted the same way? better be sure.
	my @addcountsarr = ();
	my @sladds =();
	my @fputrarr = ();
	my @eags = ();
	my @altstarts= ();

	foreach my $sladd (sort keys %sladds) { 
	    push @sladds, $sladd;
	    push @addcountsarr, $sladds{$sladd};
	    push @fputrarr, $fputrlen{$sladd};
	    push @eags, $eags{$sladd};
	    push @altstarts, $altstarts{$sladd};
	}

	my $major_eags = $eags{$major_sladd};

	print $outputfh "\t",join(",",@addcountsarr), "\t", join(",",@sladds),"\t", join(",",@fputrarr),"\t", join(",",@eags),"\t", join(",",@altstarts),"\t", $major_sladd,"\t", $major_fputr,"\t",$major_fputr_count,"\t",$major_eags, "\t",$shortest_fputr,"\t",$longest_fputr;
    } else {
	print $outputfh "\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.";
    }

    print $outputfh "\t",$description,"\n";
}
$WARNING && print STDERR "OK\n";

exit 0;

sub check_eag {  
    my $sladd = shift;
    my $gene = shift; #ref
    my $seq = shift; #ref

    my $atg_pos = ($$gene->strand == 1)?$$gene->start:$$gene->end;

    my $direction=$$gene->strand;
    if($direction == 1) {
	$first = $sladd -1;
	$last = $atg_pos - 1;
    } else {
	$first = $atg_pos + 1;
	$last = $sladd +1;
    }
    
    # longest <5kb?
    
    my $fiveputr = lc $seq->subseq($first,$last);
    ($$gene->strand == 1 ) || ($fiveputr = revcomp($fiveputr));
        
    # disable dupechecknig, check anew for each UTR.

    my $nags = 0;     
   
    if($fiveputr =~ m/ag/) {
	
	while( $fiveputr =~ m/(.*?)ag/g ) {
	    $nags++;
	}
    }

    return ($nags);
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
