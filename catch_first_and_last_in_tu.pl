#!/usr/bin/perl -w
#
# Daniel Nilsson, 2009-04-23
#

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

=head1 NAME

catch_first_and_last_in_tu.pl - pick N genes from head and tail of TUs

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@izb.unibe.ch, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2009, 2010 held by Daniel Nilsson. The package is realesed for use under the Perl Artistic License.

=head1 SYNOPSIS

USAGE: C<catch_first_and_last_in_tu.pl -g tu_gff_file [-n pickN(4)]>

=head1 DESCRIPTION

Pick a number of genes from start and end of a TU, using TU gffs from process_tu.pl, keeping
track of TU strand orientation so that "head" genes are always listed
first, and "tail" genes last.

=head1 DEPENDENCIES

BioPerl. process_tu.pl to produce the TU GFF input files.

=cut

my $DEBUG = 0;
my $WARNING = 1;

my $tu_gff_file_name ="";
my $gff_out = "";
my $output_file_name = "";

my $pickN = 4;

while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-g') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-g requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $tu_gff_file_name = $next_arg;
            }
        }

	if ($arg eq '-n') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-n requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $pickN = $next_arg;
            }
        }
	
	# TODO:
	# tail-only flag..
	# head-only flag..

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

if($tu_gff_file_name eq "") {
    print STDERR "No gene gff file name given. Bailing out.\n";
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

my @tus;

$WARNING && print STDERR "Read TU GFF file ($tu_gff_file_name): ";
open(TUGFF,$tu_gff_file_name) or die "Could not open $tu_gff_file_name\n";
while ( my $r = <TUGFF> ) {
    chomp $r;
    my $feat=parse_gff_row($r);
    if ($feat->primary eq "TU") {	
	push @tus, $feat;
    }    
}
close TUGFF;
$WARNING && print STDERR "OK\n";

my $tooshort = 0;

$WARNING && print STDERR "Write output file: ";
while ( my $tu = shift @tus ) {
    
    # genes of tu
    my (@tugenes) = ($tu->get_tag_values('gene_name'));
    my $n = scalar(@tugenes);

    if ( $n < 2*$pickN ) {
	$tooshort++;
	# skip the too short to separate units? let awk do that? nah, easy to forget. do discard them.
	next;
    }

    my @lasttugenes;
    my @firsttugenes;
    if ( ($tu->strand == 1) ) {
	(@lasttugenes)= (@tugenes[($n-$pickN)..($n-1)]);
	(@firsttugenes)= (@tugenes[0..($pickN-1)]);
    } elsif ( $tu->strand == -1 ) {
	(@firsttugenes)= (@tugenes[($n-$pickN)..($n-1)]);
	(@lasttugenes)= (@tugenes[0..($pickN-1)]);
    }

    my ($ref_name) = ($tu->get_tag_values('ref_name'));

    print $outputfh "$ref_name\t\t\t",$tu->start,"\t",$tu->end,"\t",$n,"\t",$tu->strand,"\t",join(",",@firsttugenes),"\t",join(",",@lasttugenes),"\n";
    
}
$WARNING && print STDERR "OK ($tooshort units with <",2*$pickN," genes were ignored).\n";
    
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
    
    # as we don't attach to a reference sequence, keep track anyway..

    $feat->add_tag_value("ref_name",$ref_name);

    return \$feat;
}



