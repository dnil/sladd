#!/usr/bin/perl -w

=head1 NAME

maqmaptxt2gff.pl

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@izb.unibe.ch, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2009, 2010 held by Daniel Nilsson. The package is realesed for use under the Perl Artistic License.

=head1 SYNOPSIS

USAGE: C<<maqmaptxt2.pl -l <libraryname> -f <featurename> -s <sizefilename> >>

=head1 DESCRIPTION

Convert maq mapping txt results to GFF format (version 3). For proper GFF3 output, a valid size file is required.

=cut

my $lib ="maqlibplaceholder";
my $featurename= "maq_match";
my $sizefilename = "";

while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-l') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-l requires an argument, but none given. Bailing out.\n";
                exit 1;
            } else {
                $lib = $next_arg;
            }
        }
	
        if ($arg eq '-f') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-f requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $featurename = $next_arg;
            }
        }

        if ($arg eq '-s') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-s requires an argument, but none given. Bailing out.\n";
                exit 1;
            } else {
                $sizefilename = $next_arg;
            }
        }
    }  
}

my %ref_size;
if($sizefilename ne "") {
    open SIZES, $sizefilename; 

    while(my $r=<SIZES>) {
	chomp $r;
	my @r=split(/\t+/,$r);
	my $ref_name = $r[0];
	$ref_size{$ref_name}= $r[1];
#	print $r[0]," is ",$r[1],"\n";
    }
    close SIZES;
}

print "##gff-version 3\n";
my $last_ref_name= "";

while(my $r=<STDIN>) {
    chomp $r;
    my @r=split(/\t+/, $r);

    my $strand = $r[3];
    my $ref_name = $r[1];
    if ($ref_name ne $last_ref_name) {
	my $size=$ref_size{$ref_name};
	print "##sequence-region $ref_name 1 $size\n";
	print "$ref_name\tTriTrypDB\tcontig\t1\t$size\t.\t+\t.\tName=$ref_name\n";
	$last_ref_name=$ref_name;
    }
    my $read_name = $r[0];
    my $ref_start = $r[2]; # NB BOWTIE:"0-based offset into the reference sequence where leftmost character of the alignment occurs"; appears to be the same in maq? ah, no, 1 based already!
    my $read_seq = $r[14];
    my $read_qual = $r[15];
#    my $mismatch_desc = $r[7];

    $ref_start = $ref_start;
    $ref_end = $ref_start+ length($read_seq)-1;

    if( $strand eq '+' ) {
    } elsif ($strand eq '-')  {
    } else {
	print STDERR "WARNING: weird strand assignment! Column error?\n";
    }

    #gff out

    print "$ref_name\t$lib\t$featurename\t$ref_start\t$ref_end\t.\t$strand\t.\tReadseq=$read_seq;Readqual=$read_qual";
#    defined($mismatch_desc) && print ";Mismatches=$mismatch_desc";
    print "\n";
}
