#!/usr/bin/perl -w

my $bowtielib ="bowtielibplaceholder";
my $featurename= "nucleotide_match";

while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-l') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-l requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $bowtielib = $next_arg;
            }
        }

        if ($arg eq '-s') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-s requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $sizefilename = $next_arg;
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

    my $strand = $r[1];
    my $ref_name = $r[2];
    if ($ref_name ne $last_ref_name) {
	my $size=$ref_size{$ref_name};
	print "##sequence-region $ref_name 1 $size\n";
	print "$ref_name\tTriTrypDB\tcontig\t1\t$size\t.\t+\t.\tName=$ref_name\n";
	$last_ref_name=$ref_name;
    }
    my $read_name = $r[0];
    my $ref_start = $r[3]; # NB BOWTIE:"0-based offset into the reference sequence where leftmost character of the alignment occurs"
    my $read_seq = $r[4];
    my $read_qual = $r[5];
    my $mismatch_desc = $r[7];

    $ref_start = $ref_start+1;
    $ref_end = $ref_start+ length($read_seq)-1;

    if( $strand eq '+' ) {
    } elsif ($strand eq '-')  {
    } else {
	print "WARNING: weird strand assignment! Column error?\n";
    }

    #gff out

    print "$ref_name\t$bowtielib\t$featurename\t$ref_start\t$ref_end\t.\t$strand\t.\tReadseq=$read_seq;Readqual=$read_qual";
#    defined($mismatch_desc) && print ";Mismatches=$mismatch_desc";
    print "\n";
}
