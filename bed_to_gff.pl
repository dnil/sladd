#!/usr/bin/perl -w

my $featurename ="fasterisGDG4";
my $libraryname="fasterisGDG4";

my $sizefilename = "";
while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-l') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-l requires an argument, but non given. Bailing out.\n";
                exit 1;
            } else {
                $libraryname = $next_arg;
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
	my @r=split(/\t+/, $r);
	my $ref_name = $r[0];
	$ref_size{$ref_name}= $r[1];
    }
    close SIZES;
}

print "##gff-version 3\n";
my $last_ref_name= "";
my $note_add = "";

while(my $r=<STDIN>) {
    chomp $r;
    my @r=split(/\t+/, $r);

    if($r =~ m/^browser position/) {
	next;
    }

    if( $r=~ m/track name=/ or $r=~ m/description=/ ) {	
	if($r=~m/track name="([^"]+)"/) {
	    $note_add.=$1." ";
	}
	
	if($r=~m/description="([^"]+)"/) {
	    $note_add.=$1;
	}
	next;
    }
    
    # bed file in
    my $strand = $r[5];
    my $ref_name = $r[0];
    if ($ref_name ne $last_ref_name) {
	my $size=$ref_size{$ref_name};
	print "##sequence-region $ref_name 1 $size\n";
	print "$ref_name\tTriTrypDB\tcontig\t1\t$size\t.\t+\t.\tName=$ref_name\n";
	$last_ref_name=$ref_name;
    }
    my $ref_start = $r[1]; # NB BOWTIE:"0-based offset into the reference sequence where leftmost character of the alignment occurs"

    my $ref_end = $r[2];
    my $count = $r[3];

    if( defined($strand) && $strand eq '+' ) {
    } elsif ( defined($strand) && $strand eq '-')  {
    } else {
	print STDERR "WARNING: weird strand assignment! Column error? $r\n"; 
	next;
    }

    #gff out
    print "$ref_name\t$libraryname\t$featurename\t$ref_start\t$ref_end\t$count\t$strand\t.\tNote=Count $count tags.$note_add\n";
}
