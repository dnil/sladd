#!/usr/bin/perl -w
# 
# Daniel Nilsson
# 2010-01-11
# 
# get internally spliced peptides
#

my $DEBUG = 0;

if (@ARGV == 0) {
    print STDERR "Usage: get_internally_spliced_pep.pl pepfile altspltab > out_short_pep_fasta\n";
    exit;
}

# read in and store internal splice sites (in pep lengths..)

my $pepfile = $ARGV[0];
my $altspltab = $ARGV[1];

my $pep_length_cutoff = 80;

my $inseq = 0; 
my $seq = "";
my $name = "";

# read peptides

my %seq;

open(PEPFILE, "<".$pepfile);

while (my $l = <PEPFILE>) {
    chomp $l;
    
    if ($inseq) {

	if ($l =~ m/^>(\S+)/) { 
	    $inseq=1; 

#	    print $name, "\t", length($seq), "\n";

	    $seq{$name}=$seq;
	    
	    $name = $1;
	    $seq = "";
	} else {
	    $l =~ s/\s+//g;
	    $seq .= $l;	
	}
    } elsif ($l =~ m/^>(\S+)/) { 
	# first time around..
	$inseq=1; 
	
	$name = $1;
	$seq = "";
    }
}

if($inseq==1) { 
#    print $name, "\t", length($seq), "\n";
    $seq{$name}=$seq;
}


# read alt splice site tab file, and also process for internally spliced complete peptides..

open(ALTSPLTAB, "<".$altspltab);

while(my $l = <ALTSPLTAB>) {

    if ($l=~ /^[=]+/ or $l=~/^[#]+/) {
	# comment/info line -- ignore
	next;
    } 
   
    my @l= split(/\t/,$l);
    my $id = $l[0];
    my $counts = $l[5];
    my $splicesites = $l[6];

    $counts=~ s/"//g;
    my @counts = split(/,/,$counts);

    $splicesites=~ s/"//g;
    my @splicesites = split(/,/,$splicesites);
    my @internal_splicesites =  grep { $_< 0 } @splicesites;

    if( @internal_splicesites == 0 ) {
	print STDERR "WARNING: no internal splicesites found for $id. \n";
	next;
    }

    # and check for M-starts..

    if ( ! defined($seq{$id}) ) {
	my $origid = $id;
	$id =~ s/-1//;
	
	if (! defined($seq{$id}) ){
	    $nonpsuid= $id;
	    $id = "psu|$nonpsuid";
	    if( ! defined($seq{$id}) ){
		print STDERR "WARNING: could not find sequence for id $id (nor $nonpsuid or $origid) in peptide file $pepfile.\n";
		next;
	    } else {
		# found it - reintroduce the non-psu (~tabfile) id..
		$seq{$nonpsuid} = delete $seq{$id};	    
		$id = $nonpsuid;
	    }
	} 
    }



    use POSIX;

    my $totallength=length($seq{$id});    
    my @internal_min_truncations = map { ceil($totallength - $_ )/3 } @internal_splicesites;
    
#    my $max_int_len = $internal_min_truncations[0];

    my $last_pep_length=$totallength;

    foreach my $imt (@internal_min_truncations) {
	$DEBUG && print "DEBUG: trying new imt $imt for $id\n";

	# check for stop codons..
	if ( $seq{$id} =~ m/^([^*]+)\*+/ ) {
	    # and truncate to first stop if any.
	    $seq{$id} = $1;
	    $totallength=length($seq{$id});
	}
      
	while( $seq{$id} =~ m/(M+[^M]+)/g ) {
	    
	    my $current_pep_start = $1;
	    my $current_pep_end_pos = pos $seq{$id};
	    my $current_pep_end = substr($seq{$id},$current_pep_end_pos); # the rest of it..
	    
	    my $current_pep= $current_pep_start.$current_pep_end;

	    $DEBUG && print "DEBUG: trying pep for $id, length ", length($current_pep), "\n";

	    if(length($totallength) == length($current_pep)) {
		# this is the full length
	    } elsif (length($current_pep) <= $imt && length($current_pep) != $last_pep_length)  {

		$last_pep_length = length($current_pep);

		if ($last_pep_length < $pep_length_cutoff) {
		    $DEBUG && print STDERR "DEBUG: Too short with $last_pep_length residues.\n";
		} else {
		    print ">".$id."_".$last_pep_length."_int\n";
		    # prettyprint
		    
		    my @pep = split(/ */, $current_pep);		    
		    while ( @curr = splice @pep, 0, (80>=scalar(@curr))?80:scalar(@curr)) {
			print join("",@curr),"\n";
		    }
		}
		
		# skip to next imt? yes, probably.
		last; # print no shorter peptides for this splicesite..
	    }
	}
#	if ( $imt > $max_int_len ) {
#	    $max_int_len = $imt;
#	}

    }
}

close(ALTSPLTAB);

