#!/usr/bin/perl -w
# 
# Daniel Nilsson

my $inseq = 0; 
my $seq = "";
my $name = "";

while ($l = <STDIN>) {
    chomp $l;
    
    if ($inseq) {

	if ($l =~ m/^>(\S+)/) { 
	    $inseq=1; 

	    print $name, "\t", length($seq), "\n";
	    
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
    print $name, "\t", length($seq), "\n";
}
