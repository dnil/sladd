#!/usr/bin/perl -w

my $fastafile = $ARGV[0];
my $phobius_long = $ARGV[1];

open FASTA, "<$fastafile";

my $currentseq = "";
my %seq;
$first=1;

while(my $r = <FASTA>) {
    chomp $r;
    
    if($r=~/\>(.+)/) {
	if(!$first) {
	    $seq{$currenthead} = $currentseq;
	} 
	$currenthead=$1;
	$currentseq ="";
	$first=0;
    } else {
	$currentseq .= $r;
    }           
}

$seq{$currenthead} = $currentseq;

open PHOB, "<".$phobius_long;
while(my $l = <PHOB>) {    
    chomp $l;

    if($l =~/ID\s+(\S+)/ ) {
	$id = $1;
    }
    if($l =~/FT\s+SIGNAL\s+(\d+)\s+(\d+)/) {
	$signal_start=$1;
	$signal_end =$2;
	
	my @seqs = split(/ */,$seq{$id});
	for (my $i = $signal_start-1; $i < $signal_end; $i++) {
	    $seqs[$i] = "X";
	}
	
	$seq{$id} = join("",@seqs);
    }
}

close PHOB;

foreach my $seqname (keys %seq) {
    print ">$seqname\n".$seq{$seqname}."\n";
}
