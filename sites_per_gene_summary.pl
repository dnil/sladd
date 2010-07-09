#!/usr/bin/perl -w

use POSIX;
my $sum=0;
my $n=0;
my @val=();

#my $libsize = $ARGV[0];

while(<STDIN>) {
    chomp;
    my @r=split(/\t+/);
    my $utrs=0;
    my $count_gene=0; 

    if (@r>1 && $r[3] > 0) {
#		@counts=split(/,/,$r[5]);
	my @UTRlen=split(/,/,$r[7]);

	for (my $i=0; $i<@UTRlen;$i++) { 
	    if($UTRlen[$i]>0 && $UTRlen[$i]<2000) { 
		$utrs++;
		$count_gene=1; 
	    }
	}
    }
    if($count_gene==1) {
	$n++;
	$sum +=$utrs;
	push @val, $utrs;
    }
}
    
if ($n%2==0) { 
    $median = ($val[floor($n/2)-1]+$val[floor($n/2)])/2 
} else { 
    $median = $val[$n/2-1] 
}

print "Sum\tN\tMean\tMedian\n";
print "$sum\t$n\t".sprintf("%.3f",$sum/$n)."\t$median\n";
