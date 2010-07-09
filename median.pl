#!/usr/bin/perl -w 

use POSIX; 

my $sum=0; 
my $n=0;
my @val=();

while(my $r=<STDIN>) {
    chomp $r; 
    $sum+=$r; 
    $n++; 
    push @val,$r;
} 

my $median;
if ($n%2==0){ 
    $median = ($val[floor($n/2)-1]+$val[floor($n/2)])/2;
} else { 
    $median = $val[$n/2-1];
} 

print "mean ".sprintf("%.3f",$sum/$n)." median $median\n";
