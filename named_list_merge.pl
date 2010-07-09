#!/usr/bin/perl -w  

my $row1; 
my $row2; 

while(1) { 
    if(!defined($row1)) {
	$row1 = <STDIN>; 
	chomp $row1; 
    } 

    $row2 = <STDIN>; 
    if (!defined($row2)) { 
	print $row1."\n"; 
	exit;
    } 
    chomp $row2; 
    
    ($systname1) = ($row1 =~ m/^(\S+)\s+/); 
    ($systname2) = ($row2 =~ m/^(\S+)\s+/); 

    $row1 =~ s/\s+/\t/g;
    $row2 =~ s/\s+/\t/g;
    
    if ($systname1 eq $systname2) {
	print $row1."\t".$row2."\n"; 
	undef($row1);
    } else { 
	print $row1."\n"; 
	$row1 = $row2; 
	undef($row2);
    }
}
