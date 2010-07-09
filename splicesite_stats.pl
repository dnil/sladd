#!/usr/bin/perl -w

my $DEBUG=0;

my %accdinuc; 
my $sum=0; 

my $max_poly_py_mismatches = 1;
my $min_poly_py_len = 10;
my $min_qualified_py_run = 8;

my $py=0;
my $nopy=0;
my @pylen;
my @pyagdist;
my %nf; # hash of arrays of nfs

my @nts= ('C', 'T', 'A', 'G');

while (my $str=<STDIN>) { 
    chomp $str;
    $str = uc $str;

    my $dinuc=substr $str,-8,2;
    $accdinuc{$dinuc}++;
    $sum++;
    
    my $thisseq = substr $str, 0, -7;

    my $subseqlen = length $thisseq;

   # pick the closest splicesite fulfilling our criteria and save data for it
    my $shortestdist=$subseqlen;
    my $bestmatch="";
    my $bestpylen=0;

    while($thisseq =~ /([TCYN]+)([AGRNX]?[TCYNX]*){0,$max_poly_py_mismatches}([TCYNX]+)/g) {
	my $pystart = $-[0]+1;
	my $pyend = $+[0];
	my $pylen = $pyend - $pystart + 1;

	(length($2) == 0 && (length($1)+length($3) >= $min_qualified_py_run)) 
	    || ( $pylen >= $min_poly_py_len)
	    || next;

	my $pyagdist = $subseqlen - $pyend + 1;
	if ($pyagdist<$shortestdist) {
	    $bestmatch= $&;
	    $shortestdist=$pyagdist;
	    $bestpylen=$pylen;
	}
    }
    
    if ( $shortestdist != $subseqlen ) {
	# pyagdist fig
	push @pyagdist, $shortestdist;
	push @pylen, $bestpylen;
	$py++;

	# check for AG on that stretch.. 

	my @bestmatch=split(/ */,$bestmatch);
	my %nucfreq;
	map { $nucfreq{$_}++ } @bestmatch;
	
	my %nucfreqp;

	map { defined($nucfreq{$_}) || ($nucfreq{$_}=0); $nucfreqp{$_}=$nucfreq{$_}/$bestpylen; push @{$nf{$_}}, $nucfreqp{$_} } @nts;

	if ($DEBUG == 1) {
	    my @nucfstrs=();
	
	    for my $nt (keys %nucfreqp) {
#	    defined($nucfreqp{$nt}) || ($nucfreqp{$nt} = 0);
		push @nucfstrs, $nt."=".sprintf("%.2f",$nucfreqp{$nt});
	    }
	    print $bestmatch, "\t", $shortestdist,"\t", $bestpylen, "\t",join(",",@nucfstrs), "\n";
	}
    } else {
	$nopy++;
    }

#	my %dinuc
#	my %trinuc
} 

print "*** Found $py closest sites, and $nopy entries without any criteria fulfilling sites (max_poly_py_mismatches=$max_poly_py_mismatches, min_poly_py_len=$min_poly_py_len, min_qualified_py_run=$min_qualified_py_run).\n";

#print join (" - ",@{$nf{'C'}}),"\n";
#print join (" - ",@{$nf{'T'}}),"\n";
#print join (" - ",@{$nf{'G'}}),"\n";
#print join (" - ",@{$nf{'A'}}),"\n";

my @nfstrs=();
for my $nt (@nts) {
    push @nfstrs, (map { sprintf("%.2f",$_); } central(@{$nf{$nt}})) ;
}

print join("\t\t\t","pPyAGd","pPyLen",@nts),"\n";

print join("\t",("median","mean","sd")x6 ),"\n";

print join("\t", (map { sprintf("%.2f",$_); } central(@pyagdist)), (map { sprintf("%.2f",$_); } central(@pylen)), @nfstrs ),"\n";

print "Dinuc\tNum\tFraction\n";

for my $dinuc (sort { $accdinuc{$b}<=>$accdinuc{$a}} keys %accdinuc) { 
    print $dinuc,"\t", $accdinuc{$dinuc}, "\t", sprintf("%.2f",$accdinuc{$dinuc}/$sum), "\n";
}

#
#sub central_to_string {
#    $integ =shift; # ahem, median can also be a float.. 
#    $floater =shift;
#
#    return ("$integ", sprintf("%.2f",$floater));
#}

sub central {
    use POSIX; 

    my @val = @_;

    my $sum=0;
    map { $sum+=$_ } @val;

    my $n = @val;

    my $mean=$sum/$n;
    my $median;
    
    if (($n % 2) == 0) {
	$median = ($val[floor($n/2)-1]+$val[floor($n/2)])/2;
    } else {
	$median = $val[$n/2-1];
    }

    my $squaredev=0;
    map { $squaredev += ($_ - $mean)**2} @val;

    my $variance = $squaredev / $n;
	
    my $sd=sqrt($variance);

    return ($median,$mean,$sd);
}



