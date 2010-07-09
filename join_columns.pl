#!/usr/bin/perl -w

# join_columns <files>

my $display_header = 1;

my @infiles;
while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-H') {
                $display_header = 0;
	}	
    } else {
	push @infiles,$arg;
    }
}

my %irows;
my @infileorder;

# read all files
my $n_inrows =-1;
while(my $infile = shift @infiles) {
    open IN,$infile;
    
    push @infileorder, $infile;

    while(my $row=<IN>) {
	chomp $row;
	push @{$irow{$infile}}, $row;
    }
 
    if( $n_inrows != -1 ) {
	# check against last file
	if ( $n_inrows != @{$irow{$infile}} ) {
	    print "Ooops: number of rows differ in $infile - ".scalar(@{$irow{$infile}})." rather than previously seen $n_inrows.\n";
	}
    }

    $n_inrows = @{$irow{$infile}};
    
    close IN;
}

# print header (filenames)
if( $display_header == 1 ) {
    print join("\t", @infileorder), "\n";
}

# print columns, assume same numbering
for (my $i=0; $i<$n_inrows; $i++) {

    my $first =1;
    foreach my $infile (@infileorder) {
	print "\t" unless $first;
	print ${irow{$infile}}[$i];	
	$first=0;
    }
    print "\n";    
}
