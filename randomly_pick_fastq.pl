#!/usr/bin/perl -w
# Daniel Nilsson, 09-05-22

my $DEBUG = 0;
my $WARNING = 1;

my $fastafile = "";
my $outputfile = "";
my $target_number = 100000;
my $n_entries = 500000;

while (my $arg = shift @ARGV) {

    if ($arg =~ /^-/) {	
	if ($arg eq '-n') {
	    my $next_arg = shift @ARGV; 
	    if($next_arg eq "") {
		print "-n requires an argument, but non given. Bailing out.\n";
		exit 1;
	    } else {
		$target_number = $next_arg;
	    }
	}
	
	if ($arg eq '-N') {
	    my $next_arg = shift @ARGV; 
	    if($next_arg eq "") {
		print "-N requires an argument, but non given. Bailing out.\n";
		exit 1;
	    } else {
		$n_entries = $next_arg;
# could compute as e g grep -c ^@ $fastafile
	    }
	}

	if ($arg eq '-o') {
	    my $next_arg = shift @ARGV; 
	    if($next_arg eq "") {
		print "-o requires an argument, but non given. Bailing out.\n";
		exit 1;
	    } else {
		$outputfile = $next_arg;
	    }
	}

	if($arg eq '-v') { 
	    # TO BE IMPLEMENTED!
	}

    } else {
	$fastafile = $arg;		
    }
}

if($fastafile eq "") {
    print "No fasta file to search was given. Nothing to do!\n";
    exit 1;
}

# require target > nentries
if ($target_number >= $n_entries or $target_number <= 0 or $n_entries <= 0) {
    print "Oops: target $target_number, entries $n_entries - thats not supported!\n";
    exit 1;
}

open FASTAFILE, $fastafile or die "Could not open $fastafile.\n";

my $outputfh = *STDOUT;
if($outputfile ne "") {
    local *OUTPUT;
    open OUTPUT, ">$outputfile";
    $outputfh = *OUTPUT;
} else {
    $DEBUG && print "Writing fastq out to stdout.\n";
}

$WARNING && print STDERR "Pick random entries: ";

my $init = 0;
my $not_init = 1;
my $select_limit = $target_number;

if($target_number / $n_entries > 0.5) {
    $init = 1;
    $not_init = 0;
    $select_limit = $n_entries - $target_number;
}

my @print_entries = ($init) x $n_entries;
my $selected_number = 0;

while($selected_number < $select_limit ) {
    # will get slow if $target ~ .5 * $n_entries
    my $test = int(rand($n_entries));
    if( $print_entries[$test] == $init ) {
	$print_entries[$test] = $not_init;
	$selected_number++;
    }
}

$WARNING && print STDERR "OK.\n";

$WARNING && print STDERR "Read and write fastq files: ";

my $printentry = 0;
my $current_entry = 0;

while( my $row = <FASTAFILE> ) {

    chomp $row;
    if ( $row =~ m/^\@/ ) {
	$DEBUG && print STDERR "Checking $row.\n";

	$printentry = 0;
	$current_entry++;
	
	if ( $print_entries[$current_entry-1] == 1 ) {	    
	    print $outputfh $row."\n";
	    $printentry = 1;
	}
	
    } elsif ( $printentry ) {
	print $outputfh $row."\n";
    }
}

$WARNING && print STDERR "OK.\n";
