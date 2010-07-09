package seqfunk;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(debug fatal warning revcomp iupac_revcomp);

sub debug {
    my $msg = shift;

    print STDERR "DEBUG: $msg\n";
}

sub fatal {
    my $msg = shift;

    print STDERR "FATAL: $msg\n";
    exit;
}

sub warning {
    my $msg = shift;

    print STDERR "WARNING: $msg\n";
}

sub revcomp {
    $_ = shift;
    
    tr/ATUGCatugc/TAACGtaacg/;
    return scalar(reverse($_));
}

sub iupac_revcomp {
    $_ = shift;

    tr/ATUGCRYMKBDVHatugcrymkbdvh/TAACGYRKMACTGtaacgyrkmactg/;
    return scalar(reverse($_));
}


