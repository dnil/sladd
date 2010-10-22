package seqfunk;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(debug fatal warning revcomp iupac_revcomp);

=head1 NAME
    
seqfunk - some small functions for sequnence handling

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2000-2010 held by Daniel Nilsson. The package is realesed for use under the Perl Artistic License.

=head1 SYNOPSIS
    
use seqfunk

Exports a few functions.

=head1 FUNCTIONS

=over 4

=cut

=item debug([message])

Print debug message to STDERR.

=cut 

sub debug {
    my $msg = shift;

    print STDERR "DEBUG: $msg\n";
}

=item fatal([message])

Print fatal message to STDERR and exit.

=cut 

sub fatal {
    my $msg = shift;

    print STDERR "FATAL: $msg\n";
    exit;
}

=item warning([message])

Print warning message to STDERR.

=cut

sub warning {
    my $msg = shift;

    print STDERR "WARNING: $msg\n";
}

=item revcomp(nt_sequence_string)

Reverse and complement nt sequence string (ATUGC, atugc).

=cut

sub revcomp {
    $_ = shift;
    
    tr/ATUGCatugc/TAACGtaacg/;
    return scalar(reverse($_));
}

=item iupac_revcomp(nt_sequence_string)

Reverse and complement nt sequence string, IUPAC code aware.
Will also complement IUPAC codes where relevant, e.g. Y to R.

=cut

sub iupac_revcomp {
    $_ = shift;

    tr/ATUGCRYMKBDVHatugcrymkbdvh/TAACGYRKMACTGtaacgyrkmactg/;
    return scalar(reverse($_));
}

=back

=cut
