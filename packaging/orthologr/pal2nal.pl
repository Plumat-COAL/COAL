#!/usr/bin/env perl
# A stand-in for PAL2NAL v14's pal2nal.pl that runs COAL instead.
#
# orthologr (drostlab/orthologr) runs the pal2nal.pl bundled inside the installed R package, once per
# ortholog pair:  perl <pkg>/pal2nal/pal2nal.v14/pal2nal.pl <aln> <cds> -output fasta > <out>
# It offers no option to point at another program, so this script can be copied over that file.
# It hands the unchanged PAL2NAL command line to COAL's `pal2nal` drop-in, which accepts PAL2NAL v14's
# options and writes PAL2NAL v14's output (see docs/differences.md for the documented differences).
#
# Install (R):
#   file.copy("pal2nal.pl", system.file("pal2nal/pal2nal.v14/pal2nal.pl", package = "orthologr"), overwrite = TRUE)
# Requires COAL's `pal2nal` command on PATH (pip install coal-align / conda install -c bioconda coal-align),
# or set COAL_PAL2NAL to its full path. Reinstalling orthologr restores the original file.
use strict;
use warnings;

# Which COAL command: COAL_PAL2NAL if set; else coal-cpp when installed (same output, starts in a few
# milliseconds, which matters when orthologr runs this once per gene pair); else COAL's `pal2nal` (Python).
# coal-cpp takes PAL2NAL's command line only for the options orthologr uses (the two files, -output, -codontable,
# -nogap, -nomismatch, -blockonly); anything else PAL2NAL accepts (-nostderr, -h / -help, PAL2NAL's own help) goes
# to the `pal2nal` drop-in, which understands all of them.
my @candidates;
if ($ENV{COAL_PAL2NAL}) {
    @candidates = ($ENV{COAL_PAL2NAL});
} elsif (coal_cpp_understands(@ARGV)) {
    @candidates = ("coal-cpp", "pal2nal");
} else {
    @candidates = ("pal2nal");
}
for my $exe (@candidates) {
    no warnings 'exec';
    exec($exe, @ARGV);    # returns only if it could not be started
}
print STDERR "pal2nal.pl (COAL stand-in): cannot run " . join(" or ", @candidates) . ": $!\n";
print STDERR "Install COAL (pip install coal-align, or conda install -c bioconda coal-align) so that coal-cpp or\n";
print STDERR "pal2nal is on PATH, or set COAL_PAL2NAL to the full path of one of them.\n";
exit 1;

# True when every argument is one coal-cpp takes on a PAL2NAL command line.
sub coal_cpp_understands {
    my @args = @_;
    my %flag = map { $_ => 1 } qw(-nogap -nomismatch -blockonly);
    my %with_value = map { $_ => 1 } qw(-output -codontable);
    my ($positional, $version) = (0, 0);
    while (@args) {
        my $a = shift @args;
        if ($with_value{$a}) {
            return 0 unless @args;
            my $v = shift @args;
            return 0 if $a eq "-codontable" && $v !~ /^\d+$/;
            return 0 if $a eq "-output" && $v !~ /^(fasta|clustal|paml|codon)$/;
        } elsif ($flag{$a}) {
            next;
        } elsif ($a eq "--version") {  # both engines print theirs
            $version = 1;
        } elsif ($a =~ /^-./) {    # any other option (-nostderr, -h, -help, ...): PAL2NAL-only or unknown
            return 0;
        } else {
            $positional++;
        }
    }
    return $positional >= 1 || $version;   # no files at all is a usage / help request: let pal2nal explain
}
