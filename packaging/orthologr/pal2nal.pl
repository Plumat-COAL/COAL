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
my @candidates = $ENV{COAL_PAL2NAL} ? ($ENV{COAL_PAL2NAL}) : ("coal-cpp", "pal2nal");
for my $exe (@candidates) {
    no warnings 'exec';
    exec($exe, @ARGV);    # returns only if it could not be started
}
print STDERR "pal2nal.pl (COAL stand-in): cannot run " . join(" or ", @candidates) . ": $!\n";
print STDERR "Install COAL (pip install coal-align, or conda install -c bioconda coal-align) so that coal-cpp or\n";
print STDERR "pal2nal is on PATH, or set COAL_PAL2NAL to the full path of one of them.\n";
exit 1;
