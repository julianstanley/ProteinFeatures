#!/usr/bin/perl -w

###############################################################################
# sspider_RLC_mech.pl: This script makes automatic requests to SPPIDER 2 with #
# default settings. It requires a single input argument, the path to a PDB    #
# file (e.g. C:/cygwin64/home/Roger/devel/allStruct/1bf5_A.pdb) upon which    #
# protein interaction residues are to be predicted.                           #
# Author: Roger L. Chang                                                      #
# This version created December 7, 2016                                       #
# Updated March 10, 2020                                                      #
###############################################################################

use strict;
use warnings;
use WWW::Mechanize;

# Process arguments.
@ARGV == 1 or die "usage: sspider_RLC_mech.pl <pdb_file> > <outputfile> \n";
# Open the file and get the data.
my($pdb_file) = $ARGV[0];
chomp $pdb_file;

# The server file field requires using "\" instead of "/" to denote subfolders.
# To correctly format the file name for recogition, you must use replace "/" with "\\".
$pdb_file =~ s/\//\\/g;

print "$pdb_file\n";




my $email = 'julian_stanley@hms.harvard.edu';
#my $email = 'roger.l.chang@gmail.com';

my $mech = WWW::Mechanize->new();
my $url = 'http://sppider.cchmc.org/';
$mech->get($url);

$mech->field('TQuery',1);
$mech->field('Version',2);
$mech->field('EMail',$email);
#$mech->field('PDBCode','1qav'); # For testing, "1qav" works
#$mech->field('PDBChain','A'); # For testing, "1qav" works
$mech->field('PDBFileName' => $pdb_file);
$mech->field('Trade',0.5);
$mech->field('PDBres',1);
$mech->field('BFgradient',0);
$mech->field('CASPfmt',1);
$mech->submit();

sleep 15; # A short delay may be needed between file submissions to avoid errors.

exit;
