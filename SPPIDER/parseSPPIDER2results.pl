#!/usr/bin/perl -w

###############################################################################
# parseSPPIDER2results.pl: Takes as input raw results from SPPIDER2 received  #
# as emails and concatenated into one text file and parses out the structure  #
# filename and the predicted PPI binding site residues.                       #
# Author: Roger L. Chang                                                      #
# This version created March 21, 2017                                         #
###############################################################################

use strict;
use warnings;

@ARGV == 1 or die "usage: parseSPPIDER2results.pl <SPPIDER2_file> > <outputfile> \n";

# Process arguments.
# Open the file and get the data.
my($SPPIDER2_file) = $ARGV[0];
open(SPPIDER2_FILE, "< $SPPIDER2_file")
	or die "Can't open $SPPIDER2_file for reading: $!\n";
my @data = <SPPIDER2_FILE>;

# Initialize main variables.
my $target = "";
my $siteFlag = 0;
my $bindingSite = "";

# Loop through and parse necessary data.
foreach (@data) {
	chomp $_;
	# Get target protein structure filename.
	if ($_ =~ /TARGET ./) {
		($target) = ($_ =~ /TARGET (\S*pdb)/);
		chomp $target;
		$target =~ s/\n//;
#		print "$target\t";
	}

	# Find line just prior to binding site residue data.
	if ($_ =~ /^Comment: numbering in the original PDB formatted file:/) {
		$siteFlag = 1;
		next;
	}

	# Compile binding site residue data.
	if (($siteFlag == 1) && ($_ =~ /^Comment: \D/)) {
		my $sitePart = $_;
		$sitePart =~ s/Comment: //;
		$sitePart =~ s/\s//;
		if (length($bindingSite) > 0) {
			$bindingSite .= ",";
		}
		$bindingSite .= $sitePart;
	}

	# Find end of binding site data, print site, and reinitialize variables.
	if (($siteFlag == 1) && ($_ =~ /^END/)) {
		$bindingSite = substr($bindingSite, 0, -1);
#		print "$bindingSite\n";

		# Split residues and print out site ID ("structureFilename_1-letterCodePosition")
		my(@residues) = split(',', $bindingSite);
		foreach my $residue (@residues) {
			chomp $residue;
			print "$target"."_"."$residue\tPPI\n";
		}

		$target = "";
		$siteFlag = 0;
		$bindingSite = "";
	}
}
close(SPPIDER2_FILE);

exit;
