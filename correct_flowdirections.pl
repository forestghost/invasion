#!/usr/bin/perl
use warnings;
use strict;

my $infile = 'flow_directions';
my $outfile = 'flow_directions_corrected';

open (IN, "$infile");
open (OUT, ">$outfile");

chomp (my @lines = <IN>);
for (my $i = 0; $i<=$#lines;$i++) {
	my @parts = split(' ', $lines[$i]);
	my $first = $parts[0];
	my $second = $parts[1] - 99;
	my $third = $parts[2];
	my $fourth = $parts[3] - 99;
	print OUT "$first $second $third $fourth\n";
}
close IN;
close OUT;
exit(0);