#!/usr/bin/perl
use warnings;
use strict;

my @matrix;

for (my $i = 0; $i < 720; $i++) {
	for(my $j = 0; $j<666; $j++) {
		$matrix[$i][$j] = 0;
	}
}

my $infile = 'flow_directions_corretcted_MS';
open (IN, "$infile");
chomp (my @lines = <IN>);
close(IN);

for (my $i=0; $i<=$#lines; $i++) {
	my @parts = split(' ', $lines[$i]);
	my $first = $parts[0] - $parts[2];
	my $second = $parts[1] - $parts[3];
	my $val;
	if ($first == 0) {
		if ($second == -1) {
			$val = 4;
		} elsif ($second == 0) {
			$val = 0;
		} else {
			$val = 8;
		}
	} elsif ($first == 1) {
		if ($second == -1 ) {
			$val = 3;
		} elsif ($second == 0) {
			$val = 2;
		} else {
			$val = 1;
		}
	} else {
		if ($second == -1) {
			$val = 5;
		} elsif ($second == 0) {
			$val = 6;
		} else {
			$val = 7;
		}
	}
	$matrix[$parts[0]][$parts[1]] = $val;
}

my $outfile = 'flows_visualized';
open (OUT, ">$outfile");
for (my $i = 0; $i < 720; $i++) {
	for(my $j = 0; $j<665; $j++) {
		print  OUT "$matrix[$i][$j] ";
	}
	print OUT "$matrix[$i][665]\n";
}
close(OUT);

exit(0);