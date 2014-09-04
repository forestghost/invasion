#!/usr/bin/perl
use warnings;
use strict;

my $infile = 'corrected_mississippi_flows';
open (IN, "$infile");
chomp(my @lines = <IN>);
close IN;

my %new;

for (my $i=0; $i<=$#lines;$i++) {
	my @parts = split(' ', $lines[$i]);
	my $key = $parts[0].' '.$parts[1];
	my $value = $parts[2].' '.$parts[3];
	$new{$key} = $value;
}

$infile = 'flow_directions_LASTLY';

open (IN, "$infile");
chomp(my @lines2 = <IN>);
close IN;

my $outfile = 'flow_directions_corretcted_MS';
open (OUT, ">$outfile");

for (my $i=0; $i<=$#lines2;$i++) {
	my @parts = split(' ', $lines2[$i]);
	my $key = $parts[0].' '.$parts[1];
	if (exists $new{$key}) {
		print OUT "$key $new{$key}\n";
	} else {
		print OUT "$lines2[$i]\n";
	}
}


exit(0);