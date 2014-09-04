#!/usr/bin/perl
use warnings;
use strict;


print "enter coordinates (XXX XXX): " ;
chomp (my $coords = <STDIN>);
my @parts = split(' ', $coords);

my $c1 = $parts[1] - 8;
my $c2 = $parts[0] - 108;

print "$c1 $c2\n";

exit(0);
