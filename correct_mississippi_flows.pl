#!/usr/bin/perl
use warnings;
use strict;


my $infile = 'mississippi_raster';
open (IN, "$infile");
chomp (my @lines = <IN>);
close IN;

my @ms_matrix;

for (my $i = 0; $i < 720; $i++) {
	my @parts = split(' ', $lines[$i]);
	for(my $j = 0; $j<666; $j++) {
		$ms_matrix[$i][$j] = $parts[$j];
	}
}

$infile = 'flow_directions_LASTLY';
open (IN, "$infile");
chomp (@lines = <IN>);
close IN;

my @prev_flows;

for (my $i=0; $i<=$#lines; $i++) {
	my @parts = split(' ' , $lines[$i]);
	my $xdiff = $parts[0] - $parts[2];
	my $ydiff = $parts[3] - $parts[1];
	my $change;
	if ($xdiff == 0) {
		if ($ydiff == -1) {
			$change = 2;
		} else {
			$change = 6;
		}
	} elsif ($xdiff == -1) {
		if ($ydiff == -1) {
			$change = 1;
		} elsif ($ydiff == 0) {
			$change = 8;
		} else {
			$change = 7;
		}
	} elsif ($xdiff == 1) {
		if ($ydiff == -1) {
			$change = 3;
		} elsif ($ydiff == 0) {
			$change = 4;
		} else {
			$change = 5;
		}
	}
	$prev_flows[$parts[0]][$parts[1]]  = $change;
}

my $outfile = 'corrected_mississippi_flows';
open (OUT, ">$outfile");

for (my $i =0; $i<719; $i++) {
	for (my $j=0; $j<665; $j++) {
		if ($ms_matrix[$i][$j] == 0) { #then MS cell
			my ($a, $b, $l, $r, $h, $v);
			$a=0 ; $b=0; $l=0; $r=0; $h=0; $v =0;
			my $cur = -1;
			while (1) {  #above ($a)
				if ($ms_matrix[$i][$j+$cur] ==0) {
					$a++;
					$cur--;
				} else {
					last;
				}
			}
			$cur = 1;
			while (1) { # below ($b)
				if ($ms_matrix[$i][$j+$cur] == 0) {
					$b++;
					$cur++;
				} else {
					last;
				}
			}
			$cur = -1;
			while(1) { #left ($l) 
				if ($ms_matrix[$i+$cur][$j] == 0) {
					$l++;
					$cur--;
				} else {
					last;
				}
			}
			$cur = 1;
			while(1) { #right ($r)
				if ($ms_matrix[$i+$cur][$j] == 0 ) {
					$r++;
					$cur++;
				} else {
					last;
				}
			}
			$h = 1 + $l + $r;
			$v = 1 + $a + $b;
print "$i $j $a $b $l $r $h $v\n";			
			if ($h < $v) {  # river in vertical orientation at the point analyzed
				if ($ms_matrix[$i][$j+1] == 0 &&  ( $prev_flows[$i][$j] != 2 || $prev_flows[$i][$j] != 1 || $prev_flows[$i][$j] != 3  )  ) {
					my $q = $j+1;	print OUT "$i $j $i $q\n";
print "h<v: $i $j $i $q\n"					
				}	
			}
			if ($h > $v) { # river in lateral orientation at the point analyzed
				if ( $ms_matrix[$i+1][$j] == 0  && (  $prev_flows[$i][$j] != 1 || $prev_flows[$i][$j] != 8 || $prev_flows[$i][$j] != 7) ) {
					my $q = $i+1; print OUT "$i $j $q $j\n";
print "h>v: $i $j $q $j\n"
				}
			}
			#  o.w., new movement conflicts with previous movement or $h = $v ; in these cases, no changes are made to existing flow directions
		}
	}
}


close OUT;

exit(0);