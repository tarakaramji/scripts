#!/usr/bin/perl
use strict;
use warnings;

my $minsize=$ARGV[0];

while (my $line1 = <STDIN>) {
	my $line2 = <STDIN>;
	my $line3 = <STDIN>;
	my $line4 = <STDIN>;
	
	# Remove reads shorter than $minsize
	chomp $line2;
	my $readlength=length($line2);
	if ($readlength>=$minsize) {
		print $line1,$line2,"\n",$line3,$line4;		
	}
	
}
