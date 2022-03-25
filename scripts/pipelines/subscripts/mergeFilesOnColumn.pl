#!/usr/bin/env perl 

#prints to stdout at the moment
#doesn't include column headings - need these to "merge" in order to this.
#longhand version of one-liner at:
#http://cgr.harvard.edu/cbg/scriptome/Windows/Tools/Merge.html
#"translation" done by Bela Tiwari, NEBC, March 2006

use strict;
use warnings;

if ($#ARGV < 3) { die "usage: < mergeFiles.pl filename1 filename2 columnNo1 columnNo2 >\n\n"; }

my ($file1,$file2,$col1_uncorr,$col2_uncorr) = @ARGV;

my $col1 = $col1_uncorr - 1;
my $col2 = $col2_uncorr - 1;

open (FILE1, $file1) or die "I can't open $file1\n";

my %line1 = ();

while (<FILE1>) 
{
     s/\r?\n//; 
     my @F=split /\t/, $_; 
     $line1{$F[$col1]} .= "$_\n";

}

warn "\nJoining $file1 column $col1_uncorr with $file2 column $col2_uncorr\n$file1: $. lines\n";

open (FILE2, $file2) or die "I can't open $file2\n";

my $merged = 0;

while (<FILE2>) 
{
     s/\r?\n//; 
     my @F=split /\t/, $_; 
     my $x = $line1{$F[$col2]};
     if ($x) { $x =~ s/\n/\t$_\n/g; 
             print $x; 
             $merged++; }
}

warn "$file2: $. lines\nMerged file: $merged lines\n";
