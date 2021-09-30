#!/usr/bin/perl -w
use strict;

my @files = glob "*sn00113*";

foreach ( @files ) {
   # (1) Set the original file name to $f
   my($f) = $_;
   # (2) Set the new file name to $_
   s/sn/fn/;
   # (3) Espace special characters 
   #$din =~ s/,/\\,/g;
   #$din =~ s/&/\\&/g;
   #$din =~ s/ /\\ /g;
   #$_   =~ s/,/\\,/g;
   #$_   =~ s/&/\\&/g;
   #$_   =~ s/ /\\ /g;
   # (4) Rename
   if ($f ne $_) {
      print "input file  = $f\n";
      print "output file = $_\n";
      system "mv $f $_";
      #print "mv $f $_\n";
   }
}

