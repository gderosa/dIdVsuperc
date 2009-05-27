#!/usr/bin/perl -w

# Copyright 2007, "andye", http://www.perlmonks.org/?node_id=41758

# Adapted from http://www.perlmonks.org/?node_id=612779 by Guido De Rosa

# transforms "X Y1 Y2 Y3 ..."
# into: "X Y_mean sigma_Y JB-test-result-on-Y-distrib"

# usage: ./Jarque-Bera.pl < input.txt > output.dat

use PDL;
use strict;

my ($mean,$std_dev,$median,$min,$max,$adev,$rms);
my $source;
my $skewness;
my $exs_kurtosis;
my $JB;

my $x;
my @numbers;

while(<STDIN>) { 
    @numbers = split;

    $x = shift(@numbers);
    $source = pdl(@numbers);

    ($mean,$std_dev,$median,$min,$max,$adev,$rms) = stats($source);
        
    $skewness = sclr(
        sum(
            ($source - $mean)**3
        )  /  
        ( 
            ( nelem($source) - 1 ) * $std_dev**3 
        ) 
    );
            
    $exs_kurtosis = sclr(
        sum(
            ($source - $mean)**4
        ) / 
        ( ( nelem($source)-1 ) * $std_dev**4 ) -3 
     );
            
    $JB = ( nelem($source) / 6 ) * ($skewness**2 + ($exs_kurtosis**2 / 4) );
    
    print "$x\t$mean\t$std_dev\t$JB \n";
}
