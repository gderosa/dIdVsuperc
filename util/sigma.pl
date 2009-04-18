#!/usr/bin/perl -w

# transform "X Y1 Y2 Y3 ..."
# into: "X Y_mean sigma_Y"

# usage: ./sigma.pl < input.txt > output.dat

use strict;

main();

sub main {
    my @numbers = ();
    my ($x, $ymean, $sigmay);

    while(<STDIN>) { 
        @numbers = split;               # @numbers = (X, Y_1, Y_2, ..., YN);

        ($x, $ymean, $sigmay) = x_ymean_sigmay(@numbers); 
                                        # outputs: X Y_mean sigma_Y

        print "$x \t $ymean \t $sigmay \n";
    }
}    

sub x_ymean_sigmay {
    my $x       = shift;
    my @y       = @_;
    my $ymean   = mean(@y);
    my $sigmay  = sqrt( mean(squares(@y)) - $ymean**2 ); 
    
    return ($x, $ymean, $sigmay);
}

sub mean {
    my @y       = @_;
    my $n       = @y;                   # array length
    my $sum     = 0;
    my $i       = 0;

    for ($i=0;$i<$n;$i++) {
        $sum += $y[$i];
    }

    return $sum/$n;
}

sub squares {
    my @retval  = ();

    foreach (@_) {
        push @retval, $_**2;
    }

    return @retval;
}
