#!/usr/bin/perl

$datadir="../../data";

opendir DATADIR,$datadir or die "Couldn't open $datadir";

open GNUPLOT,"|gnuplot";

while($f=readdir DATADIR) {
    if ($f =~ m/\.dat$/) {
        $f =~ s/(.*)\.dat$/$1/; 
        print GNUPLOT "name=\"$f\"\n";
        print GNUPLOT "load \"plot.gpi\"\n";
    }
}

close GNUPLOT;
close DATADIR;

