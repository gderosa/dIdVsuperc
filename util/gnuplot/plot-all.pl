#!/usr/bin/perl

use File::Basename;

$datadir = dirname($0)."/../../data";

opendir DATADIR,$datadir or die "Couldn't open $datadir";

open GNUPLOT,"|gnuplot";

print GNUPLOT "cwd=\"".dirname($0)."\"\n";

while($f=readdir DATADIR) {
    if ($f =~ m/\.dat$/) {
        $f =~ s/(.*)\.dat$/$1/; 
        print GNUPLOT "name=\"$f\"\n";
        print GNUPLOT "load \"".dirname($0)."/plot.gpi\"\n";
    }
}

close GNUPLOT;
close DATADIR;

