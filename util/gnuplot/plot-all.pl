#!/usr/bin/perl

use File::Basename;

$datadir = dirname($0)."/../../data";

opendir DATADIR,$datadir or die "Couldn't open $datadir";

open GNUPLOT,"|gnuplot";

print GNUPLOT "cwd=\"".dirname($0)."\"\n";

while($f=readdir DATADIR) {
    if ($f =~ m/\.dat$/ and -f "$datadir/$f.fit") {
        $f =~ s/(.*)\.dat$/$1/; 
        print "Plotting \"$f\"... ";
        print GNUPLOT "name=\"$f\"\n";
        print GNUPLOT "load \"".dirname($0)."/plot.gpi\"\n";
        print "done.\n";
    }
}

close GNUPLOT;
close DATADIR;

