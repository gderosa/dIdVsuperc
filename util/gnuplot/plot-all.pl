#!/usr/bin/perl

use File::Basename;

$datadir = dirname($0)."/../../data";

opendir DATADIR,$datadir or die "Couldn't open $datadir";

open GNUPLOT,"|gnuplot";

print GNUPLOT "cwd=\"".dirname($0)."\"\n";

while($f=readdir DATADIR) {
    if ($f =~ m/\.dat$/) {
        $n = $f;
        $n =~ s/(.*)\.dat$/$1/;
        print "Jarque-Bera graph for \"$n\"... ";
        print GNUPLOT "name=\"$n\"\n";
        print GNUPLOT "load \"".dirname($0)."/common.gpi\"\n";
        print GNUPLOT "load \"".dirname($0)."/Jarque-Bera.gpi\"\n";
        print "done.\n";
        
        if (-f "$datadir/$f.fit") {
            print "Plotting \"$n\"... ";
            print GNUPLOT "load \"".dirname($0)."/plot.gpi\"\n";
            print "done.\n";
        }
    }
}

close GNUPLOT;
close DATADIR;

