#!/usr/bin/perl

use File::Basename;

$GNUPLOT = 'gnuplot';
# if you use Windows, uncomment the following, editing acccording to
# your needs (of course, you have to install Gnuplot binaries for Win32).
# You have to cd into this directory and issue 'perl plot-all.pl' from
# a command prompt
#
# $GNUPLOT = 'C:/gnuplot/bin/pgnuplot.exe';

$datadir = dirname($0)."/../../data";

opendir DATADIR,$datadir or die "Couldn't open $datadir";

open GNUPLOT,"|$GNUPLOT";

print GNUPLOT "cwd=\"".dirname($0)."\"\n";

while($f=readdir DATADIR) {
    if ($f =~ m/\.dat$/) {
        $name = $f;
        $name =~ s/(.*)\.dat$/$1/;
        print "Jarque-Bera graph for \"$name\"... ";
        print GNUPLOT "name=\"$name\"\n";
        print GNUPLOT "load \"".dirname($0)."/common.gpi\"\n";
        print GNUPLOT "load \"".dirname($0)."/Jarque-Bera.gpi\"\n";
        print "done.\n";
    }
    if ($f =~ m/^(.*)(\.dat)(.*)(\.fit)$/) {
        $name = $1;
        $fitfilename = $f;
        $mode = $3;
        print "Plotting $name$mode ... ";
        print GNUPLOT "name=\"$name\"\n";
        print GNUPLOT "fitfilename=\"$fitfilename\"\n";
        print GNUPLOT "mode=\"$mode\"\n";
        print GNUPLOT "load \"".dirname($0)."/common.gpi\"\n";
        print GNUPLOT "load \"".dirname($0)."/plot.gpi\"\n";
        print "done.\n";
    }
}

close GNUPLOT;
close DATADIR;

