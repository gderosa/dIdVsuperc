#!/usr/bin/perl

use File::Basename;

$GNUPLOT = 'gnuplot';

# if you use Windows, uncomment the following, editing acccording to
# your needs (of course, you have to install Gnuplot binaries for Win32).
# $GNUPLOT = 'C:/gnuplot/bin/pgnuplot.exe';

chdir dirname($0);

$datadir = "../../data";

opendir DATADIR,$datadir or die "Couldn't open data dir!";

open GNUPLOT,"|$GNUPLOT";

while($f=readdir DATADIR) {
    if ($f =~ m/^(.*)(\.dat)(.*)(\.fit)$/) {
        $name = $1;
        $fitfilename = $f;
        $mode = $3;
	$mode =~ s/^\.//; 
        print "Plotting best fit for \"$name\" ($mode) ... \n";
        print GNUPLOT "name=\"$name\"\n";
        print GNUPLOT "fitfilename=\"$fitfilename\"\n";
        print GNUPLOT "mode=\"$mode\"\n";
        print GNUPLOT "load \"common.gpi\"\n";
        print GNUPLOT "load \"plot.gpi\"\n";
	#print "done.\n";
    }
}

close GNUPLOT;
close DATADIR;

# Uncomment the following lines if you "point and click" to this script
# and just don't want the command window to disappear...
#
# print "Done. Hit ENTER to exit.\n";
# read STDIN, $dummy, 1;
