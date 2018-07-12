#!/usr/bin/perl
# Copyright 2009 M. Goerz
#
# This file is part of qdyn.
#
# qdyn is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# qdyn is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with qdyn.  If not, see <http://www.gnu.org/licenses/>.
use strict;
use warnings;
use CGI::Pretty; # core module, see http://perldoc.perl.org/CGI.html
use Storable; # core module, see http://perldoc.perl.org/Storable.html

if (@ARGV != 2) {
    my $usage = "Usage: doc_indexgen.pl INDEXFILE OUTFILE\n";
    die($usage);
}
my $index_file = $ARGV[-2];
my $outfile = $ARGV[-1];
print "Creating index from $index_file, storing in $outfile\n";

my $global_index = {}; # public symbol (including from other files) => link

if (-f $index_file){
    $global_index = retrieve($index_file);
} else {
    die ("Can't find $index_file\n");
}

my $cgi = new CGI;

#  return $p1 with those leading parts of the path removed that are also
#  present in $p2. E.g. 
#  fix_path('doc/reference/oct.html', 'doc/lib_index.html')
#  returns 'reference/oct.html'
sub fix_path{
    my $p1 = shift;
    my $p2 = shift;
    my @parts1 = split(/\//, $p1);
    my @parts2 = split(/\//, $p2);
    my $dir1 = $parts1[0];
    my $dir2 = $parts2[0];
    while ($dir1 eq $dir2){
        $dir1 = shift(@parts1);
        $dir2 = shift(@parts2);
    }
    my $result = $dir1.'/'.join("/", @parts1);
    return $result;
}

open(OUTFILE, ">$outfile") or die ("Couldn't open $outfile\n");
print OUTFILE $cgi->start_html(
    -title=>"Library Index",
    -style=>{'src'=>'../style.css'}
);
print OUTFILE $cgi->h1("Index");
print OUTFILE "<ul>\n";
foreach my $symbol (sort(keys(%{$global_index}))){
    my $entry = $global_index->{$symbol};
    my ($kind, $url) = split(/;/, $entry);
    $url = fix_path($url, $outfile);
    $url =~ s'/$'';
    print OUTFILE " <li>$kind <a href=\"$url\">$symbol</a></li>\n";
}
print OUTFILE "</ul>\n";

print OUTFILE $cgi->end_html();
close OUTFILE;
