#!/home/ivan/bin/perl -s

use blib;
use strict;
use warnings;
use Chemistry::File::SMARTS;
use Chemistry::File::SMILES;

our $debug ||= 0;
our $permute ||= 0;
our $overlap ||= 1;

$Chemistry::File::SMARTS::DEBUG = $debug;
$Chemistry::Pattern::Debug = $debug;

my $smarts = $ARGV[0] || 'C(OC)C';
my %options = (permute => $permute, overlap => $overlap);

#print "$smarts\n";

print "Pattern: $smarts\n";
print "Options: ", join(" ", %options), "\n";

my $patt = Chemistry::Pattern->parse($smarts, format => "smarts");
die "Invalid SMARTS" unless $patt;
$patt->options(%options);
#print "pattern compiled\n";

# Test matching on a smiles molecule
my $smiles = $ARGV[1] || "COCC";
print "Mol: $smiles\n";
my $mol = Chemistry::Mol->parse($smiles, format => 'smiles');

my @ret;
while ($patt->match($mol) ) {
    @ret = $patt->atom_map;
    print "Matched: (@ret)\n";
}
print "Matched: ()\n";



