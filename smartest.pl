#!/home/ivan/bin/perl -s

use blib;
use strict;
use warnings;
use Chemistry::File::SMARTS;

our $Debug ||= 0;
our $permute ||= 0;
our $overlap ||= 1;

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
use Chemistry::Smiles;
my $mol_parser = new Chemistry::Smiles();
my $mol;
my $smiles = $ARGV[1] || "COCC";
print "Mol: $smiles\n";
$mol_parser->parse($smiles, $mol = Chemistry::Mol->new);
my @ret;
while ($patt->match($mol) ) {
    @ret = $patt->atom_map;
    print "Matched: (@ret)\n";
}
print "Matched: ()\n";



