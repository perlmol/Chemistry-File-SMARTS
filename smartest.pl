#!/home/ivan/bin/perl

use blib;
use strict;
use warnings;
use Chemistry::File::SMARTS;

#Parse::RecDescent->Precompile($grammar, "SmilesGrammar");
#use SmilesGrammar;
#my $parser = SmilesGrammar->new() or die;

my $s = $ARGV[0] || 'C(OC)C';

print "$s\n";

my $patt = Chemistry::Pattern->parse($s, format => "smarts");
die "Invalid SMARTS" unless $patt;
print "pattern compiled\n";

# Test matching on a smiles molecule
use Chemistry::Smiles;
my $mol_parser = new Chemistry::Smiles();
my $mol;
$mol_parser->parse($ARGV[1] || "COCC", $mol = Chemistry::Mol->new);
my @ret;
while ($patt->match($mol) ) {
    @ret = $patt->atom_map;
    print "Matched: (@ret)\n";
}
print "Matched: ()\n";



