#!/home/ivan/bin/perl -s

use blib;
use strict;
use warnings;
use Chemistry::File::SMARTS;
use Chemistry::File::SMILES;
use Chemistry::Ring 'aromatize_mol';

our ($debug, $debug_pattern);
our $permute ||= 0;
our $overlap ||= 1;

$Chemistry::File::SMARTS::DEBUG = $debug;
$Chemistry::Pattern::DEBUG = $debug_pattern;

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
my $mol = Chemistry::Mol->parse($smiles, format => 'smiles');
print "Mol: $smiles\n";
my $asmiles = $mol->print(format => 'smiles', aromatic => 1);
my $ksmiles = $mol->print(format => 'smiles');
#print "AMol: $asmiles\n";
#print "KMol: $ksmiles\n";
aromatize_mol($mol);

my @ret;
while ($patt->match($mol) ) {
    @ret = $patt->atom_map;
    print "Matched: (@ret)\n";
}
print "Matched: ()\n";



