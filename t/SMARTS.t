use Test::More;
use Chemistry::Smiles;
use Chemistry::Mol;
use Chemistry::Pattern;
use strict;

my @files;

BEGIN { 
    @files = glob "t/*.pat";
    plan tests => 1 + @files;

    use_ok('Chemistry::File::SMARTS');
};


my $mol_parser = new Chemistry::Smiles();

for my $file (@files) {
    open F, $file or die "couldn't open $file\n";   
    my ($patt_str, $options, $mol_str, @expected_matches) = map { /: ([^\n\r]*)/g } <F>;
    
    my ($mol, $patt);
    Chemistry::Atom->reset_id;
    $patt = Chemistry::Pattern->parse($patt_str, format => "smarts");
    $mol_parser->parse($mol_str, $mol = Chemistry::Mol->new);
    $patt->options(split " ", $options);

    my @matches;
    while ($patt->match($mol) ) {
        my @ret = $patt->atom_map;
        push @matches, "(@ret)";
    }
    push @matches, "()";

    is_deeply(\@matches, \@expected_matches, "$file: $mol_str =~ /$patt_str/");
}
