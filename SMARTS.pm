package Chemistry::File::SMARTS;

$VERSION = "0.01";

use 5.006001;
use strict;
use warnings;
use Chemistry::Pattern;
use base "Chemistry::File";
use Parse::RecDescent;
use Carp;
use Data::Dumper;

=head1 NAME

Chemistry::SMARTS - SMARTS parser

=head1 SYNOPSYS


    #!/usr/bin/perl
    use Chemistry::File::SMARTS;

    my $smarts = "CC(O)C";
    my $patt = Chemistry::Pattern->parse($smarts, format => "smarts")
        or die "invalid SMARTS!\n";
    
    # somehow get a molecule into $mol ...
    $mol = Chemistry::Mol->read("myfile.mol");

    while ($patt->match($mol) ) {
        my @atoms = $patt->atom_map;
        print "Matched: (@atoms)\n";
    }

=head1 DESCRIPTION

This module is used for reading SMARTS patterns, returning Chemistry::Pattern
objects.

=cut


# Initialization
Chemistry::Mol->register_format(smarts => __PACKAGE__);
my $parser = init_parser();
my $Debug = 0;

# Chemistry::File interface

sub parse_string {
    my ($self, $s, %options) = @_;
    my $result;
    my $patt;
    if (defined($result = $parser->smiles($s))) {
        # print "It's a valid SMARTS!\n";
        print Dumper($result) if $Debug;
        $patt = compile($result);
    } 
    $patt;
}

# Private stuff

sub init_parser {
    my $grammar = <<'END_GRAMMAR';
        # MAIN STRUCTURE
        smiles: atom_or_branch chain_or_branch(s?) /\z/ {[$item[1], @{$item[2]}]}
        atom_or_branch: atomd | branch
        chain_or_branch: nextatom | branch
        nextatom: bond_expr atomd {[nextatom => $item[1], $item[2]]}

        branch: "(" chain_or_branch(s) ")" {[branch => $item[2]]}

        # ATOMS
        atomd: atom bdigit(s?) {[@item]}
        atom: simple_atom 
            { [atom_expr => [atom_or_expr => [ atom_and_expr => $item[1]]]] }
            | complex_atom
        simple_atom: organic {["", symbol => $item[1]]} 
            | aromatic {["", symbol => $item[1]]}

        # CHEMICAL ELEMENTS
        # Order matters! Otherwise P would be found before Pt...
        element: inorganic | organic   
        # Cl must be before C and Br before B
        organic: "Cl" | "Br" | "B" | "C" | "N" | "O" | "S" | "P" | "F" | "I" 
        inorganic: "Ru" | "Pt" | "H" | "Na"
        aromatic: "c" | "n" | "o" | "s"

        # ATOMIC EXPRESSIONS
        complex_atom: "[" atom_expr "]" {$item[2]}
        atom_expr: atom_or_expr(s /;/) {[$item[0],@{$item[1]}]} 
        atom_or_expr: atom_and_expr(s /,/) {[$item[0],@{$item[1]}]}
        atom_and_expr: atom_not_expr(s /&?/) {[$item[0],@{$item[1]}]}
        atom_not_expr: /!?/ atomic_primitive {[$item[1],@{$item[2]}]}
        #atom_and_expr: atomic_primitive(s /&?/) {[$item[0],@{$item[1]}]}

        # BOND EXPRESSIONS
        bond_expr: bond_or_expr(s /;/) {[$item[0],@{$item[1]}]} 
        bond_or_expr: bond_and_expr(s /,/) {[$item[0],@{$item[1]}]}
        bond_and_expr: bond_not_expr(s /&?/) {[$item[0],@{$item[1]}]} 
            | "" {[$item[0], ["", ""] ]}
        bond_not_expr: /!?/ bond_primitive {[$item[1], $item[2]]} 
        bond_primitive: "-" | "=" | "#" | ":" | "." | "/" | "\\\\"

        # ATOMIC PRIMITIVES
        atomic_primitive: hydrogens | isotope | charge 
            | element { [symbol => $item[1]]} 
        hydrogens: "H" /\d?/ {[ hcount => $item[2] || 1]}
        isotope: /\d+/ {[iso => $item[1]]}
        charge: /((?:\+|-)+)(\d*)/ {
            my $chg;
            if(length($1) > 1) {
                $chg = ($2 eq '') ? substr($1, 0, 1).(length($1)) : undef;
            } else {
                $chg = "${1}1" * ($2 ne '' ? $2 : 1);
            }
            [charge => $chg];
        }

        # DIGITS
        bdigit: bond_expr digit   {[@item]}
        digit: /(\d)/ | /%(\d\d)/ {$1}

END_GRAMMAR


    $Parse::RecDescent::skip = "";
    my $parser = Parse::RecDescent->new($grammar) or die;
    $parser;
}


# compiles a parse tree and returns a Chemistry::Pattern object.
# The tree itself is destroyed!
sub compile {
    my ($tree) = @_;
    my @stack;
    my @last_atom;
    my $patt = Chemistry::Pattern->new;
    my %digits;

    while (1) {
        my $tok = shift @$tree;
        unless ($tok) {
            last unless @stack;
            $tree = pop @stack;
            pop @last_atom;
            next;
        }
        my $type = $tok->[0];
        print "tok: $type\n" if $Debug;

        if ($type eq 'branch') {
            push @stack, $tree;
            $tree = $tok->[1];
            push @last_atom, $last_atom[-1];
        } elsif ($type eq 'atomd') {
            my $new_atom = $patt->new_atom(test_sub => atom_sub($tok->[1]));
            push @last_atom, $new_atom;
            for my $digit (@{$tok->[2]}) {
                do_digit($patt, \%digits, $digit, $new_atom);
            }
        } elsif ($type eq "nextatom") {
            print "a bond!\n" if $Debug;
            #print Dumper($tok);
            my $old_atom = pop @last_atom;
            my $new_atom = $patt->new_atom(test_sub => atom_sub($tok->[2][1]));
            push @last_atom, $new_atom;
            $patt->new_bond(atoms => [$old_atom, $new_atom], 
                test_sub => bond_sub($tok->[1]));
            for my $digit (@{$tok->[2][2]}) {
                do_digit($patt, \%digits, $digit, $new_atom);
            }
        } elsif (ref $tok eq "ARRAY") {
            die;
            push @stack, $tree;
            $tree = $tok;
        } else {
            die;
        }
    }
    #print Dumper ($patt);
    $patt;
}

sub do_digit {
    my ($patt, $digits, $dig, $atom) = @_;
    my ($bnd, $d) = ($dig->[1], $dig->[2]);

    if ($digits->{$d}) {
        print "closing ring $d\n" if $Debug;
        $patt->new_bond(atoms=>[$atom, $digits->{$d}{atom}], 
            test_sub => bond_sub($bnd));
        delete $digits->{$d};
    } else {
        print "opening ring $d\n" if $Debug;
        $digits->{$d} = {bond => $bnd, atom => $atom};
    }
}

sub atom_sub {
    my $atom = shift;
    my $expr = atom_expr($atom);

    eval <<SUB
    sub {
        my (\$self, \$atom) = \@_;
        $expr;
    }
SUB
}

sub atom_expr {
    my $atom = shift;
    die unless $atom->[0] eq 'atom_expr';
    my $expr;
    if (@$atom >= 2) {
        my @or_exprs = map {atom_or_expr($_)} @{$atom}[1 .. $#$atom];
        $expr = join(") && (", @or_exprs);
        $expr = "($expr)";
    } else {
        die;
    }
    print "expr: $expr\n" if $Debug;
    $expr;
}

sub atom_or_expr {
    my $term = shift;
    die unless $term->[0] eq 'atom_or_expr';
    my $expr;
    if (@$term >= 2) {
        my @exprs = map {atom_and_expr($_)} @{$term}[1 .. $#$term];
        $expr = join " || ", @exprs;
    } else {
        die;
    }
    $expr;
}


sub atom_and_expr {
    my $term = shift;
    die unless $term->[0] eq 'atom_and_expr';
    my $expr;
    if (@$term >= 2) {
        my @exprs = map {atom_primitive_expr($_)} @{$term}[1 .. $#$term];
        $expr = join " && ", @exprs;
    } else {
        die;
    }
    $expr;
}

sub atom_primitive_expr {
    my $prim = shift;

    my $expr = "\$atom->$prim->[1] eq '$prim->[2]'";
    $expr = "!($expr)" if $prim->[0];
    $expr;
}

sub bond_sub {
    my $bond = shift;
    my $expr = bond_expr($bond);

    eval <<SUB
    sub {
        my (\$self, \$bond) = \@_;
        $expr;
    }
SUB
}

sub bond_expr {
    my $bond = shift;
    die unless $bond->[0] eq 'bond_expr';
    my $expr;
    if (@$bond >= 2) {
        my @or_exprs = map {bond_or_expr($_)} @{$bond}[1 .. $#$bond];
        $expr = join(") && (", @or_exprs);
        $expr = "($expr)";
    } else {
        die;
    }
    print "expr: $expr\n" if $Debug;
    $expr;
}

sub bond_or_expr {
    my $term = shift;
    die unless $term->[0] eq 'bond_or_expr';
    my $expr;
    if (@$term >= 2) {
        my @exprs = map {bond_and_expr($_)} @{$term}[1 .. $#$term];
        $expr = join " || ", @exprs;
    } else {
        die;
    }
    $expr;
}


sub bond_and_expr {
    my $term = shift;
    die unless $term->[0] eq 'bond_and_expr';
    my $expr;
    if (@$term >= 2) {
        my @exprs = map {bond_primitive_expr($_)} @{$term}[1 .. $#$term];
        $expr = join " && ", @exprs;
    } else {
        die;
    }
    $expr;
}

sub bond_primitive_expr {
    my $prim = shift;

    $prim =~ s/\\/\\\\/;
    my $expr = "\$bond->type eq '$prim->[1]'";
    $expr = "!($expr)" if $prim->[0];
    $expr;
}

1;

=head1 BUGS

=head1 SEE ALSO

For more information about SMARTS, see the SMARTS Theory Manual at
http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2003 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

