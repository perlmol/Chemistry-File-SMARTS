package Chemistry::File::SMARTS;

$VERSION = "0.01";

use 5.006;
use strict;
use warnings;
use Chemistry::Pattern;
use base "Chemistry::File";
use Carp;
#use Data::Dumper;
use List::Util 'sum';
#use Text::Balanced qw(extract_multiple extract_bracketed);

# Initialization
Chemistry::Mol->register_format(smarts => __PACKAGE__);
our $DEBUG = 0;

# Chemistry::File interface

sub parse_string {
    my ($self, $s, %options) = @_;
    my $patt = parse_smarts($s);
    $patt;
}

sub parse_smarts {
    my ($s) = @_;
    my %digits;

    my @toks = tokenize($s);
    my $tok = shift(@toks);
    my $patt = Chemistry::Pattern->new();
    my @atom_stack;
    my $current_atom = parse_atom($patt, $tok);
    while (defined ($tok = shift @toks)) {
        print "tok: $tok\n" if $DEBUG;
        if ($tok eq '(') {
            push @atom_stack, $current_atom;
        } elsif ($tok eq ')') {
            $current_atom = pop @atom_stack;
        } else {  # bond, atom
            my $next_tok = shift @toks;
            if ($next_tok =~ /^\d+$/) {  # digit
                if ($digits{$next_tok}) {  # close ring
                    parse_bond($patt, $tok, $current_atom, 
                        $digits{$next_tok});
                    $digits{$next_tok} = undef;
                } else { # open ring
                    $digits{$next_tok} = $current_atom;
                }
            } else {
                my $next_atom = parse_atom($patt, $next_tok);
                parse_bond($patt, $tok, $current_atom, $next_atom);
                $current_atom = $next_atom;
            }
        }
    }
    $patt;
}

sub parse_atom {
    my ($patt, $s) = @_;
    
    my $expr = 
        join " and ", map { 
            join " or ", map { 
                join ' && ', map {
                    parse_atomic_primitive($_);
                } split '&', $_;
            } split ',', $_;
        } split ';', $s;

    print "atom expr: $expr\n" if $DEBUG;
    my $sub = eval <<SUB;
        sub {
            my (\$patt, \$atom) = \@_;
            $expr;
        };
SUB
    my $atom = Chemistry::Pattern::Atom->new(test_sub => $sub);
    $patt->add_atom($atom);
    $atom;
}


# missing primitives: R, r, @, @@
sub parse_atomic_primitive {
    local ($_) = @_;
    my @terms;
    no warnings 'uninitialized';
    s/(!?)\*// &&           # wildcard
        push @terms, "${1}1";
    s/(!?)D(\d?)// &&       # explicit connections (shouldn't count h)
        push @terms, "$1(\$atom->bonds == " . (defined $2 ? $2 : 1) . ')';
    s/(!?)a// &&            # aromatic
        push @terms, "$1\$atom->aromatic";
    s/(!?)A// &&            # aliphatic
        push @terms, "$1(!\$atom->aromatic)";
    s/(!?)X(\d?)// &&       # total connections (should add implicit H)
        push @terms, "$1(\$atom->bonds == " . (defined $2 ? $2 : 1) . ')';
    s/(!?)v(\d?)// &&       # valence
        push @terms, "$1(sum(map {\$_->order} \$atom->bonds) == " 
            . (defined $2 ? $2 : 1) . ')';
    s/(!?)[Hh](\d?)// &&    # H-count
        push @terms, "$1(sum(map {\$_->symbol eq 'H'} \$atom->neighbors) == " 
            . (defined $2 ? $2 : 1) . ')';
    s/(!?)#(\d+)// &&       # atomic number
        push @terms, "$1(\$atom->Z == $2)";
    s/(!?)([+-]\d+)// &&    # numerical charge 
        push @terms, "$1(\$atom->formal_charge == $2)";
    s/(!?)(\++)// &&        # positive charge
        push @terms, "$1(\$atom->formal_charge == " . length $2 . ')';
    s/(!?)(-+)// &&         # negative charge 
        push @terms, "$1(\$atom->formal_charge == -" . length $2 . ')';
    s/(!?)(\d+)// &&        # mass
        push @terms, "$1(\$atom->mass == $2)";
    s/(!?)([cnosp])// &&    # aromatic symbol
        push @terms, "$1(\$atom->symbol eq '$2' && \$atom->aromatic)";
    s/(!?)([A-Z][a-z]?)// &&    # aliphatic symbol
        push @terms, "$1(\$atom->symbol eq '$2' && ! \$atom->aromatic)";
    join ' && ', @terms;
}

sub parse_bond {
    my ($patt, $s, @atoms) = @_;

    my $expr;
    
    if ($s) {
        $expr = 
            join " and ", map { 
                join " or ", map { 
                    join ' && ', map {
                        parse_bond_primitive($_);
                    } split '&', $_;
                } split ',', $_;
            } split ';', $s;
    } else {
        $expr = '($bond->order == 1 || $bond->aromatic)';
    }

    print "bond expr: $expr\n" if $DEBUG;
    my $sub = eval <<SUB;
        sub {
            my (\$patt, \$bond) = \@_;
            $expr;
        };
SUB
    my $bond = Chemistry::Pattern::Bond->new(test_sub => $sub, 
        atoms => \@atoms);
    $patt->add_bond($bond);
    $bond;
}

sub parse_bond_primitive {
    local ($_) = @_;
    my @terms;
    s/(!?)~// &&        # wildcard
        push @terms, "${1}1";
    s/(!?)-// &&        # single
        push @terms, "$1(\$bond->order == 1)";
    s/(!?)=// &&        # double
        push @terms, "$1(\$bond->order == 2)";
    s/(!?)#// &&        # triple
        push @terms, "$1(\$bond->order == 3)";
    s/(!?):// &&        # triple
        push @terms, "$1\$bond->aromatic";
    join ' && ', @terms;
}

my %ORGANIC_ELEMS = (
    Br => 1, Cl => 1, B => 1, C => 1, N => 1, O => 1, P => 1, S => 1, 
    F => 1, I => 1, s => 1, p => 1, o => 1, n => 1, c => 1, b => 1,
);

sub tokenize {
    my ($s) = @_;
    my $state = 3;
    my $paren_depth = 0;
    my $atom;
    my $digit;
    my $rec_smart;
    my @rec_smarts;
    my $bond = '';
    my $symbol;
    my @toks;
    
    # 0 expects atom or branch
    # after atom: atom or bond or branch
    # after bond: atom or branch or /,;/
    # token types: atom, bond, (, )
    my @chars = split '', $s;
    my $char;
    while (defined ($char = shift @chars)) {
        print "char: $char\n" if $DEBUG;
        if ($state == 0) { # expect atom or branch (not used!)
            push(@toks,  $char), $state = 2, next if ($char =~ /[()]/);
            $state = 3, redo;
        } elsif ($state == 1) { # in complex atom
            if ($char eq ']') {
                push @toks, $atom, @rec_smarts;
                @rec_smarts = ();
                $state = 4, next;
            }
            $atom .= $char;
            $state = 5 if ($char eq '$'); # recursive smarts
            next;
        } elsif ($state == 2) { # expect bond
            if ($char =~ /[-=:~\\\/#@?,;&]/) {
                $bond .= $char;
                next;
            } else {
                push @toks, $bond;
                $bond = '';
                $state = 3, redo;
            }
        } elsif ($state == 3) {  # expect atom
            $state = 1, $atom = '', next if $char eq '[';
            if ($char eq '%') {
                $state = 7, next;
            }
            $symbol = $char;
            push(@toks, $symbol), last unless @chars;
            $char = shift @chars;
            if ($ORGANIC_ELEMS{$symbol.$char}) {
                push @toks, $symbol.$char;
                $state = 4, next;
            } else {
                push @toks, $symbol;
                $state = 4, redo;
            }
        } elsif ($state == 4) {  # expect atom or bond or branch
            push(@toks, $char), $state = 2, next if ($char =~ /[()]/); # branch
            $state = 2, redo; # bond
        } elsif ($state == 5) {  # expect left paren
            croak "expected (" unless $char eq '(';
            $rec_smart = '';
            $paren_depth++, $state = 6, next;
        } elsif ($state == 6) {  # within recursive smarts
            $paren_depth++ if $char eq '(';
            $paren_depth-- if $char eq ')';
            unless ($paren_depth) {
                push @rec_smarts, $rec_smart;
                $state = 1, next;
            }
            $rec_smart .= $char;
        } elsif ($state == 7) {  # double digit
            $digit = $char . (shift @chars || die "expected second digit");
            push @toks, $digit;
            $state = 2;
        } else {
            die "shouldn't be here";
        }
    }
    #print Dumper \@toks if $DEBUG;
    @toks;
}


package Chemistry::File::SMARTS::simple_atom;

sub parse {
    my ($s) = shift;
    # should do some validation here...
    Chemistry::Pattern::Atom->new(symbol => $s);
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

