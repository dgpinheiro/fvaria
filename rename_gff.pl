#!/usr/bin/perl
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2012  Universidade de São Paulo
#
#  Universidade de São Paulo
#  Laboratório de Biologia do Desenvolvimento de Abelhas
#  Núcleo de Bioinformática (LBDA-BioInfo)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://zulu.fmrp.usp.br/bioinfo
#
# $Id$

=head1 NAME

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Arguments:

        -h/--help   Help
        -l/--level  Log level [Default: FATAL] 
            OFF
            FATAL
            ERROR
            WARN
            INFO
            DEBUG
            TRACE
            ALL

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2012 Universidade de São Paulo

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Readonly;
use Getopt::Long;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $outfile, $prefix, $delim);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "o|outfile=s"=>\$outfile,
            "p|prefix=s"=>\$prefix,
            "d|delim=s"=>\$delim
    ) or &Usage();


if ($level) {
    my %LEVEL = (   
    'OFF'   =>$OFF,
    'FATAL' =>$FATAL,
    'ERROR' =>$ERROR,
    'WARN'  =>$WARN,
    'INFO'  =>$INFO,
    'DEBUG' =>$DEBUG,
    'TRACE' =>$TRACE,
    'ALL'   =>$ALL);
    $LOGGER->logdie("Wrong log level ($level). Choose one of: ".join(', ', keys %LEVEL)) unless (exists $LEVEL{$level});
    Log::Log4perl->easy_init($LEVEL{$level});
}

use FileHandle;

my $outfh;

if ($outfile) {
    $outfh = FileHandle->new(">$outfile");
    if (!defined $outfh) {
        $LOGGER->logdie("Output file ($outfile) handle not defined!");
    }
} else {
    $outfh = \*STDOUT;
}

use Storable;

$prefix||="";
$delim||="";

unless (-e 'rename.dump') {
    $LOGGER->logdie("Missing rename.dump");
}

my %name = %{ retrieve('rename.dump') };

open(IN, "<", $infile) or $LOGGER->logdie($!);

my %featcount;
while(<IN>) {
    chomp;
    if ($_=~/^#/) {
        print $_,"\n";
        next;
    }

    my @F=split(/\t/, $_);
    $F[0]=~s/^(NODE_\d+).*/$1/;
    
    $LOGGER->logdie(">>>$_") unless ($F[8]);

    my ($ID) = $F[8]=~/ID=([^;:]+);?/;

    if (exists $name{ $ID }) {
        no warnings;
        $F[8]=~/ID=([^;:]+)(\:(?:exon|cds|three_prime_utr|five_prime_utr)(?:\:\d+)?)?;?/;
        my $feat_id=$name{$ID}->{'name'}.($2||'');
        if (($2 eq ':cds')||($2 eq ':five_prime_utr')||($2 eq ':three_prime_utr')) {
            unless (exists $featcount{ $feat_id }) {
                $featcount{ $feat_id } = 1;
            } else {
                $featcount{ $feat_id }++;
            }
            $feat_id.=':'.$featcount{$feat_id};
        }            
        $F[8]=~s/ID=([^;:]+)(\:(?:exon|cds|three_prime_utr|five_prime_utr)(?:\:\d+)?)?;?/ID=$feat_id;/
    } else {
        $LOGGER->logdie("Not found $ID ($F[8])");
    }

    if ($F[8] =~/Parent=([^;:]+);?/) {
        my $Parent = $1;
        if (exists $name{ $Parent }) {
            $F[8]=~s/Parent=([^;:]+);?/Parent=$name{$Parent}->{'name'};/
        } else {
            $LOGGER->logdie("Not found $Parent ($F[8])");
        }
    }
    if (($F[2] ne 'exon')&&($F[2] ne 'CDS')) {
        if ( $name{$ID}->{'annot'} ) {
            $name{$ID}->{'annot'} =~ s/\"//g;
            if ($F[8]!~/;$/) {
                $F[8].=';';
            }
            $F[8] .= 'Note="' . $name{$ID}->{'annot'} . '";';
        }
    } 
    
    print { $outfh } join("\t", @F),"\n";
}
close(IN);

$outfh->close() if (defined $outfh);

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help      Help
        -l      --level     Log level [Default: FATAL]
        -i      --infile    Input file
        -o      --outfile   Output file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

sub lpad {
    my ( $num, $len ) = @_;
    return '0' x ( $len - length $num ) . $num;
}
