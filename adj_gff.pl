#!/usr/bin/env perl
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
use Bio::SeqIO;

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $emannot);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "e|emannot=s"=>\$emannot
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

$LOGGER->logdie("Missing input gff file") unless ($infile);
$LOGGER->logdie("Wrong input gff file ($infile)") unless (-e $infile);

$LOGGER->logdie("Missing emapper.py output file") unless ($emannot);
$LOGGER->logdie("Wrong emapper.py output file ($emannot)") unless (-e $emannot);

my %annot;
my %gene;
open(EM, "<", $emannot) or $LOGGER->logdie($!);
while(<EM>) {
    chomp;
    my (@c) = split(/\t/, $_);
    $annot{$c[0]} = $c[21]||'hypothetical protein';
    my $g = $c[0];
    $g=~s/mRNA://;
    $gene{$g}=$c[5];
#    print STDERR $c[0],"\t",$g,"\t",$gene{$g},"\n";
#    print STDERR $c[0],"\t",$annot{$c[0]},"\n";
}
close(EM);


my %trna;
my %AA =( 'Ala'=>0, 'Arg'=>0, 'Asn'=>0, 'Asp'=>0, 'Asx'=>0, 'Cys'=>0, 'Glu'=>0, 'Gln'=>0, 'Glx'=>0, 'Gly'=>0, 'His'=>0, 'Ile'=>0, 'Leu'=>0, 'Lys'=>0, 'Met'=>0, 'Phe'=>0, 'Pro'=>0, 'Ser'=>0, 'Thr'=>0, 'Trp'=>0, 'Tyr'=>0, 'Val'=>0);

open(GFF, "<", $infile) or $LOGGER->logdie($!);

while(<GFF>) {
	chomp;
	next if ($_=~/^#/);
	my @FIELD=split(/\t/, $_);
	if ($FIELD[2] eq 'tRNA') {
		my $trnaid;
		if ($FIELD[8]=~/Parent=([^;]+)/) {
			$trnaid=$1;
		} else {
			$LOGGER->logdie("Wrong GFF pattern: $_");
		}
		my $aminoacid;
		my $codon;
		if ($FIELD[8]=~/OldName=trnascan-[^-]+-noncoding-(\S{3})_(\S{3})/) {
			$aminoacid=$1;
			$codon=$2;
		} else {
			$LOGGER->logdie("Wrong GFF pattern: $_");
		}
		$LOGGER->logdie("Not found aminoacid $aminoacid in table ") unless (exists $AA{$aminoacid});
		$trna{$trnaid}={name=>'tRNA-'.$aminoacid,codon=>$codon};
		$AA{$aminoacid}++;
#print $trnaid,"\t",$trna{$trnaid}->{'name'},"\n";
		$FIELD[8]=~s/;$//;
		$FIELD[8].=';product='.$trna{$trnaid}->{'name'};
	} elsif ($FIELD[2] eq 'gene') {
		my ($locus_tag) = $FIELD[8]=~/ID=([^;]+)/;
		$FIELD[8]=~s/;$//;
		$FIELD[8].=';locus_tag='.$locus_tag;
        if ($gene{$locus_tag}) {
            my $name=$gene{$locus_tag};
            $FIELD[8]=~s/Name=[^;]+/Name=$name/;
        }
	} elsif ($FIELD[2] eq 'mRNA') {
		my ($transcript_id) = $FIELD[8]=~/ID=([^;]+)/;
		$FIELD[8]=~s/;$//;
		$FIELD[8].=';transcript_id=gnl|ncbi|'.$transcript_id.(($annot{$transcript_id}) ? ';product='.$annot{$transcript_id} : '');
	} elsif ($FIELD[2] eq 'exon') {
		my ($transcript_id) = $FIELD[8]=~/Parent=([^;]+)/;
		$FIELD[8]=~s/;$//;
		$FIELD[8].=';transcript_id=gnl|ncbi|'.$transcript_id;
	} elsif ($FIELD[2] eq 'CDS') {
		my ($protein_id) = $FIELD[8]=~/Parent=([^;]+)/;
		$protein_id=~s/mRNA/protein/;
		$FIELD[8]=~s/;$//;
		$FIELD[8].=';protein_id=gnl|ncbi|'.$protein_id;
	}
	print join("\t", @FIELD),"\n";
}
close(GFF);


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
        -i      --infile    Input file (.gff)
        -e      --emannot   eggNOG annotation (eggNOG-mapper output)

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

