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
use FileHandle;
use POSIX 'isatty';

use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile1, $infile2);

#Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i1|infile1=s"=>\$infile1,
            "i2|infile2=s"=>\$infile2
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

$LOGGER->logdie("Missing input file 1") unless ($infile1);
$LOGGER->logdie("Wrong input file 1 ($infile1)") unless (-e $infile1);
$LOGGER->logdie("Missing input file 2") unless ($infile2);
$LOGGER->logdie("Wrong input file 2 ($infile2)") unless (-e $infile2);

my %CROSS;

my %dato;
my %gffo;
my %rowo;
my @rooto;

open(GFF1, "<$infile1") or $LOGGER->logdie($!);

my %mRNAcount;
while(<GFF1>) {
    chomp;
    next if ($_=~/^#/);
    last if ($_=~/^>/);
    my @dt=split(/\t/, $_);
    $LOGGER->logdie($_) unless ($dt[8]);
    my ($id) = $dt[8]=~/ID=([^;]+)/;
    my ($parent) = $dt[8]=~/Parent=([^;]+)/;
    
    $gffo{$id} = \@dt;
    $rowo{$id} = undef;
    if ($id =~/RNA-\d+/) {
	$mRNAcount{ $parent } = 0 unless (exists $mRNAcount{ $parent });
	$mRNAcount{ $parent }++;
	$CROSS{$dt[0]}->{$dt[2]}->{$dt[3]}->{$dt[4]}->{'o'}->{ $mRNAcount{ $parent } }=$id;
    } 

    if ($parent) {
	foreach my $p (split(/,/, $parent)) {
	        push(@{ $dato{$p} }, $id);
	}
    } else {
        push(@rooto, $id);
    }
}

close(GFF1);

my %datc;
my %gffc;
my %rowc;
my @rootc;

open(GFF2, "<$infile2") or $LOGGER->logdie($!);

while(<GFF2>) {
    chomp;
    next if ($_=~/^#/);
    last if ($_=~/^>/);
    my @dt=split(/\t/, $_);
    $LOGGER->logdie($_) unless ($dt[8]);
    my ($id) = $dt[8]=~/ID=([^;]+)/;
    my ($parent) = $dt[8]=~/Parent=([^;]+)/;
    
    $gffc{$id} = \@dt;
    $rowc{$id} = undef;
	
 	if ($id =~/mRNA:Fvar\d+-(\d+)\.2/) {
		$CROSS{$dt[0]}->{$dt[2]}->{$dt[3]}->{$dt[4]}->{'c'}->{$1}=$id;
	}

    if ($parent) {
        push(@{ $datc{$parent} }, $id);
    } else {
        push(@rootc, $id);
    }
}

close(GFF2);

print '##gff-version 3',"\n";
# gene
foreach my $id ( sort { $gffc{$a}->[0] cmp $gffc{$b}->[0] or $gffc{$a}->[3] <=> $gffc{$b}->[3] or $gffc{$a}->[4] <=> $gffc{$b}->[4] } @rootc ) {
    print join("\t", @{ $gffc{ $id } }),"\n";
    delete($rowc{$id});
    if ($datc{$id}) {
        # rna
        foreach my $subid ( sort { $gffc{$a}->[0] cmp $gffc{$b}->[0] or $gffc{$a}->[3] <=> $gffc{$b}->[3] or $gffc{$a}->[4] <=> $gffc{$b}->[4] } @{ $datc{$id} }) {
            print join("\t", @{ $gffc{ $subid } }),"\n";
            delete($rowc{$subid});
            # exon e cds
            if ($datc{ $subid }) {
		my ($i) = $subid=~/-(\d+)\.2/;
		my $subido = $CROSS{ $gffc{ $subid }->[0]  }->{ $gffc{ $subid }->[2] }->{ $gffc{ $subid }->[3] }->{ $gffc{ $subid }->[4] }->{'o'}->{$i};
		$LOGGER->logdie("$subid\t[$i]>>>NOT FOUND CROSS{ $gffc{ $subid }->[0]  }->{ $gffc{ $subid }->[2] }->{ $gffc{ $subid }->[3] }->{ $gffc{ $subid }->[4] }->{'o'}") unless ($subido);
		
		if ($dato{$subido}) {
			my $parent = $subid;
			my $c=1;
			foreach my $subsubid ( sort { $gffo{$a}->[0] cmp $gffo{$b}->[0] or $gffo{$a}->[3] <=> $gffo{$b}->[3] or $gffo{$a}->[4] <=> $gffo{$b}->[4] } @{ $dato{$subido} }) {
				if ($gffo{ $subsubid }->[2] eq 'exon') {
					my $id = 'exon:'.$parent.':'.$c++;
					$gffo{ $subsubid }->[8]=~s/ID=[^;]+;/ID=$id;/;
					$gffo{ $subsubid }->[8]=~s/Parent=[^;]+/Parent=$parent;/;
					print join("\t",@{ $gffo{ $subsubid } }),"\n";
				}
			}
		} else {
			die "Not found $subido / but found $subid";
		}
		
                foreach my $subsubid ( sort { $gffc{$a}->[0] cmp $gffc{$b}->[0] or $gffc{$a}->[3] <=> $gffc{$b}->[3] or $gffc{$a}->[4] <=> $gffc{$b}->[4] } @{ $datc{$subid} }) {
		    if ($gffc{ $subsubid }->[2] eq 'CDS') {
	                    print join("\t", @{ $gffc{ $subsubid } }),"\n";
		    }
                    delete($rowc{$subsubid});
                }
            } else {
                $LOGGER->logdie("Not found children for $subid");
            }
        }
    } else {
        $LOGGER->logdie("Not found children for $id");
    }        
}

foreach my $id (keys %rowc) {
    $LOGGER->logwarn("Not found parent/child relation for $id");
}

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
        -i1     --infile1   Input file 1 (GFF) 
        -i2     --infile2   Input file 2 (GFF) 

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

