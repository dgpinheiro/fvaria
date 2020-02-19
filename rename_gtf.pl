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

my ($level, $infile, $outfile, $prefix, $delim, $b2gcsvfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "i|infile=s"=>\$infile,
            "o|outfile=s"=>\$outfile,
            "p|prefix=s"=>\$prefix,
            "d|delim=s"=>\$delim,
            "b|b2gcsvfile=s"=>\$b2gcsvfile
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
use Storable;

my $outfh;

if ($outfile) {
    $outfh = FileHandle->new(">$outfile");
    if (!defined $outfh) {
        $LOGGER->logdie("Output file ($outfile) handle not defined!");
    }
} else {
    $outfh = \*STDOUT;
}

$prefix||="";
$delim||="";

my %annot;

if ($b2gcsvfile) {
    open(B2G, "<", $b2gcsvfile) or $LOGGER->logdie($!);
    while(<B2G>) {
        chomp;
        next if ($_ =~ /^SeqName/ );
        my ($seqname, $description) = split(/\t/, $_);
        $annot{$seqname} = $description;
    }
    close(B2G);
}

my %name;
my $count=1;
open(IN, "<", $infile) or $LOGGER->logdie($!);
while(<IN>) {
    chomp;
    my @F=split(/\t/, $_);
    $F[0]=~s/^(NODE_\d+).*/$1/;
    
    my ($gene_name) = $F[8]=~/gene_name "([^"]+)"/;
    
    $LOGGER->logdie("gene_name not found here: $_") unless ($gene_name);
    
    my ($transcript_id) = $F[8]=~/transcript_id "([^"]+)"/;
    $LOGGER->logdie("transcript_id not found here: $_") unless ($transcript_id);
    
    if ($transcript_id !~ /$F[0]/) {
        $LOGGER->logdie("transcript_id ($transcript_id) doen't have seq_name ($F[0])");
    }
    
    if ($transcript_id !~ /^$gene_name/) {
        $LOGGER->logdie("transcript_id ($transcript_id) doen't begin with gene_name ($gene_name)");
    }
    
    
    my ($mrna) = $transcript_id=~/-mRNA-(\d+)/;
    $LOGGER->logdie("mRNA number not found here: $transcript_id") unless ($mrna);
    
    unless (exists $name{$gene_name}) {
        $name{$gene_name}->{'name'} = $prefix.$delim.&lpad($count,5);
        $name{$gene_name}->{'annot'} = $annot{ $transcript_id }||'';
        $count++;
    }
    unless (exists $name{$transcript_id}) {
        $name{$transcript_id}->{'annot'} = $annot{ $transcript_id }||'';
        $name{$transcript_id}->{'name'} = $name{$gene_name}->{'name'}.'-'.$mrna;
    }
    
    $F[8] =~s/transcript_id "[^"]+"/transcript_id "$name{$transcript_id}->{'name'}"/;
    $F[8] =~s/gene_name "[^"]+"/gene_name "$gene_name"/;
    $F[8] =~s/gene_id "[^"]+"/gene_id "$name{$gene_name}->{'name'}"/;
    if ($name{$gene_name}->{'annot'}) {
        $name{$gene_name}->{'annot'}=~s/\"//g;
        if ($F[8]!~/;$/) {
            $F[8].=';';
        }
        $F[8].='Note "'.$name{$gene_name}->{'annot'}.'";';
    }
#    print $transcript_id,"\t",$new_transcript_id,"\t",$gene_name,"\n";   
    print { $outfh } join("\t", @F),"\n";
}
close(IN);

$outfh->close() if (defined $outfh);

store \%name, 'rename.dump';

# Subroutines

sub Usage {
    my ($msg) = @_;
	Readonly my $USAGE => <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -i      --infile        Input file
        -o      --outfile       Output file
        -p      --prefix        Prefix [Default: BLANK]
        -d      --delim         Delimiter [Default: BLANK]
        -b      --b2gcsvfile    Blast2GO .csv file

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
