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
use Getopt::Long;
use File::Temp qw/ tempfile /;
use Algorithm::Combinatorics qw(combinations);


use vars qw/$LOGGER/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $blocks_file, $aligns_file);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
            "b|blocks=s"=>\$blocks_file,
            "a|alignments=s"=>\$aligns_file
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

my %read;
my %block;
open(IN, "<", $blocks_file) or $LOGGER->logdie("$! $blocks_file");
while(<IN>) {
    chomp;
    my (@f) = split(/\t/, $_);
    my $b=$f[3];
    my ($tmpfh, $tmpfilename) = tempfile();
    print { $tmpfh } $_;
    $tmpfh->close();
    $LOGGER->info("Processing block $b ...");
    if ( ! -e "block_$b.sam" ) {
        `samtools view -L $tmpfilename $aligns_file > block_$b.sam`;
    }        
    open(SAM, "<", "block_$b.sam") or $LOGGER->logdie($!);
    while(<SAM>) {
        my (@sam) = split(/\t/, $_);
        $read{$sam[0]}->{$b} = $sam[1];
        push(@{ $block{$b} }, $sam[0]);
    }
    close(SAM);
}
close(IN);

# Discarding reads mapped in more than two blocks
foreach my $r (keys %read) {
    if (scalar(keys %{ $read{$r} }) > 2) {
        delete($read{$r});
    }
}

my %count;
my $iter = combinations([keys %block], 2);

while (my $c = $iter->next) {
        $LOGGER->info("Adjacency analysis between blocks: $c->[0] and $c->[1] ...");
        
        unless ((exists $count{ $c->[0] }->{ $c->[1] }) ||
                (exists $count{ $c->[1] }->{ $c->[0] })) {
            $count{ $c->[0] }->{ $c->[1] } = 0;
            $count{ $c->[1] }->{ $c->[0] } = 0;
        }

        # bloco a
        foreach my $r (@{ $block{$c->[0]} }) {
            # se read ainda está sendo considerada
            next unless (exists $read{$r});
            
            # se estão nos dois blocos
            if ($read{$r}->{$c->[1]}) {
                my $a_revflag=reverse(sprintf("%012b", $read{$r}->{$c->[0]} ));
                my @a_flag=split(//,$a_revflag);
                my $a_type;
                # se reads mapeadas de modo apropriado
                if ($a_flag[1] == 1) {
                    # se é a primeira read R1
                    if ($a_flag[6] == 1) {
                        $a_type = 1;
                    # se é a segunda read R2                        
                    } elsif ($a_flag[7] == 1 ) {
                        $a_type = 2;
                    # caso imprevisto                        
                    } else {
                        $a_type = 0;
                    }

                    # bloco b
                    my $b_revflag=reverse(sprintf("%012b", $read{$r}->{$c->[1]} ));
                    my @b_flag=split(//,$b_revflag);
                    my $b_type;
                    # se reads mapeadas de modo apropriado
                    if ($b_flag[1] == 1) {
                        # se é a primeira read R1
                        if ($b_flag[6] == 1) {
                            $b_type = 1;
                        # se é a segunda read R2                        
                        } elsif ($b_flag[7] == 1 ) {
                            $b_type = 2;
                        # caso imprevisto                        
                        } else {
                            $b_type = 0;
                        }
                        # se a==1 e b==2 ou a==2 e b==1
                        if (($a_type)&&($b_type)&&($a_type ne $b_type)) {
                            $count{ $c->[0] }->{ $c->[1] }++;        
                            $count{ $c->[1] }->{ $c->[0] }++;
                        }
                    }
                } 
            } 
        }
}

my @order=keys %count;
print join("\t",'Block/Block', @order),"\n";
foreach my $a (@order) {
    print $a;
    foreach my $b (@order) {
        print "\t",$count{$a}->{$b}||0;
    }
    print "\n";
}


# Subroutines

sub Usage {
    my ($msg) = @_;
    my $USAGE = "
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

        -h      --help          Help
        -l      --level         Log level [Default: FATAL]
        -b      --blocks        BED file with blocks' coordinates
        -a      --alignments    BAM file of PE reads' alignments

";

    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . qq[\n];
	print STDERR $USAGE;
    exit(1);
}

