#!/usr/bin/perl
#
#  last change Time-stamp: <2020-05-26 11:33:44 adonath>
#
#  Copyright (C) 2020  Bernhard Misof, Alexander Donath
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
#  AUTHORS: Original implementation by B. Misof
#           Reimplementation by Alexander Donath <a.donath@leibniz-zfmk.de>

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use FileHandle;
use Benchmark qw(:all);
use POSIX ();
use vars qw($USAGE);

$|++;

BEGIN {
    $USAGE = qq{
    Usage:    $0 -a alignment(s) [further options]

    Options:
	-a    Path to alignment file(s)
	      [Required]

	-s    Path to file with taxa of interest to test

	-r    Path to reference taxa file

	-k    Ensure 1KITE compatible input and output

	-m    Substitution matrix [possible values: b62|b80|b62n|b80n]
	      [Default: b62 (= original BLOSUM62 matrix)]

	      WARNING: b62n and b80n are not tested.

	-l    Name of logfile
	      [Default: log.txt]

	-o1   Name of first output file
	      [Default: outlier_1.txt]

	      Contains for each alignment a list of outlier taxa and
	      their distance score.

	-o2   Name of second output file
	      [Default: outlier_2.txt]

	      Contains for each alignment a list of outlier FASTA
	      headers.

	-h    Print this help message

	-v    Print version

    Purpose:  Identification of outliers in multiple sequence alignments.

};
}

my ( $opt_v, $opt_h, $opt_m, $opt_k, $opt_l ) = ( 0, 0, 'b62', 0, 'log.txt' );
my @opt_r;
my @opt_s;
my @opt_fl;
my $opt_o1  = 'outlier_1.txt';
my $opt_o2  = 'outlier_2.txt';
my $VERSION = 'Version 2.0.0b';

usage()
  unless GetOptions(
    'a=s{1,}' => \@opt_fl,
    'r=s{1,}' => \@opt_r,
    's=s{1,}' => \@opt_s,
    'm:s'     => \$opt_m,
    'l:s'     => \$opt_l,
    'o1:s'    => \$opt_o1,
    'o2:s'    => \$opt_o2,
    'k'       => \$opt_k,
    'v'       => \$opt_v,
    'h'       => \$opt_h
  );

die $USAGE if ($opt_h);

if ($opt_v) {
    print {*STDOUT} "checker_complete.2.pl $VERSION\n";
    exit();
}

# We need at least one alignment.
if ( !@opt_fl ) {
    print {*STDERR} "[ERROR] No alignment provided.\n";
    die $USAGE;
}

# We don't want more than one reference file
if ( @opt_r > 1 ) {
    die "[ERROR] More than one reference file provided.\n";
}

# We don't want more than one subject file
if ( @opt_s > 1 ) {
    die "[ERROR] More than one subject file provided.\n";
}

# 1KITE or default mode?
if ( !$opt_k ) {

    print {*STDOUT} "[INFO] Default mode.\n";

    # One global reference file?
    if ( @opt_r == 1 ) {

        print {*STDERR} "[INFO] Global reference file provided. @opt_r\n";

        # One specific subject file and one alignment
        if ( @opt_s == 1 && @opt_fl == 1 ) {

            print {*STDERR} "[INFO] One alignment and subject file provided. Good.\n";

        }

        # If we have more than one alignment, we need more than one subject file
        elsif ( @opt_s == 1 && @opt_fl > 1 ) {
            die "[ERROR] One reference file, more than one alignment but only one subject file provided.\n";
        }

        # No specific subject file
        elsif ( !@opt_s ) {

            print {*STDERR} "[INFO] No subject file provided. Inferring subjects from alignments automatically.\n";

        }

        # Some combination we didn't account for?
        else {
            die "som'things dodgy.\n";
        }
    }

    # No reference files provided
    else {

        print {*STDERR} "[INFO] No reference file given. Checking for *.ref files.\n";

        my $num_refs = 0;
        for my $f (@opt_fl) {

            my $fr = $f . '.refs';

            unless ( -f $fr ) {
                printf {*STDERR} "[ERROR] No file with the name %s present.\nPlease provide a reference sequence file for the alignment '%s' or a global reference file.\n",
                  $fr, $f;
                exit();
            }
            else {
                print {*STDERR} '[INFO] ', $f, ' => ', $fr, " present. Good.\n";
                $num_refs++;
            }
        }

        print {*STDOUT} "[INFO] All ($num_refs) reference file(s) present. Good.\n";

    }
}

# 1KITE mode
else {

    print {*STDOUT} "[INFO] 1kite enforced.\n";

    unless ( @opt_s == 1 && @opt_r == 0 ) {
        die "[ERROR] 1KITE mode needs a global subject file and no reference file.\n";
    }
}

#if (!$opt_r && $opt_k) {
#    print STDERR "[INFO] No explicit reference file given. 1KITE enforced. get reference taxa from files.\n";
#} elsif ($opt_r && $opt_k) {
#    print STDERR "[ERROR] Contradictory settings. Reference file given and 1kite enforced. \n";
#    die $USAGE;
#} elsif ($opt_r && !$opt_k) {
#    print STDERR "[INFO] Reference file given and 1kite not enforced.\n";
#}

if ( -f $opt_l ) {
    print {*STDERR} "[ERROR] File \'$opt_l\' already exists.\n";
    die "Stopped.\n";
} elsif (-f $opt_o1 ) {
    print {*STDERR} "[ERROR] File \'$opt_o1\' already exists.\n";
    die "Stopped.\n";
} elsif (-f $opt_o2) {
    print {*STDERR} "[ERROR] File \'$opt_o2\' already exists.\n";
    die "Stopped.\n";
}

# get starting time
my $start_time = localtime;

# print starting time
print {*STDERR} '[INFO] Starting. ', $start_time, "\n";

### MAIN ###

open my $log, '>', $opt_l;
open my $outlier1, '>', $opt_o1;
open my $outlier2, '>', $opt_o2;

print {$outlier1} 'genes with outliers', "\n";
print {$outlier2} 'genes with outliers', "\n";

my $file_counter  = 0;
my $genes_outlier = 0;

if ($opt_k) {

    # read list of taxa of interest
    # these taxa of interests have to be assembled in flat lists, for
    # example a set of our transcriptome taxa could be in this list
    # for 1KITE-compatible output, this corresponds to the taxon name
    # between the second an the third pipe symbol "|"
    my @subjects = &arrayFile( $opt_s[0] );

    for my $file (@opt_fl) {

        unless ($file) {
            print {*STDERR} "[ERROR] Couldn\'t open file \'$file\'.\n";
            die $USAGE;
        }

        # read fasta file
        my ( $ref_FASTA, $ref_al_length ) = &hashFasta($file);
        $file_counter++;

        # local list of taxa of interest (transcriptomes) as keys and
        # specific reference taxa to this taxon of interest as value
        my %subject_local = ();

        # local list of all reference taxa from which median BLOSUM
        # distances will be calculated
        my @reference_taxa = ();

        # local reference taxa drawn from the fasta file to which
        # reference taxa should be compared in BLOSUM distances
        # OBSOLET: always a subset of the @reference_taxa list
        my @query_local = ();

        my $ntaxa = keys %$ref_FASTA;

        my $logtext = sprintf "file: %-20.20s\talignment length: %-10.10s\n",
          $file,
          $$ref_al_length;

        print {*STDOUT} $logtext;
        print {$log} $logtext;

        # pushes all taxa which do not belong to the list of taxa of
        # interest into a reference list, with this list, a median
        # distance value and variance of the expected BLOSUM distances
        # will be calculated
        for my $reference_taxa ( keys %$ref_FASTA ) {
            push @reference_taxa, $reference_taxa
              if 2 == $reference_taxa =~ tr/\|/\|/;
        }

        # reads full names of taxa of interest from the fasta file and
        # stores these full names in a hash as keys only for the present
        # file with the name of the reference taxon to be compared with as
        # value of this key, if a taxon of interest is not present in the
        # fasta file, it will not be stored in the subject_local file,
        # nothing will be stored for this taxon, also no key
        @query_local = &localQueryList( $ref_FASTA, \@reference_taxa );
        %subject_local = &local1KiteSubjectList( $ref_FASTA, \@subjects, \@query_local );

        my $transcriptomes = $ntaxa - @query_local;

        for my $taxon (@subjects) {
            if ( !grep /$taxon/, keys %subject_local ) {
                $logtext = sprintf "%-30.30s\t%-10.10s\n", $taxon, 'missing';
                print {*STDOUT} $logtext;
                print {$log} $logtext;
            }
        }

        $logtext = sprintf "%-30.30s\t%-10.10s\n", 'Number of taxa in alignment:', $ntaxa;

        print {*STDOUT} $logtext;
        print {$log} $logtext;

        $logtext = sprintf "%-30.30s\t%-10.10s\n", 'Transcriptomes:',
          $transcriptomes;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        # projects coordinates of taxon of interest onto reference taxon,
        # goes through all taxa of interests and derives the section of
        # consensus overlap for all taxa of interest; corresponds to the
        # maximal extent of contigs only for this section of the
        # alignment, distance calculations for all comparisons will be
        # performed, this sequence length should be given in the log file
        # and on the screen
        my ( $start, $end ) = ( $$ref_al_length, 0 );

        for my $taxon ( keys %subject_local ) {
            my ( $s, $e ) = &indices( $taxon, $ref_al_length, $ref_FASTA );
            $start = $s < $start ? $s : $start;
            $end   = $e > $end   ? $e : $end;
        }

        $logtext = sprintf "%-30.30s\t%-10.10s\n\n", 'overlap length: ',
          $end - $start;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        $logtext = sprintf "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n", 'min',
          'max', 'median', 'Q1', 'Q3';
        print {*STDOUT} $logtext;
        print {$log} $logtext;

# calculates BLOSUM distances for all reference taxa in the alignment within the
# section of overlap; indel and X positions are ignored in each pair separately!
# works for nucleotide and aminoacid data in case of nucleotide data
# match/missmatch distances are calculated
        my $ref_scoring     = &BLOSUM_scoring($opt_m);
        my %BLOSUM_dist_of  = ();
        my @dissimilarities = ();

        for (@query_local) {
            $BLOSUM_dist_of{$_} = 1000;
        }

        while ( 1 < @query_local ) {

            my $taxon    = shift @query_local;
            my $dist_min = 1000;

            for my $r (@query_local) {

                my $dist = &BLOSUM_dist( $taxon, $r, $opt_m, [ $start .. $end ],
                    $ref_FASTA, $ref_scoring );
                push @dissimilarities, $dist;
                $dist_min           = $dist if $dist < $dist_min;
                $BLOSUM_dist_of{$r} = $dist if ( $BLOSUM_dist_of{$r} > $dist );
            }
            $BLOSUM_dist_of{$taxon} = $dist_min
              if ( $BLOSUM_dist_of{$taxon} > $dist_min );
        }

        # uses distances to calculate median, min, max, and quartiles
        @dissimilarities = sort { $a <=> $b } @dissimilarities;

        my $median_dissi = &quantile( 2, \@dissimilarities );
        $median_dissi = 1 if $median_dissi == 0;

        my $dist_min = $dissimilarities[0];
        my $dist_max = $dissimilarities[-1];

        my $Q1 = &quantile( 1, \@dissimilarities );
        my $Q3 = &quantile( 3, \@dissimilarities );

        $logtext = sprintf "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n\n",
          $dist_min, $dist_max, $median_dissi, $Q1, $Q3;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        # calculates outlier cutoff from reftaxa distances the upper
        # whisker of reftaxa distances
        my $outlier_cutoff = &outlier( \@dissimilarities );

        $logtext = sprintf "%-30.30s\t%-10.10s\n\n",
          'Transcriptome outlier dist: >', $outlier_cutoff;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        $logtext = sprintf "%-30.30s\t%-6.6s\t%-6.6s\n", 'taxon', 'dist',
          're.dist';
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        if ( $dist_max > $outlier_cutoff ) {
            $logtext = "some reftaxon distances outlier\n";
            print {*STDOUT} $logtext;
            print {$log} $logtext;
        }

        # calculates BLOSUM distance for closest reftaxon for each transcriptome
        # contig only for consensus
        # overlap length of the alignment and collects them in a hash
        for my $taxon ( keys %subject_local ) {

            $BLOSUM_dist_of{$taxon} =
              &BLOSUM_dist( $taxon, $subject_local{$taxon}, $opt_m,
                [ $start .. $end ],
                $ref_FASTA, $ref_scoring );
        }

        my ( $flag_outlier, $count_outlier );
        ( $flag_outlier, $count_outlier ) = &print1KiteResults(
            $outlier_cutoff, $dist_min, $dist_max, $median_dissi,
            $Q1,             $Q3,       $file,     $log,
            $outlier1,       $outlier2, %BLOSUM_dist_of
        );

        $genes_outlier++ if $flag_outlier == 1;
        $logtext = sprintf "\noutliers counted: %i\n\n", $count_outlier;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

    }
}

# Default mode
else {

    # read list of taxa of interest
    # these taxa of interests have to be assembled in flat lists, for
    # example a set of our transcriptome taxa could be in this list
    # for 1KITE-compatible output, this corresponds to the taxon name
    # between the second an the third pipe symbol "|"
    my @subjects = ();

    # local list of taxa of interest (transcriptomes) as keys and
    # specific reference taxa to this taxon of interest as value
    my %subject_local = ();

    # One global subject list
    my $globalsubject = 0;
    if ( @opt_s == 1 ) {
        @subjects      = &arrayFile( $opt_s[0] );
        $globalsubject = 1;
    }

    # local list of all reference taxa from which median BLOSUM
    # distances will be calculated
    my @reference_taxa = ();

    # pushes all taxa which do not belong to the list of taxa of
    # interest into a reference list, with this list, a median
    # distance value and variance of the expected BLOSUM distances
    # will be calculated
    if ( @opt_r == 1 ) {
        @reference_taxa = &arrayFile( $opt_r[0] );
    }

    for my $file (@opt_fl) {

        unless ($file) {
            print {*STDERR} "[ERROR] Couldn\'t open file \'$file\'.\n";
            die $USAGE;
        }

        # read fasta file
        my ( $ref_FASTA, $ref_al_length ) = &hashFasta($file);
        $file_counter++;

        # local reference taxa drawn from the fasta file to which
        # reference taxa should be compared in BLOSUM distances
        # OBSOLET: always a subset of the @reference_taxa list
        my @query_local = ();

        my $ntaxa = keys %$ref_FASTA;

        my $logtext = sprintf "file: %-20.20s\talignment length: %-10.10s\n",
          $file,
          $$ref_al_length;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        if ( @opt_r == 0 ) {
            my $tmp_ref_file = $file . '.refs';
            @reference_taxa = &arrayFile($tmp_ref_file);
        }

        # we need this here already
        my $ref_scoring = &BLOSUM_scoring($opt_m);

        # reads full names of taxa of interest from the fasta file and
        # stores these full names in a hash as keys only for the present
        # file with the name of the reference taxon to be compared with as
        # value of this key, if a taxon of interest is not present in the
        # fasta file, it will not be stored in the subject_local file,
        # nothing will be stored for this taxon, also no key
        @query_local = &localQueryList( $ref_FASTA, \@reference_taxa );

# A custom subject list, infered from the alignment, based on the reference file(s).
        if ( !@opt_s ) {
            %subject_local =
              &localSubjectList( $ref_FASTA, $ref_al_length, $opt_m,
                $ref_scoring, \@query_local );

            for my $subject ( keys %subject_local ) {
                push @subjects, $subject;
            }

        }
        elsif ( $globalsubject && @opt_s == 1 ) {
            %subject_local = &globalSubjectList(
                $ref_FASTA,   $ref_al_length, $opt_m,
                $ref_scoring, \@subjects,     \@query_local
            );
        }

        my $transcriptomes = keys %subject_local;

        for my $taxon (@subjects) {
            if ( !grep /\Q$taxon\E/, keys %subject_local ) {
                $logtext = sprintf "%-30.30s\t%-10.10s\n", $taxon, 'missing';
                print {*STDOUT} $logtext;
                print {$log} $logtext;
            }
        }

        $logtext = sprintf "%-30.30s\t%-10.10s\n",
          'Number of taxa in alignment:', $ntaxa;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        $logtext = sprintf "%-30.30s\t%-10.10s\n", 'Transcriptomes:',
          $transcriptomes;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        # projects coordinates of taxon of interest onto reference taxon,
        # goes through all taxa of interests and derives the section of
        # consensus overlap for all taxa of interest; corresponds to the
        # maximal extent of contigs only for this section of the
        # alignment, distance calculations for all comparisons will be
        # performed, this sequence length should be given in the log file
        # and on the screen
        my ( $start, $end ) = ( $$ref_al_length, 0 );

        for my $taxon ( keys %subject_local ) {
            my ( $s, $e ) = &indices( $taxon, $ref_al_length, $ref_FASTA );
            $start = $s < $start ? $s : $start;
            $end   = $e > $end   ? $e : $end;
        }

        $logtext = sprintf "%-30.30s\t%-10.10s\n\n", 'overlap length: ',
          $end - $start;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        $logtext = sprintf "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n", 'min',
          'max', 'median', 'Q1', 'Q3';
        print {*STDOUT} $logtext;
        print {$log} $logtext;

# calculates BLOSUM distances for all reference taxa in the alignment within the
# section of overlap; indel and X positions are ignored in each pair separately!
# works for nucleotide and aminoacid data in case of nucleotide data
# match/missmatch distances are calculated
        my %BLOSUM_dist_of  = ();
        my @dissimilarities = ();

        for (@query_local) {
            $BLOSUM_dist_of{$_} = 1000;
        }

        while ( 1 < @query_local ) {

            my $taxon    = shift @query_local;
            my $dist_min = 1000;

            for my $r (@query_local) {

                my $dist = &BLOSUM_dist( $taxon, $r, $opt_m, [ $start .. $end ],
                    $ref_FASTA, $ref_scoring );
                push @dissimilarities, $dist unless ( $dist == 1000 );   #modded
                $dist_min           = $dist if $dist < $dist_min;
                $BLOSUM_dist_of{$r} = $dist if ( $BLOSUM_dist_of{$r} > $dist );
            }
            $BLOSUM_dist_of{$taxon} = $dist_min
              if ( $BLOSUM_dist_of{$taxon} > $dist_min );
        }

        # uses distances to calculate median, min, max, and quartiles
        @dissimilarities = sort { $a <=> $b } @dissimilarities;

        my $median_dissi = &quantile( 2, \@dissimilarities );
        $median_dissi = 1 if $median_dissi == 0;

        my $dist_min = $dissimilarities[0];
        my $dist_max = $dissimilarities[-1];

        my $Q1 = &quantile( 1, \@dissimilarities );
        my $Q3 = &quantile( 3, \@dissimilarities );

        $logtext = sprintf "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n\n",
          $dist_min, $dist_max, $median_dissi, $Q1, $Q3;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        # calculates outlier cutoff from reftaxa distances the upper
        # whisker of reftaxa distances
        my $outlier_cutoff = &outlier( \@dissimilarities );

        $logtext = sprintf "%-30.30s\t%-10.10s\n\n",
          'Transcriptome outlier dist: >', $outlier_cutoff;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        $logtext = sprintf "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n", 'out?', 'dist',
          're.dist', 'taxon';
        print {*STDOUT} $logtext;
        print {$log} $logtext;

        #	if ($dist_max > $outlier_cutoff) {
        #	    $logtext = "some reftaxon distances outlier\n";
        #	    print STDOUT $logtext;
        #	    print {$log} $logtext;
        #	}

        # calculates BLOSUM distance for closest reftaxon for each transcriptome
        # contig only for consensus
        # overlap length of the alignment and collects them in a hash
        for my $taxon ( keys %subject_local ) {
            $BLOSUM_dist_of{$taxon} =
              &BLOSUM_dist( $taxon, $subject_local{$taxon}, $opt_m,
                [ $start .. $end ],
                $ref_FASTA, $ref_scoring );
        }

        my ( $flag_outlier, $count_outlier );
        ( $flag_outlier, $count_outlier ) = &printResults(
            $outlier_cutoff, $dist_min, $dist_max, $median_dissi,
            $Q1,             $Q3,       $file,     $log,
            $outlier1,       $outlier2, %BLOSUM_dist_of
        );

        $genes_outlier++ if $flag_outlier == 1;
        $logtext = sprintf "\noutliers counted: %i\n\n", $count_outlier;
        print {*STDOUT} $logtext;
        print {$log} $logtext;

    }
}

print {$log}
  "\n$file_counter files (genes) checked\n$genes_outlier files with outliers\n";
close $log;
close $outlier1;
close $outlier2;

print <<STAT;

total number of files checked: $file_counter
files with outliers:           $genes_outlier

STAT

# get ending time
my $end_time = localtime;

# print ending time
print {*STDERR} '[INFO] Done. ', $end_time, "\n";

### SUBS ###

# Title    : hashFasta
# Usage    :
# Function : read in a fasta file as hash
# Args     : file = the fasta file
# Returns  : hash containing the file
sub hashFasta {

    my ($file) = @_;
    my (%fastafile);
    my ( $id, $description ) = ( '', '' );
    open (my $fh, '<', $file) or die "[ERROR] Couldn\'t open \'$file\'.\n", $USAGE;

    # read file and store in hash
    while (<$fh>) {
        my $line = $_;
        chomp $line;

        if ( $line =~ /^>(\S+)(.*)$/ ) {
            $id                           = $1;
            $description                  = $2;
            $fastafile{"$id$description"} = '';
        }
        elsif ( $line =~ /^\n/ ) {
            next;
        }
        elsif ( $line =~ /^\s/ ) {
            next;
        }
        else {
            $fastafile{"$id$description"} .= uc($line);
        }
    }
    close $fh;

    # check whether all sequences have an equal length
    my $aln_length = -1;
    for my $seqs ( values %fastafile ) {
        if ( $aln_length == -1 ) {
            $aln_length = length $seqs;
            next;
        }
        elsif ( $aln_length != length $seqs ) {
            die "[ERROR] Sequences of unequal length.\n";
        }
    }

    # mask leading and trailing gaps
    for my $seqs ( values %fastafile ) {
        for ( 0 .. $aln_length - 1 ) {
            substr( $seqs, $_, 1 ) =~ /\-/ ? substr( $seqs, $_, 1, 'X' ) : last;
        }
        for ( 0 .. $aln_length - 1 ) {
            substr( $seqs, $aln_length - ( 1 + $_ ), 1 ) =~ /\-/
              ? substr( $seqs, $aln_length - ( 1 + $_ ), 1, 'X' )
              : last;
        }
    }

    return ( \%fastafile, \$aln_length );
}

# Title    : arrayFile
# Usage    :
# Function : Read in a file as an array. Duplicated lines are not allowed.
#            The i'th non-empty line is stored in $a[i].
# Args     : file = the file
# Returns  : reference to array containing the file
sub arrayFile {

    my ($file) = @_;

    my @array = ();
    my %seen  = ();

    open (my $fh, '<', $file) or die "[ERROR] Couldn\'t open \'$file\'.\n", $USAGE;

    # read file and store in hash
    while (<$fh>) {

        my $line = $_;
        chomp $line;

        if ( $line =~ /\S+/ ) {

            die "[ERROR] Duplicate entry in \'$file\'.\n"
              if defined $seen{$line};
            $seen{$line}++;
            push @array, $line;

        }
        else {
            next;
        }
    }

    close $fh;

    return @array;
}

# Title    : local1KiteSubjectList
# Usage    :
# Function :
# Args     : ref_FASTA =
#            subjects =
#            reference_taxa =
# Returns  : hash containing subjects as key and respective reference
#            taxon as value
sub local1KiteSubjectList {

    my ( $ref_FASTA, $subjects, $reference_taxa ) = @_;

    my %hash = ();

    for my $taxon (@$subjects) {

        for my $k ( keys %$ref_FASTA ) {

            next if 2 == $k =~ tr/\|/\|/;    #ignore reference species

            my @array = split /\Q|\E/, $k;

            if ( $taxon eq $array[2] ) {

                my $taxon_full = $k;

                for my $k2 ( keys %$ref_FASTA ) {

                    if ( 2 == $k2 =~ tr/\|/\|/ and $k2 =~ /$array[1]/ ) {

                        push @$reference_taxa, $k2
                          if !grep /$array[1]/, @$reference_taxa;
                        $hash{$taxon_full} = $k2;
                        last;
                    }
                }
            }
        }
    }

    return %hash;
}

# Title    : localSubjectList
# Usage    :
# Function :
# Args     : ref_FASTA =
#            ref_al_length =
#            matrix =
#            ref_scoring =
#            references =
# Returns  : hash containing subjects as key and closest taxon as value
sub localSubjectList {

    my ( $ref_FASTA, $ref_al_length, $matrix, $ref_scoring, $references ) = @_;

    my %hash = ();

    my ( $start, $end ) = ( $ref_al_length, 0 );

    my %reference_taxa = ();
    $reference_taxa{$_}++ for (@$references);

    my %subject_taxa = ();

    for my $taxon ( keys %$ref_FASTA ) {
        if ( !exists $reference_taxa{$taxon} && defined $$ref_FASTA{$taxon} ) {
            $subject_taxa{$taxon}++;
        }
        else {
            print {*STDERR} "[WARN] Subject taxon $taxon not present in Alignment.\n";
        }
    }

    for my $taxon ( keys %subject_taxa ) {
        my ( $s, $e ) = &indices( $taxon, $ref_al_length, $ref_FASTA );
        $start = $s < $start ? $s : $start;
        $end   = $e > $end   ? $e : $end;
    }

    for my $taxon ( keys %subject_taxa ) {

        my $dist_min = 1001;
        my $ref_min;

        for my $current_reference (@$references) {

            next if ( $taxon eq $current_reference );

            my $dist = &BLOSUM_dist(
                $taxon,             $current_reference, $matrix,
                [ $start .. $end ], $ref_FASTA,         $ref_scoring
            );

            if ( $dist < $dist_min ) {
                $dist_min = $dist;
                $ref_min  = $current_reference;
            }
        }

        $hash{$taxon} = $ref_min;
    }

    return %hash;
}

# Title    : globalSubjectList
# Usage    :
# Function :
# Args     : ref_FASTA =
#            ref_al_length =
#            matrix =
#            ref_scoring =
#            subjects =
#            references =
# Returns  : hash containing subjects as key and closest taxon as value
sub globalSubjectList {

    my (
        $ref_FASTA,   $ref_al_length, $matrix,
        $ref_scoring, $subjects,      $references
    ) = @_;

    my %hash = ();

    my ( $start, $end ) = ( $ref_al_length, 0 );

    my %subject_taxa = ();

    for my $taxon (@$subjects) {
        if ( defined $$ref_FASTA{$taxon} ) {
            $subject_taxa{$taxon}++;
        }
        else {
            print {*STDERR} "[WARN] Subject taxon $taxon not present in Alignment.\n";
        }
    }

    for my $taxon ( keys %subject_taxa ) {
        my ( $s, $e ) = &indices( $taxon, $ref_al_length, $ref_FASTA );
        $start = $s < $start ? $s : $start;
        $end   = $e > $end   ? $e : $end;
    }

    for my $taxon ( keys %subject_taxa ) {

        my $dist_min = 1001;
        my $ref_min;

        for my $current_reference (@$references) {

            next if ( $taxon eq $current_reference );

            my $dist = &BLOSUM_dist(
                $taxon,             $current_reference, $matrix,
                [ $start .. $end ], $ref_FASTA,         $ref_scoring
            );

            if ( $dist < $dist_min ) {
                $dist_min = $dist;
                $ref_min  = $current_reference;
            }
        }

        $hash{$taxon} = $ref_min;
    }

    return %hash;
}

# Title    : localQueryList
# Usage    :
# Function :
# Args     : ref_FASTA =
#            list =
# Returns  : array containing all present reference header
sub localQueryList {

    my ( $ref_FASTA, $list ) = @_;

    my @array = ();

    for my $taxon (@$list) {
        if ( defined $$ref_FASTA{$taxon} ) {
            push @array, $taxon;
        }
        else {
            print {*STDERR} "[WARN] Reference taxon $taxon not present in Alignment.\n";
        }
    }

    return @array;
}

# Title    : quantile
# Usage    :
# Function : Calculates the ith quantile using quantile algorithm Q[7]
#            as defined in http://doi.org/10.2307/2684934
#            Q[7] is the default algorithm in R, Excel and Statistics::Descriptive
# Args     : quantile = an integer between 0 and 4, where
#              Q0 = min value,
#              Q1 = 25th percentile (lower quartile),
#              Q2 = 50th percentile (median),
#              Q3 = 75th percentile (upper quartile), and
#              Q4 = max value.
#            array = array reference which contains sorted elements
# Returns  : the ith quantile of array
sub quantile {

    my ( $quantile, $array ) = @_;

    my $k = ( $quantile / 4 ) * ( @$array - 1 ) + 1;
    my $f = $k - POSIX::floor($k);
    $k = POSIX::floor($k);

    my $kq = @$array[ $k - 1 ];
    return $kq if ( $f == 0 );

    my $kpq = @$array[$k];

    my $q = $kq + ( $f * ( $kpq - $kq ) );

    return $q;
}

# Title    : outlier
# Usage    :
# Function : Calculates ...
# Args     : array = array reference which contains sorted elements
# Returns  : upper whisker
sub outlier {

    my ($ref_array) = @_;

    my ( $q1, $median, $q3, $w1, $w3 );

    @$ref_array = sort { $a <=> $b } @$ref_array;

    my $dist_max = $$ref_array[-1];

    $q1     = &quantile( 1, $ref_array );
    $q3     = &quantile( 3, $ref_array );
    $median = &quantile( 2, $ref_array );

    my $iqr = $q3 - $q1;
    my $aw1 = $q1 - ( 1.5 * $iqr );
    my $aw3 = $q3 + ( 1.5 * $iqr );

    return ($aw3);
}

# Title    : BLOSUM_scoring
# Usage    :
# Function :
# Args     : matrix = the requested substitution matrix
# Returns  : the scoring matrix
sub BLOSUM_scoring {

    my ($matrix) = @_;

    my %scoring;

    my ( $blosum, $aminoacids ) = &blosumMatrix($matrix);

    my $aa_lead = 0;
    for my $line (@$blosum) {

        for ( my $i = 0 ; $i <= $#{$line} ; $i++ ) {

            $scoring{ $$aminoacids[$aa_lead] . $$aminoacids[$i] } = $line->[$i];

        }
        $aa_lead++;
    }

    my $indel_score = $scoring{ 'A' . '*' };
    for (@$aminoacids) {
        $scoring{ $_ . '-' }  = $indel_score;
        $scoring{ '-' . $_ }  = $indel_score;
        $scoring{ '-' . '-' } = $indel_score;
    }

    return \%scoring;
}

# Title    : printHash
# Usage    :
# Function :
# Args     : hash = the hash that should be printed
# Returns  : nothing
sub printHash {

    my ($hash) = @_;

    foreach my $key ( sort keys %$hash ) {
        print {*STDOUT} "$key : $hash->{$key}\n";
    }
}

# Title    : print1KiteResults
# Usage    :
# Function :
# Args     : outlier_cutoff
#            dist_min
#	     	 dist_max
#	     	 median_dissi
#	     	 Q1
#            Q3
#	     	 file
#	     	 log
#	     	 outlier1
#	     	 outlier2
#	     	 BLOSUM_dist_of
# Returns  : nothing
sub print1KiteResults {

    my (
        $outlier_cutoff, $dist_min, $dist_max, $median_dissi,
        $Q1,             $Q3,       $file,     $log,
        $outlier1,       $outlier2, %BLOSUM_dist_of
    ) = @_;

    my $flag_outlier  = 0;
    my $count_outlier = 0;
    my $logtext;

    # checks whether BLOSUM distance of transcriptome to closest
    # reftaxon is above outlier cutoff
    for my $taxon ( keys %BLOSUM_dist_of ) {

        my @taxon = split /\Q|\E/, $taxon;

        if (    ( $BLOSUM_dist_of{$taxon} > $outlier_cutoff )
            and ( $flag_outlier == 0 ) )
        {

            $flag_outlier = 1;
            $count_outlier++;

            printf {$outlier1} "\n$file:\n";
            printf {$outlier2} "\n$file:\n";
            printf {$outlier1} "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n",
              'min', 'max', 'median', 'Q1', 'Q3';
            printf {$outlier1} "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n\n",
              $dist_min, $dist_max, $median_dissi, $Q1, $Q3;

            if ( 4 == @taxon ) {

                $logtext = sprintf "%-30.30s\t%-6.6s\t%-6.6s\t%-6.6s\n",
                  "$taxon[2] = ",
                  $BLOSUM_dist_of{$taxon},
                  $BLOSUM_dist_of{$taxon} / $median_dissi,
                  'out!';

                print {*STDOUT} $logtext;
                print {$log} $logtext;
                print {$outlier1} $logtext;
                print {$outlier2} "$taxon\n";
            }
            else {
                $logtext = sprintf "%-30.30s\t%-6.6s\t%-6.6s\t%-6.6s\n",
                  "$taxon[1] = ",
                  $BLOSUM_dist_of{$taxon},
                  $BLOSUM_dist_of{$taxon} / $median_dissi,
                  'out!';
                print {*STDOUT} $logtext;
                print {$log} $logtext;
                print {$outlier1} $logtext;
                print {$outlier2} "$taxon\n";
            }
        }
        elsif ( ( $BLOSUM_dist_of{$taxon} > $outlier_cutoff )
            and ( $flag_outlier == 1 ) )
        {

            $count_outlier++;

            if ( 4 == @taxon ) {
                $logtext = sprintf "%-30.30s\t%-6.6s\t%-6.6s\t%-6.6s\n",
                  "$taxon[2] = ",
                  $BLOSUM_dist_of{$taxon},
                  $BLOSUM_dist_of{$taxon} / $median_dissi,
                  'out!';
                print {*STDOUT} $logtext;
                print {$log} $logtext;
                print {$outlier1} $logtext;
                print {$outlier2} "$taxon\n";
            }
            else {
                $logtext = sprintf "%-30.30s\t%-6.6s\t%-6.6s\t%-6.6s\n",
                  "$taxon[1] = ",
                  $BLOSUM_dist_of{$taxon},
                  $BLOSUM_dist_of{$taxon} / $median_dissi,
                  'out!';
                print {*STDOUT} $logtext;
                print {$log} $logtext;
                print {$outlier1} $logtext;
                print {$outlier2} "$taxon\n";
            }
        }
        else {
            if ( @taxon == 4 ) {
                $logtext = sprintf "%-30.30s\t%-6.6s\t%-6.6s\n",
                  "$taxon[2] = ",
                  $BLOSUM_dist_of{$taxon},
                  $BLOSUM_dist_of{$taxon} / $median_dissi;
                print {*STDOUT} $logtext;
                print {$log} $logtext;
            }
            else {
                $logtext = sprintf "%-30.30s\t%-6.6s\t%-6.6s\n",
                  "$taxon[1] = ",
                  $BLOSUM_dist_of{$taxon},
                  $BLOSUM_dist_of{$taxon} / $median_dissi;
                print {*STDOUT} $logtext;
                print {$log} $logtext;
            }
        }
    }

    return ( $flag_outlier, $count_outlier );
}

# Title    : printResults
# Usage    :
# Function :
# Args     : outlier_cutoff
#            dist_min
#	     	 dist_max
#	     	 median_dissi
#	     	 Q1
#            Q3
#	     	 file
#	     	 log
#	     	 outlier1
#	     	 outlier2
#	     	 BLOSUM_dist_of
# Returns  : nothing
sub printResults {

    my (
        $outlier_cutoff, $dist_min, $dist_max, $median_dissi,
        $Q1,             $Q3,       $file,     $log,
        $outlier1,       $outlier2, %BLOSUM_dist_of
    ) = @_;

    my $flag_outlier  = 0;
    my $count_outlier = 0;
    my $logtext;

    # checks whether BLOSUM distance of transcriptome to closest
    # reftaxon is above outlier cutoff
    for my $taxon ( keys %BLOSUM_dist_of ) {

        #my @taxon = split/\Q|\E/,$taxon;

        if (    ( $BLOSUM_dist_of{$taxon} > $outlier_cutoff )
            and ( $flag_outlier == 0 ) )
        {

            $flag_outlier = 1;
            $count_outlier++;

            printf {$outlier1} "\n$file:\n";
            printf {$outlier2} "\n$file:\n";
            printf {$outlier1} "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n",
              'min', 'max', 'median', 'Q1', 'Q3';
            printf {$outlier1} "%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\t%-6.6s\n\n",
              $dist_min, $dist_max, $median_dissi, $Q1, $Q3;

            $logtext = sprintf "%-6.6s\t%-6.6s\t%-6.6s\t%s\n",
              'out!',
              $BLOSUM_dist_of{$taxon},
              $BLOSUM_dist_of{$taxon} / $median_dissi,
              " = $taxon";

            print {*STDOUT} $logtext;
            print {$log} $logtext;
            print {$outlier1} $logtext;
            print {$outlier2} "$taxon\n";
        }
        elsif ( ( $BLOSUM_dist_of{$taxon} > $outlier_cutoff )
            and ( $flag_outlier == 1 ) )
        {

            $count_outlier++;

            $logtext = sprintf "%-6.6s\t%-6.6s\t%-6.6s\t%s\n",
              'out!',
              $BLOSUM_dist_of{$taxon},
              $BLOSUM_dist_of{$taxon} / $median_dissi,
              " = $taxon";
            print {*STDOUT} $logtext;
            print {$log} $logtext;
            print {$outlier1} $logtext;
            print {$outlier2} "$taxon\n";

        }
        else {
            $logtext = sprintf "%14.6s\t%-6.6s\t%s\n",
              $BLOSUM_dist_of{$taxon},
              $BLOSUM_dist_of{$taxon} / $median_dissi,
              " = $taxon";
            print {*STDOUT} $logtext;
            print {$log} $logtext;
        }
    }

    return ( $flag_outlier, $count_outlier );
}

# Title    : expectedScore
# Usage    :
# Function :
# Args     : matrix = the requested substitution matrix
# Returns  : the expected score for scaling
sub expectedScore {

    my ($matrix) = @_;

    for ($matrix) {

        if    (/^b62$/)  { return -0.5209 }
        elsif (/^b80$/)  { return -0.7442 }
        elsif (/^b62n$/) { return -0.5209 }               # probably not correct
        elsif (/^b80n$/) { return -0.7442 }               # probably not correct
        else             { die '[ERROR] Unknown matrix.' }
    }
}

# Title    : BLOSUM_dist
# Usage    :
# Function : Calculates BLOSUM distance according to the Scoredist approach of
#            Sonnhammer & Hollich BMC Bioinformatics 2005, 6:108
#            doi:10.1186/1471-2105-6-108.
#            Distances are not calibrated like in Sonnhammer & Hollich 2005.
# Args     :
# Returns  : The raw (uncalibrated) BLOSUM distance of two sequences.
#            If the overlap is less than 30 residues or in case of model violation,
#            the BLOSUM distance will be set to 1000. This prevents a strong
#            influence of only short overlaps between proteins.
sub BLOSUM_dist {

    my ( $taxon, $query, $scoreMatrix, $ref_range, $ref_FASTA, $ref_scoring ) =
      @_;

    my $overlap        = @$ref_range;
    my $score          = 0;
    my $score_taxon    = 0;
    my $score_query    = 0;
    my $expected_score = 0;
    my $BLOSUM_dist    = 0;

    my ( $overlap_I, $identical ) =
      &compare( $taxon, $query, $ref_range, $ref_FASTA );

    # No overlap.
    if ( $overlap_I == 0 ) {

        print {*STDERR} "[WARN] No sequence overlap between $taxon and $query. Setting BLOSUM distance to 1000.\n";
        return $BLOSUM_dist = 1000;
    }

    # Multiply overlap with scaling factor.
    $expected_score = &expectedScore($scoreMatrix) * $overlap_I;

    # Expected score shouldn't be positive.
    if ( 0 <= $expected_score ) {

        print {*STDERR} "[WARN] Expected score shouldn't be positive. Setting BLOSUM distance to 1000.\n";
        return $BLOSUM_dist = 1000;
    }

    for (@$ref_range) {
        if (   substr( $ref_FASTA->{$taxon}, $_, 1 ) =~ /X|\-|\*/
            or substr( $ref_FASTA->{$query}, $_, 1 ) =~ /X|\-|\*/ )
        {

            $overlap-- and next;

        }

        # sigma(s1,s2)
        $score += $ref_scoring->{ substr( $ref_FASTA->{$taxon}, $_, 1 )
              . substr( $ref_FASTA->{$query}, $_, 1 ) };

        # sigma(s1,s1)
        $score_taxon += $ref_scoring->{ substr( $ref_FASTA->{$taxon}, $_, 1 )
              . substr( $ref_FASTA->{$taxon}, $_, 1 ) };

        # sigma(s2,s2)
        $score_query += $ref_scoring->{ substr( $ref_FASTA->{$query}, $_, 1 )
              . substr( $ref_FASTA->{$query}, $_, 1 ) };
    }

    # Expected score for two random sequences.
    # Null model. = lower limit.
    # sigma^r(l) = sigma(0) * l
    $expected_score = &expectedScore($scoreMatrix) * $overlap;

    # Expected score between two random sequences shouldn't be larger
    # than score between two non-random sequences, i.e.
    # sigma(s1,s2) < sigma^r(l) is a violation of the lower boundary.
    if ( $score <= $expected_score ) {

        print {*STDERR} "[WARN] $taxon VERSUS $query\n";
        print {*STDERR} "[WARN] Expected score of two random sequences shouldn't be larger than expected\n";
        print {*STDERR} "[WARN] score of two non-random sequences. Setting BLOSUM distance to 1000.\n";

        return $BLOSUM_dist = 1000;
    }

    # normalized score.
    # sigma(N) = sigma(s1,s2) - sigma^r(l)
    my $normalized_score = $score - $expected_score;

    # normalized upper score.
    # sigma(UN) = sigma(U) - sigma^r(l)
    my $normalized_upper_score =
      ( $score_taxon + $score_query ) / 2 - $expected_score;

    # raw (uncalibrated) BLOSUM distance.
    # d_r = -ln( sigma(N) / sigma(UN) ) * 100
    $BLOSUM_dist = -log( $normalized_score / $normalized_upper_score ) * 100;

    return $BLOSUM_dist;
}

# Title    : compare
# Usage    :
# Function : Calculates overlap of two sequences.
# Args     :
# Returns  : overlap = overlap of two sequences
#            identical = number of identical sequence positions
#
#            If overlap is less than 30 residues it will be set back to 0.
#            This prevents a strong influence of only short overlaps between
#            sequences.
sub compare {

    my ( $taxon, $query, $ref_range, $ref_FASTA ) = @_;

    my $overlap   = 0;
    my $identical = 0;

    for (@$ref_range) {

        next
          if ( substr( $ref_FASTA->{$taxon}, $_, 1 ) =~ /X|\-|\*/
            or substr( $ref_FASTA->{$query}, $_, 1 ) =~ /X|\-|\*/ );

        $overlap++;

        $identical++
          if (
            (
                substr( $ref_FASTA->{$taxon}, $_, 1 ) eq
                substr( $ref_FASTA->{$query}, $_, 1 )
            )
          );

    }
    $overlap = 0 if $overlap < 30;

    return ( $overlap, $identical );
}

# Title    : indices
# Usage    :
# Function :
# Args     :
# Returns  : First and last position in the alignment that is not a gap or X.
sub indices {

    my ( $taxon, $ref_al_length, $ref_FASTA ) = @_;

    $ref_FASTA->{$taxon} =~ /[^X\-]/;

    # offsets of the start of the last successful match
    my ($start) = @-;

    my $rev = reverse $ref_FASTA->{$taxon};
    $rev =~ /[^X\-]/;

    # the ends of the last successful submatches
    my ($end) = @+;

    return ( $start, $$ref_al_length - $end ),;
}

# see also
# https://biology.stackexchange.com/questions/20938/where-to-find-the-corrected-blosum-matrices/20986
# http://web.mit.edu/bamel/blosum/RBLOSUM64
# http://web.mit.edu/bamel/blosum/RBLOSUM62
sub blosumMatrix {

    my ($matrix) = @_;

    my @aminoacids = (
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
        'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'
    );

    #  source: ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62
    #  Matrix made by matblas from blosum62.iij
    #  * column uses minimum score
    #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
    #  Blocks Database = /data/blocks_5.0/blocks.dat
    #  Cluster Percentage: >= 62
    #  Entropy =   0.6979, Expected =  -0.5209
    #  NOTE: This is the original Henikoff BLOSUM62 matrix.
    #  IMPLICATIONS:
    #   1) Matrix is not correct but seems to work.
    #      See doi:10.1038/nbt0308-274
    #   2) J is not included.
    my @blosum62 = (

       #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
        [
            4,  -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1,
            -1, -2, -1, 1,  0, -3, -2, 0, -2, -1, 0,  -4
        ],
        [
            -1, 5,  0,  -2, -3, 1,  0,  -2, 0,  -3, -2, 2,
            -1, -3, -2, -1, -1, -3, -2, -3, -1, 0,  -1, -4
        ],
        [
            -2, 0,  6,  1, -3, 0,  0,  0,  1, -3, -3, 0,
            -2, -3, -2, 1, 0,  -4, -2, -3, 3, 0,  -1, -4
        ],
        [
            -2, -2, 1,  6, -3, 0,  2,  -1, -1, -3, -4, -1,
            -3, -3, -1, 0, -1, -4, -3, -3, 4,  1,  -1, -4
        ],
        [
            0,  -3, -3, -3, 9,  -3, -4, -3, -3, -1, -1, -3,
            -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4
        ],
        [
            -1, 1,  0,  0, -3, 5,  2,  -2, 0, -3, -2, 1,
            0,  -3, -1, 0, -1, -2, -1, -2, 0, 3,  -1, -4
        ],
        [
            -1, 0,  0,  2, -4, 2,  5,  -2, 0, -3, -3, 1,
            -2, -3, -1, 0, -1, -3, -2, -2, 1, 4,  -1, -4
        ],
        [
            0,  -2, 0,  -1, -3, -2, -2, 6,  -2, -4, -4, -2,
            -3, -3, -2, 0,  -2, -2, -3, -3, -1, -2, -1, -4
        ],
        [
            -2, 0,  1,  -1, -3, 0,  0, -2, 8, -3, -3, -1,
            -2, -1, -2, -1, -2, -2, 2, -3, 0, 0,  -1, -4
        ],
        [
            -1, -3, -3, -3, -1, -3, -3, -4, -3, 4,  2,  -3,
            1,  0,  -3, -2, -1, -3, -1, 3,  -3, -3, -1, -4
        ],
        [
            -1, -2, -3, -4, -1, -2, -3, -4, -3, 2,  4,  -2,
            2,  0,  -3, -2, -1, -2, -1, 1,  -4, -3, -1, -4
        ],
        [
            -1, 2,  0,  -1, -3, 1,  1,  -2, -1, -3, -2, 5,
            -1, -3, -1, 0,  -1, -3, -2, -2, 0,  1,  -1, -4
        ],
        [
            -1, -1, -2, -3, -1, 0,  -2, -3, -2, 1,  2,  -1,
            5,  0,  -2, -1, -1, -1, -1, 1,  -3, -1, -1, -4
        ],
        [
            -2, -3, -3, -3, -2, -3, -3, -3, -1, 0,  0,  -3,
            0,  6,  -4, -2, -2, 1,  3,  -1, -3, -3, -1, -4
        ],
        [
            -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1,
            -2, -4, 7,  -1, -1, -4, -3, -2, -2, -1, -2, -4
        ],
        [
            1,  -1, 1,  0, -1, 0,  0,  0,  -1, -2, -2, 0,
            -1, -2, -1, 4, 1,  -3, -2, -2, 0,  0,  0,  -4
        ],
        [
            0,  -1, 0,  -1, -1, -1, -1, -2, -2, -1, -1, -1,
            -1, -2, -1, 1,  5,  -2, -2, 0,  -1, -1, 0,  -4
        ],
        [
            -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3,
            -1, 1,  -4, -3, -2, 11, 2,  -3, -4, -3, -2, -4
        ],
        [
            -2, -2, -2, -3, -2, -1, -2, -3, 2,  -1, -1, -2,
            -1, 3,  -3, -2, -2, 2,  7,  -1, -3, -2, -1, -4
        ],
        [
            0, -3, -3, -3, -1, -2, -2, -3, -3, 3,  1,  -2,
            1, -1, -2, -2, 0,  -3, -1, 4,  -3, -2, -1, -4
        ],
        [
            -2, -1, 3,  4, -3, 0,  1,  -1, 0, -3, -4, 0,
            -3, -3, -2, 0, -1, -4, -3, -3, 4, 1,  -1, -4
        ],
        [
            -1, 0,  0,  1, -3, 3,  4,  -2, 0, -3, -3, 1,
            -1, -3, -1, 0, -1, -3, -2, -2, 1, 4,  -1, -4
        ],
        [
            0,  -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -2, 0,  0,  -2, -1, -1, -1, -1, -1, -4
        ],
        [
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1
        ]

    );

    #  source: ftp://ftp.ncbi.nih.gov/repository/blocks/unix/blosum/blosum.tar.Z
    #  Matrix made by matblas from blosum80.iij
    #  * column uses minimum score
    #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
    #  Blocks Database = /data/blocks_5.0/blocks.dat
    #  Cluster Percentage: >= 80
    #  Entropy =   0.9868, Expected =  -0.7442
    #  NOTE: This is the original Henikoff BLOSUM80 matrix.
    #  IMPLICATIONS:
    #   1) Matrix is not correct but seems to work.
    #      See doi:10.1038/nbt0308-274
    #   2) J is not included.
    my @blosum80 = (

        # A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
        [
            5,  -2, -2, -2, -1, -1, -1, 0, -2, -2, -2, -1,
            -1, -3, -1, 1,  0,  -3, -2, 0, -2, -1, -1, -6
        ],
        [
            -2, 6,  -1, -2, -4, 1,  -1, -3, 0,  -3, -3, 2,
            -2, -4, -2, -1, -1, -4, -3, -3, -2, 0,  -1, -6
        ],
        [
            -2, -1, 6,  1, -3, 0,  -1, -1, 0, -4, -4, 0,
            -3, -4, -3, 0, 0,  -4, -3, -4, 4, 0,  -1, -6
        ],
        [
            -2, -2, 1,  6,  -4, -1, 1,  -2, -2, -4, -5, -1,
            -4, -4, -2, -1, -1, -6, -4, -4, 4,  1,  -2, -6
        ],
        [
            -1, -4, -3, -4, 9,  -4, -5, -4, -4, -2, -2, -4,
            -2, -3, -4, -2, -1, -3, -3, -1, -4, -4, -3, -6
        ],
        [
            -1, 1,  0,  -1, -4, 6,  2,  -2, 1, -3, -3, 1,
            0,  -4, -2, 0,  -1, -3, -2, -3, 0, 3,  -1, -6
        ],
        [
            -1, -1, -1, 1, -5, 2,  6,  -3, 0, -4, -4, 1,
            -2, -4, -2, 0, -1, -4, -3, -3, 1, 4,  -1, -6
        ],
        [
            0,  -3, -1, -2, -4, -2, -3, 6,  -3, -5, -4, -2,
            -4, -4, -3, -1, -2, -4, -4, -4, -1, -3, -2, -6
        ],
        [
            -2, 0,  0,  -2, -4, 1,  0, -3, 8,  -4, -3, -1,
            -2, -2, -3, -1, -2, -3, 2, -4, -1, 0,  -2, -6
        ],
        [
            -2, -3, -4, -4, -2, -3, -4, -5, -4, 5,  1,  -3,
            1,  -1, -4, -3, -1, -3, -2, 3,  -4, -4, -2, -6
        ],
        [
            -2, -3, -4, -5, -2, -3, -4, -4, -3, 1,  4,  -3,
            2,  0,  -3, -3, -2, -2, -2, 1,  -4, -3, -2, -6
        ],
        [
            -1, 2,  0,  -1, -4, 1,  1,  -2, -1, -3, -3, 5,
            -2, -4, -1, -1, -1, -4, -3, -3, -1, 1,  -1, -6
        ],
        [
            -1, -2, -3, -4, -2, 0,  -2, -4, -2, 1,  2,  -2,
            6,  0,  -3, -2, -1, -2, -2, 1,  -3, -2, -1, -6
        ],
        [
            -3, -4, -4, -4, -3, -4, -4, -4, -2, -1, 0,  -4,
            0,  6,  -4, -3, -2, 0,  3,  -1, -4, -4, -2, -6
        ],
        [
            -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1,
            -3, -4, 8,  -1, -2, -5, -4, -3, -2, -2, -2, -6
        ],
        [
            1,  -1, 0,  -1, -2, 0,  0,  -1, -1, -3, -3, -1,
            -2, -3, -1, 5,  1,  -4, -2, -2, 0,  0,  -1, -6
        ],
        [
            0,  -1, 0,  -1, -1, -1, -1, -2, -2, -1, -2, -1,
            -1, -2, -2, 1,  5,  -4, -2, 0,  -1, -1, -1, -6
        ],
        [
            -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4,
            -2, 0,  -5, -4, -4, 11, 2,  -3, -5, -4, -3, -6
        ],
        [
            -2, -3, -3, -4, -3, -2, -3, -4, 2,  -2, -2, -3,
            -2, 3,  -4, -2, -2, 2,  7,  -2, -3, -3, -2, -6
        ],
        [
            0, -3, -4, -4, -1, -3, -3, -4, -4, 3,  1,  -3,
            1, -1, -3, -2, 0,  -3, -2, 4,  -4, -3, -1, -6
        ],
        [
            -2, -2, 4,  4, -4, 0,  1,  -1, -1, -4, -4, -1,
            -3, -4, -2, 0, -1, -5, -3, -4, 4,  0,  -2, -6
        ],
        [
            -1, 0,  0,  1, -4, 3,  4,  -3, 0, -4, -3, 1,
            -2, -4, -2, 0, -1, -4, -3, -3, 0, 4,  -1, -6
        ],
        [
            -1, -1, -1, -2, -3, -1, -1, -2, -2, -2, -2, -1,
            -1, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -6
        ],
        [
            -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,
            -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1
        ]

    );

    my @aminoacids_new = (
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
        'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*'
    );

# Source: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
# Entries for the BLOSUM62 matrix at a scale of ln(2)/2.0.
    my @blosum62_new = (

        #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
        [
            4,  -1, -2, -2, 0,  -1, -1, 0,  -2, -1, -1, -1, -1, -2,
            -1, 1,  0,  -3, -2, 0,  -2, -1, -1, -1, -4
        ],
        [
            -1, 5,  0,  -2, -3, 1,  0,  -2, 0, -3, -2, 2, -1, -3,
            -2, -1, -1, -3, -2, -3, -1, -2, 0, -1, -4
        ],
        [
            -2, 0, 6, 1,  -3, 0,  0, 0,  1, -3, -3, 0, -2, -3,
            -2, 1, 0, -4, -2, -3, 4, -3, 0, -1, -4
        ],
        [
            -2, -2, 1,  6,  -3, 0,  2, -1, -1, -3, -4, -1, -3, -3,
            -1, 0,  -1, -4, -3, -3, 4, -3, 1,  -1, -4
        ],
        [
            0,  -3, -3, -3, 9,  -3, -4, -3, -3, -1, -1, -3, -1, -2,
            -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4
        ],
        [
            -1, 1, 0,  0,  -3, 5,  2, -2, 0, -3, -2, 1, 0, -3,
            -1, 0, -1, -2, -1, -2, 0, -2, 4, -1, -4
        ],
        [
            -1, 0, 0,  2,  -4, 2,  5, -2, 0, -3, -3, 1, -2, -3,
            -1, 0, -1, -3, -2, -2, 1, -3, 4, -1, -4
        ],
        [
            0,  -2, 0,  -1, -3, -2, -2, 6,  -2, -4, -4, -2, -3, -3,
            -2, 0,  -2, -2, -3, -3, -1, -4, -2, -1, -4
        ],
        [
            -2, 0,  1,  -1, -3, 0,  0, -2, 8, -3, -3, -1, -2, -1,
            -2, -1, -2, -2, 2,  -3, 0, -3, 0, -1, -4
        ],
        [
            -1, -3, -3, -3, -1, -3, -3, -4, -3, 4,  2, -3, 1, 0,
            -3, -2, -1, -3, -1, 3,  -3, 3,  -3, -1, -4
        ],
        [
            -1, -2, -3, -4, -1, -2, -3, -4, -3, 2,  4, -2, 2, 0,
            -3, -2, -1, -2, -1, 1,  -4, 3,  -3, -1, -4
        ],
        [
            -1, 2, 0,  -1, -3, 1,  1, -2, -1, -3, -2, 5, -1, -3,
            -1, 0, -1, -3, -2, -2, 0, -3, 1,  -1, -4
        ],
        [
            -1, -1, -2, -3, -1, 0, -2, -3, -2, 1,  2, -1, 5, 0,
            -2, -1, -1, -1, -1, 1, -3, 2,  -1, -1, -4
        ],
        [
            -2, -3, -3, -3, -2, -3, -3, -3, -1, 0,  0, -3, 0, 6,
            -4, -2, -2, 1,  3,  -1, -3, 0,  -3, -1, -4
        ],
        [
            -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,
            7,  -1, -1, -4, -3, -2, -2, -3, -1, -1, -4
        ],
        [
            1,  -1, 1, 0,  -1, 0,  0, 0,  -1, -2, -2, 0, -1, -2,
            -1, 4,  1, -3, -2, -2, 0, -2, 0,  -1, -4
        ],
        [
            0,  -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2,
            -1, 1,  5, -2, -2, 0,  -1, -1, -1, -1, -4
        ],
        [
            -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1,
            -4, -3, -2, 11, 2,  -3, -4, -2, -2, -1, -4
        ],
        [
            -2, -2, -2, -3, -2, -1, -2, -3, 2,  -1, -1, -2, -1, 3,
            -3, -2, -2, 2,  7,  -1, -3, -1, -2, -1, -4
        ],
        [
            0,  -3, -3, -3, -1, -2, -2, -3, -3, 3,  1, -2, 1, -1,
            -2, -2, 0,  -3, -1, 4,  -3, 2,  -2, -1, -4
        ],
        [
            -2, -1, 4,  4,  -3, 0,  1, -1, 0, -3, -4, 0, -3, -3,
            -2, 0,  -1, -4, -3, -3, 4, -3, 0, -1, -4
        ],
        [
            -1, -2, -3, -3, -1, -2, -3, -4, -3, 3,  3, -3, 2, 0,
            -3, -2, -1, -2, -1, 2,  -3, 3,  -3, -1, -4
        ],
        [
            -1, 0, 0,  1,  -3, 4,  4, -2, 0, -3, -3, 1, -1, -3,
            -1, 0, -1, -2, -2, -2, 0, -3, 4, -1, -4
        ],
        [
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4
        ],
        [
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1
        ]

    );

# Source: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM80
# Entries for the BLOSUM80 matrix at a scale of ln(2)/2.0.
    my @blosum80_new = (

        #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
        [
            5,  -2, -2, -2, -1, -1, -1, 0,  -2, -2, -2, -1, -1, -3,
            -1, 1,  0,  -3, -2, 0,  -2, -2, -1, -1, -6
        ],
        [
            -2, 6,  -1, -2, -4, 1,  -1, -3, 0, -3, -3, 2, -2, -4,
            -2, -1, -1, -4, -3, -3, -1, -3, 0, -1, -6
        ],
        [
            -2, -1, 6, 1,  -3, 0,  -1, -1, 0, -4, -4, 0, -3, -4,
            -3, 0,  0, -4, -3, -4, 5,  -4, 0, -1, -6
        ],
        [
            -2, -2, 1,  6,  -4, -1, 1, -2, -2, -4, -5, -1, -4, -4,
            -2, -1, -1, -6, -4, -4, 5, -5, 1,  -1, -6
        ],
        [
            -1, -4, -3, -4, 9,  -4, -5, -4, -4, -2, -2, -4, -2, -3,
            -4, -2, -1, -3, -3, -1, -4, -2, -4, -1, -6
        ],
        [
            -1, 1, 0,  -1, -4, 6,  2, -2, 1, -3, -3, 1, 0, -4,
            -2, 0, -1, -3, -2, -3, 0, -3, 4, -1, -6
        ],
        [
            -1, -1, -1, 1,  -5, 2,  6, -3, 0, -4, -4, 1, -2, -4,
            -2, 0,  -1, -4, -3, -3, 1, -4, 5, -1, -6
        ],
        [
            0,  -3, -1, -2, -4, -2, -3, 6,  -3, -5, -4, -2, -4, -4,
            -3, -1, -2, -4, -4, -4, -1, -5, -3, -1, -6
        ],
        [
            -2, 0,  0,  -2, -4, 1,  0,  -3, 8, -4, -3, -1, -2, -2,
            -3, -1, -2, -3, 2,  -4, -1, -4, 0, -1, -6
        ],
        [
            -2, -3, -4, -4, -2, -3, -4, -5, -4, 5,  1, -3, 1, -1,
            -4, -3, -1, -3, -2, 3,  -4, 3,  -4, -1, -6
        ],
        [
            -2, -3, -4, -5, -2, -3, -4, -4, -3, 1,  4, -3, 2, 0,
            -3, -3, -2, -2, -2, 1,  -4, 3,  -3, -1, -6
        ],
        [
            -1, 2,  0,  -1, -4, 1,  1,  -2, -1, -3, -3, 5, -2, -4,
            -1, -1, -1, -4, -3, -3, -1, -3, 1,  -1, -6
        ],
        [
            -1, -2, -3, -4, -2, 0, -2, -4, -2, 1,  2, -2, 6, 0,
            -3, -2, -1, -2, -2, 1, -3, 2,  -1, -1, -6
        ],
        [
            -3, -4, -4, -4, -3, -4, -4, -4, -2, -1, 0, -4, 0, 6,
            -4, -3, -2, 0,  3,  -1, -4, 0,  -4, -1, -6
        ],
        [
            -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4,
            8,  -1, -2, -5, -4, -3, -2, -4, -2, -1, -6
        ],
        [
            1,  -1, 0, -1, -2, 0,  0, -1, -1, -3, -3, -1, -2, -3,
            -1, 5,  1, -4, -2, -2, 0, -3, 0,  -1, -6
        ],
        [
            0,  -1, 0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2,
            -2, 1,  5, -4, -2, 0,  -1, -1, -1, -1, -6
        ],
        [
            -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2, 0,
            -5, -4, -4, 11, 2,  -3, -5, -3, -3, -1, -6
        ],
        [
            -2, -3, -3, -4, -3, -2, -3, -4, 2,  -2, -2, -3, -2, 3,
            -4, -2, -2, 2,  7,  -2, -3, -2, -3, -1, -6
        ],
        [
            0,  -3, -4, -4, -1, -3, -3, -4, -4, 3,  1, -3, 1, -1,
            -3, -2, 0,  -3, -2, 4,  -4, 2,  -3, -1, -6
        ],
        [
            -2, -1, 5,  5,  -4, 0,  1, -1, -1, -4, -4, -1, -3, -4,
            -2, 0,  -1, -5, -3, -4, 5, -4, 0,  -1, -6
        ],
        [
            -2, -3, -4, -5, -2, -3, -4, -5, -4, 3,  3, -3, 2, 0,
            -4, -3, -1, -3, -2, 2,  -4, 3,  -3, -1, -6
        ],
        [
            -1, 0, 0,  1,  -4, 4,  5, -3, 0, -4, -3, 1, -1, -4,
            -2, 0, -1, -3, -3, -3, 0, -3, 5, -1, -6
        ],
        [
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6
        ],
        [
            -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,
            -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1
        ],

    );

    for ($matrix) {

        if (/^b62$/) {
            print {*STDERR} "[INFO] Using BLOSUM62 matrix.\n";
            return \@blosum62, \@aminoacids;
        }
        elsif (/^b80$/) {
            print {*STDERR} "[INFO] Using BLOSUM80 matrix\n";
            return \@blosum80, \@aminoacids;
        }
        elsif (/^b62n$/) {
            print {*STDERR} "[INFO] Using new BLOSUM62 matrix.\n";
            print {*STDERR}
              "[WARN] New matrix is still an experimental feature.\n";
            return \@blosum62_new, \@aminoacids_new;
        }
        elsif (/^b80n$/) {
            print {*STDERR} "[INFO] Using new BLOSUM80 matrix.\n";
            print {*STDERR}
              "[WARN] New matrix is still an experimental feature.\n";
            return \@blosum80_new, \@aminoacids_new;
        }
        else { die "[ERROR] Unknown matrix." }
    }
}

__END__