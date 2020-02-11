#!/usr/bin/perl
use strict;
use warnings;
######################################################################
#this code will go through the output from the
#COG_proteinID_and_COG_category.pl script. which basically lists the
#COG categry, start and end positions. It will then split up any of
#the lines in that file that have multiple COG categories per protein
#into separate data points. so the output of this script is the COG
#category, start pos, end pos with only one letter in the COG cat
#column.
#############
#to run: COG_proteinID_and_COG_category.pl
#*_COG_category_and_positions.txt >
#*_COG_category_and_positions_multiple_COGs_dealt_with.txt
#############
#output: COG category	start pos	end pos
#######################################################################

my @line = ();
my @COG_cat = ();
my @start = ();
my @end = ();
my @new_COG_cat = ();
my @new_start = ();
my @new_end = ();

#open COG category and position file
open(FILE, '<:encoding(UTF-8)', $ARGV[0]) || die "can not open file: $!";
while(<FILE>) {
    @line = split /\t/, $_;
    $line[1] =~ s/ // ;
    $line[1] =~ s/\n// ;
    push (@start, $line[1]);
    $line[2] =~ s/ // ;
    $line[2] =~ s/\n// ;
    push (@end, $line[2]);
    push (@COG_cat, $line[0]);
}#while
close(FILE);


#make new data pt if there is more than one COG category in each line
#of the original COG cat and position file
for(my $i =0; $i < scalar(@COG_cat); $i++) {
    if($COG_cat[$i] =~ m/[A-Z]{2}/) {
	my @tmp = split //, $COG_cat[$i];
	for(my $j =0; $j < scalar(@tmp); $j++) {
	    push (@new_COG_cat, $tmp[$j]);
	push (@new_start, $start[$i]);
	push (@new_end, $end[$i]);
    }#for
    } else {
	push (@new_COG_cat, $COG_cat[$i]);
	push (@new_start, $start[$i]);
	push (@new_end, $end[$i]);
    }#else
}#for

#print the data into a file
#COG category 	start	end
for(my $n=0; $n < scalar(@new_COG_cat); $n++) {
    print "$new_COG_cat[$n]\t$new_start[$n]\t$new_end[$n]\n";
}#for

