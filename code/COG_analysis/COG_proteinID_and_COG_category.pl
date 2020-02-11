#!/usr/bin/perl
use strict;
use warnings;
######################################################################
#this code will go through the output from the
#get_COG_protein_positions.pl wich has the COG proteinID and its
#start and ending positions. It will then match each of the COG IDs to
#the corresponding borad COG category and print out the start and end
#position as well as the broad COG category 
#to run: COG_proteinID_and_COG_category.pl
#*_COG_proteinIDs_and_positions.txt cognames2003-2014.tab >
#*_COG_category_and_positions.txt 
#output: COG category	start pos	end pos
#######################################################################

my @line = ();
my @COG_protID = ();
my @COG_ID_org = (); #from orginism
my @COG_ID = ();
my @COG_cat = ();
my @start = ();
my @end = ();
my @line2 = ();

#open COG proteinID and pos file
open(FILE, '<:encoding(UTF-8)', $ARGV[0]) || die "can not open file: $!";
while(<FILE>) {
    @line = split /\t/, $_;
    $line[1] =~ s/ // ;
    $line[1] =~ s/\n// ;
    push (@start, $line[1]);
    $line[2] =~ s/ // ;
    $line[2] =~ s/\n// ;
    push (@end, $line[2]);
    push (@COG_protID, $line[0]);
    $line[3] =~ s/\n// ;
    push (@COG_ID_org, $line[3]);
}#while
close(FILE);


#open COG ID names and broad categories file store info
open(FH, '<', $ARGV[1]) || die "can not open file: $!";
while(<FH>) {
    if ($_ =~ /^C/) {
	@line2 = split /\t/, $_;
	$line2[0] =~ s/\n//;
	$line2[0] =~ s/ //;
	push (@COG_ID, $line2[0]);
	$line2[1] =~ s/\n//;
	$line2[1] =~ s/ //;
	push (@COG_cat, $line2[1]);
    }#if
}#while
close(FH);


#search the COG IDs to match with the category 
my @COG_cat_array = ();
for(my $i =0; $i < scalar(@COG_ID_org); $i++) {
    for(my $j =0; $j < scalar(@COG_ID); $j++) {
	if ($COG_ID_org[$i] eq $COG_ID[$j]) {
	    push (@COG_cat_array, $COG_cat[$j]);
	}#if
    }#for
}#for


#print the data into a file
#COG category 	start	end
for(my $n=0; $n < scalar(@start); $n++) {
    print "$COG_cat_array[$n]\t$start[$n]\t$end[$n]\n";
}#for

