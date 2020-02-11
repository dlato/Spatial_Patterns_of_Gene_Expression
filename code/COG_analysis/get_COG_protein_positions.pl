#!/usr/bin/perl
use strict;
use warnings;
######################################################################
#this code will go through the .csv file from the COG database which
#contains the genome name and COG protein ID. It will also go through
#the .tab file which contains the COG protein ID and accociated
#acession numbers. It will then go through the gbk file and find these
#protein IDs and associated genomic positions
#to run: get_COG_protein_positions.pl *COG_proteinIDs.csv
#COG_protein_and_refIDs.tab *gbk_prot_and_pos.txt 
#
#output: NCBI protein ID	Start Pos	End Pos		COG ID
#######################################################################

my @line = ();
my @COG_posID = ();
my @COG_ID = ();
my @line2 = ();
my @COG_posID_to_search = ();
my @refID = ();
my @gbkID = ();

#open COG protein IDs file
open(FILE, '<:encoding(UTF-8)', $ARGV[0]) || die "can not open file: $!";
while(<FILE>) {
    @line = split /,/, $_;
    $line[2] =~ s/\n//;
    $line[2] =~ s/ //;
    push (@COG_posID, $line[2]);
    #below will grab the COG ID which is needed later on to determine
    #the broad COG category
    push (@COG_ID, $line[6]);
}#while
close(FILE);
#print "@COG_posID" . "\n";
#print "@COG_ID" . "\n";

#open COG protein IDs and refIDs file
open(FILE, '<:encoding(UTF-8)', $ARGV[1]) || die "can not open file: $!";
while(<FILE>) {
    @line2 = split /\t/, $_;
    $line2[0] =~ s/\n//;
    $line2[0] =~ s/ //;
    push (@COG_posID_to_search, $line2[0]);
    $line2[1] =~ s/\n//;
    $line2[1] =~ s/ //;
    push (@refID, $line2[1]);
    $line2[2] =~ s/\n//;
    $line2[2] =~ s/ //;
    push (@gbkID, $line2[2]);
}#while
close(FILE);
#print "@COG_posID_to_search" . "\n";
#print "@refID" . "\n";
#print "@gbkID" . "\n";

#search the COG protein IDs to match with refIDs
my @protein_array_pos = ();
for(my $i =0; $i < scalar(@COG_posID); $i++) {
    for(my $j =0; $j < scalar(@COG_posID_to_search); $j++) {
	if ($COG_posID[$i] == $COG_posID_to_search[$j]) {
	    push (@protein_array_pos, $j);
	}#if
    }#for
}#for
#
#match the COG protein ID with corresponding ref and genome IDs
my @protein_refID =();
while (my $element=shift(@protein_array_pos)) {
    push (@protein_refID, $refID[$element]);
}#while
#print "protein_refID: " . "@protein_refID" . "\n";

#searching the prot and pos files
my @gbk_pos = ();
my @gbk_prot = ();
my @gbk_start = ();
my @gbk_end = ();
open(FILE, '<:encoding(UTF-8)', $ARGV[2]) || die "can not open file: $!";
while(<FILE>) {
    if ($_ =~ "CDS") {
	my $pos = $_;
	$pos =~ m/(\d+)..(\d+)/;
	my $start = $1;
	my $end = $2;
	my $sum = $start + $end;
#	my $half = $sum/2;
#	$half = int($half);
#	push(@gbk_pos, $half);
	push(@gbk_start, $start);
	push(@gbk_end, $end);
    }#if
    if ($_ =~ "protein_id") {
	my $prot = $_;
	$prot =~ m/"/;
	$prot =~ s/$`//;
	$prot =~ m/\./;
	$prot =~ s/$'//;
	$prot =~ s/"//;
	$prot =~ s/\.//;
	push(@gbk_prot, $prot);
    }#if
}#while
close(FILE);

#print "gbk_start: " . "@gbk_start" . "\n";
#print "gbk_end: " . "@gbk_end" . "\n";
#print "gbk_prot: " . "@gbk_prot" . "\n";


#make coordinates for proteins
my @protein = ();
my @position_start = ();
my @position_end = ();
for(my $i =0; $i < scalar(@gbk_prot); $i++) {
    for(my $j =0; $j < scalar(@protein_refID); $j++) {
	if ($gbk_prot[$i] eq $protein_refID[$j]) {
	    push (@protein, $protein_refID[$j]);
	    push (@position_start, $gbk_start[$i]);
	    push (@position_end, $gbk_end[$i]);
	}#if
    }#for
}#for

#print the data into a file
#protein 	start	end	COG ID
for(my $n=0; $n < scalar(@protein); $n++) {
    print "$protein[$n]\t$position_start[$n]\t$position_end[$n]\t$COG_ID[$n]\n";
}#for

