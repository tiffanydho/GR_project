
#################################################################################
# This file generates a Markov model from the sequences in a fasta file.
# Input:
#     $ARGV[0] = the fasta file 
#     $ARGV[1] = the order of the MM
#     $ARGV[2] = the output file
#################################################################################

use strict;

if (scalar(@ARGV) != 3) {
   print "Syntax: perl generateMM.pl file.fasta n file.out\n";
   print "    file.fasta = a fasta file\n";
   print "    n = the order of the Markov model to be generated\n";
   print "    file.out = the output file\n";
   exit;
}

# generate a hash with each possible word as a key.
# give 1 pseudocount for each word.
my %words = ("", 1);
my $order = $ARGV[1];
for (my $i=1; $i<=$order+1; $i++)
{
   my @keys = keys(%words);
   foreach my $key (@keys) 
   {
      $words{$key."A"} = 1;
      $words{$key."C"} = 1;
      $words{$key."G"} = 1;
      $words{$key."T"} = 1;
   }
}
delete $words{""};

# count the number of occurrences for each word
my ($hashref,$keysref) = getfasta_ordered($ARGV[0]);
foreach my $name (@$keysref) 
{
   my $seq = $$hashref{$name};
   # for each size of an word, go through all the words and count
   for (my $wsize=1; $wsize<=$order+1; $wsize++) 
   {
      for (my $i=0; $i<length($seq)-$wsize+1; $i++) {
         $words{uc(substr($seq, $i, $wsize))}++;
      }
   }
}

# generate the MM
open(OUT,">$ARGV[2]");
my @keys = keys(%words);
my @previouslocalkeys = ("");
for (my $i=1; $i<=$order+1; $i++)
{
   # store the keys of size $i
   my @localkeys = ();
   foreach my $key (@keys) {
      if ((length($key) == $i) && ($key =~ /^[ACGT]+$/)) {
         push(@localkeys, $key);
      }
   }

   foreach my $key (@previouslocalkeys) {
      my $sum = $words{$key."A"} + $words{$key."C"} + $words{$key."G"} + $words{$key."T"};
      print OUT $words{$key."A"}/$sum, "\n";
      print OUT $words{$key."C"}/$sum, "\n";
      print OUT $words{$key."G"}/$sum, "\n";
      print OUT $words{$key."T"}/$sum, "\n";
   }

   @previouslocalkeys = sort(@localkeys);
}
close(OUT);



#######################################################################
# Given a fasta input file, return a reference to a hash with the
# name of each sequence as a key for each actual sequence. Thus, if you
# call '($hashref,$keysref) = getfasta($file);' you can see the sequence 
# 'read1' (given that is exists) by typing 'print $$hashref{read1};'.
# $keysref is an ordered list of the actual keys.
####################################################################### 
sub getfasta_ordered
{
  my($filename) = @_;
  chomp($filename);
  my $trying = "Unable to open $filename.\n";
  # Create a new file
  open(INPUTFILE, $filename) or die($trying);
  my @filedata = <INPUTFILE>;
  close INPUTFILE;
  my %seqs;
  my @keys = ();
  my $seqname;
  my $seq = '';
  foreach my $line (@filedata) {
    if($line =~/^\s*$/) { next; }     # Ignore blank lines
    elsif($line =~/^\s*#/) { next; }  # Ignore comments
    elsif($line =~ /^\>\s*(.*)\s*$/) {
      my $temp = $1;
      #$temp =~ s/\s//sg; 
      if(length($seq) > 0) {
        $seq =~ s/\s//sg;
        $seqs{$seqname} = $seq;
        push(@keys, $seqname);
        $seq = '';
      }
      $seqname = $temp;
	if (defined $seqs{$seqname}) {
		print "Error: Duplicate sequence $seqname in file $filename\n"; }
      next;
    } else { $seq .= $line; }
  }
  $seq =~ s/\s//sg; 
  $seqs{$seqname} = $seq;
  push(@keys, $seqname);
  return (\%seqs, \@keys);
}
