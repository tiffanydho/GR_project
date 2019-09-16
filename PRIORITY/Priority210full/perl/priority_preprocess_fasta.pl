
#################################################################################
# This file preprocesses an input fasta file to eliminate any sequence that 
# contains only masked wmers.
# It also writes the fasta file in the format expected by Priority
# Input:
#     $ARGV[0] = the fasta file 
#     $ARGV[1] = w (the word length)
#     $ARGV[2] = the output fasta file
#################################################################################

use strict;

if (scalar(@ARGV) != 3) {
   print "Syntax: perl priority_preprocess_fasta.pl file.fasta w fileout.fasta\n";
   print "    file.fasta = a fasta file\n";
   print "    w = the word length\n";
   print "    fileout.fasta = the output fasta file\n";
   exit;
}


my ($hashref,$keysref) = getfasta_ordered($ARGV[0]);
my $w = $ARGV[1];
open(OUT,">$ARGV[2]") || die("Unable to create output file \"$ARGV[2]\"\n");

my $total = 0;
my $masked = 0;
foreach my $probe (@$keysref)
{
   my $seq = $$hashref{$probe};
   $total++;
   if ($seq =~ /[acgtACGT]{$w}/) {
      print OUT ">$probe\n";
      print OUT $seq,"\n";
   }
   else {
      $masked++;
   }
}
close(OUT);

print "Done. $masked out of $total sequences have been deleted.\n";

#######################################################################
# Given a fasta input file, return a reference to a hash with the
# name of each sequence as a key for each actual sequence. Thus, if you
# call '($hashref,$keysref) = getfasta_ordered($file);' you can see the 
# sequence 'read1' (given that it exists) by typing 'print 
# $$hashref{'read1'};'.
# $keysref is an ordered list of the actual keys.
####################################################################### 
sub getfasta_ordered
{
  my($filename) = @_;
  chomp($filename);
  my $trying = "Unable to open file $filename.\n";
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
