
#################################################################################
# This Perl script generates discriminative priors.
# Input arguments:
#     1. size of the motif
#     2. name of a fasta file containing the positive sequences
#     3. name of a fasta file containing the negative sequences
#     4. name of the output file
# Output:
#     a fasta-like prior file that can be used by Priority.
#
# Example run:
#     perl discriminative_from_pos_and_neg.pl 8 TFpos.fasta TFneg.fasta TF.prior
#################################################################################


use strict;
use prioritylib;

if (scalar(@ARGV) != 4) {
   print "Syntax: perl discriminative_from_pos_and_neg.pl n TFpos.fasta TFneg.fasta TF.prior\n";
   print "   n = size of the motif\n";
   print "   TFpos.fasta = a fasta file containing the positive sequences\n";
   print "   TFneg.fasta = a fasta file containing the negative sequences\n";
   print "   TF.prior    = the output file\n";
   exit;
}


my $w = $ARGV[0];
my $positive_fasta_file = $ARGV[1];
my $negative_fasta_file = $ARGV[2];
my $output_file = $ARGV[3];

# generate the counts:
my @array_positive = get_counts_simple($w, $positive_fasta_file);
my @array_negative = get_counts_simple($w, $negative_fasta_file);
my $number_words_pos = sum(@array_positive);
my $number_words_neg = sum(@array_negative);

my $pos_pseudo = 1;  # the number of pseudocounts for the positives
my $all_pseudo = ($number_words_pos + $number_words_neg)/$number_words_pos;

#get the bound probes for this TF
my ($hashref,$keysref) = getfasta_ordered($positive_fasta_file);

print "$output_file\n";
open(OUT, ">$output_file");
foreach my $probe (@$keysref) 
{
   my $seq = $$hashref{$probe};
   print OUT ">$probe\n";
   my $wmer; my $index;
   for (my $pos = 0; $pos<length($seq)-$w+1; $pos++) {
      $wmer = substr($seq, $pos, $w);
      $index = get_index($wmer);
      print OUT ($array_positive[$index] + $pos_pseudo)/
                ($array_positive[$index] + $array_negative[$index] + $all_pseudo), " ";
    }
    for (my $pos = length($seq)-$w+1; $pos<length($seq); $pos++) {
       print OUT 0, " ";
    }
        
    $seq = get_reverse($seq);
    for (my $pos = 0; $pos<length($seq)-$w+1; $pos++) {
       $wmer = substr($seq, $pos, $w);
       $index = get_index($wmer);
       print OUT ($array_positive[$index] + $pos_pseudo)/
                ($array_positive[$index] + $array_negative[$index] + $all_pseudo), " ";
     }
     for (my $pos = length($seq)-$w+1; $pos<length($seq); $pos++) {
        print OUT 0, " ";
     }      
     print OUT "\n";
}
close(OUT);
