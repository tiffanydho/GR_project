
#####################################################################################
# This Perl script generates discriminative informative priors.
# Input arguments:
#     1. size of the motif
#     2. name of a fasta file containing the positive sequences
#     3. name of a fasta file containing the negative sequences
#     4. name of a fasta-like file containing the information for positive sequences
#     5. name of a fasta-like file containing the information for negative sequences 
#     6. name of the output file
# Output:
#     a fasta-like prior file that can be used by Priority.
#
# Example run:
#     perl discriminative_INFO_from_pos_and_neg.pl 8 TFpos.fasta TFneg.fasta 
#                                                    TFpos.info TFneg.info TF.prior
#####################################################################################

use strict;
use prioritylib;


my $w = $ARGV[0];
my $positive_fasta_file = $ARGV[1];
my $negative_fasta_file = $ARGV[2];
my $positive_fastalike_file = $ARGV[3];
my $negative_fastalike_file = $ARGV[4];
my $output_file = $ARGV[5];

# generate the counts:
my ($ref_counts_pos,$ref_array_pos) = get_counts_info($w, $positive_fasta_file, $positive_fastalike_file);
my ($ref_counts_neg,$ref_array_neg) = get_counts_info($w, $negative_fasta_file, $negative_fastalike_file);

my @array_positive = @$ref_array_pos;
my @array_negative = @$ref_array_neg;

my $avg_info = (sum(@array_positive) + sum(@array_negative)) / (sum(@$ref_counts_pos) + sum(@$ref_counts_neg));

my $pos_pseudo = 1 * $avg_info;  # the number of pseudocounts for the positives
my $all_pseudo = (sum(@$ref_counts_pos) + sum(@$ref_counts_neg)) / sum(@$ref_counts_pos) * $avg_info;


#get the bound probes for this TF
my ($hashref,$keysref) = getfasta_ordered($positive_fasta_file);

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

