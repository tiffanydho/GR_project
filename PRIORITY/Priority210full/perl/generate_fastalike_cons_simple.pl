
#####################################################################################
# This Perl script generates simple conservation priors.
# Input arguments:
#     1. name of the fasta file
#     2. name of the folder contained fasta files with homologs
#     3. size of the motif
#     4. name of the output (fasta-like) file
# Output:
#     a fasta-like prior file that can be used by Priority.
#
# Example run:
#     perl generate_fastalike_cons_simple.pl TF.fasta homologs/ 8 TF.prior
# For each sequence in TF.fasta, there must be a fasta file in the homologs 
# folder with the same name as the sequence header. For example, given the fasta
# file:
# >seq1
# .....
# there must be a file homologs/seq1.fasta that contains all sequences homologous
# to seq1 (including seq1).
#####################################################################################

use strict;

my ($fastafile, $cons_path, $w, $priorfile) = @ARGV;
print "Running with params: $fastafile, $cons_path, $w, $priorfile\n";

my ($hashref, $keysref) = getfasta_ordered($fastafile);
# go through the sequences once and compute average conservation
my ($sum_cons,$cnt_cons,$max_hom) = (0,0,0);
foreach my $probe (@$keysref)
{  
   if (!(-e $cons_path.$probe)) {
      print "Warning: did not find homologs file $cons_path$probe. Skipped.\n";
      next;
   }
   # open the fasta file with the homologs and read the homologs
   my ($hashref_h, $keysref_h) = getfasta_ordered($cons_path.$probe);
   my @homseqs = ();
   foreach my $key (@$keysref_h)
   {
      push @homseqs, $$hashref_h{$key};
   }
   shift @homseqs; # because the first is the reference organism
      
   # get the sequence in the species of reference, and append its reverse complement;
   my $seq = $$hashref{$probe};
   $seq = $seq . getreverse($seq);
   
   # for each $k-mer in the seq of reference see if it appears in the other species
   my ($wmer, $rwmer, $count, $rcount, $score);
   for (my $pos = 0; $pos<length($seq)-$w+1; $pos++) 
   {      
      $wmer = substr($seq, $pos, $w);
      $rwmer = getreverse($wmer);
       
      $count = 0;    # number of homologs that have the direct wmer
      $rcount = 0;   # number of homologs that do not have the direct wmer, but have the reverse compl
      foreach my $homseq (@homseqs) 
      {
         if (index($homseq, $wmer) >= 0) {
            $count++;
         } elsif (index($homseq, $rwmer) >= 0) { 
            $rcount++;
         }
      } 
      $score = ($count+$rcount)/scalar(@homseqs);       
      $sum_cons += $score; $cnt_cons++;
      if (scalar(@homseqs)>$max_hom) { $max_hom = scalar(@homseqs); } 
   }
}
print "Avg cons: ",$sum_cons/$cnt_cons,"\n";

#average conservation score. Consider this when no homolog is available.
my $avg_cons = $sum_cons/$cnt_cons;

# go through the sequences a second time and compute the prior
open(OUT, ">$priorfile");
foreach my $probe (@$keysref)
{  
   # get the sequence in the species of reference, and append its reverse complement;
   my $seq = $$hashref{$probe};
   $seq = $seq . getreverse($seq);

   print OUT ">$probe\n";
   if (!(-e $cons_path.$probe)) {
      #print "Warning: did not find homologs file $cons_path$name. Skipped.\n";
      for (my $pos = 0; $pos<length($seq)-$w+1; $pos++) 
   	{
   		if (($pos < length($seq)/2-$w+1) || ($pos >= length($seq)/2)) {
       	  print OUT sprintf("%.5f ",$avg_cons);
      	}
      	else {
         	print OUT sprintf("%.5f ",0);
      	}
   	}
   	for (my $pos=length($seq)-$w+1; $pos<length($seq); $pos++) {
      	print OUT sprintf("%.5f ",0);
   	}
   	print OUT "\n";
      next;
   }
   # open the fasta file with the homologs and read the homologs
   my ($hashref_h, $keysref_h) = getfasta_ordered($cons_path.$probe);
   my @homseqs = ();
   foreach my $key (@$keysref_h)
   {
      push @homseqs, $$hashref_h{$key};
   }
    
   shift @homseqs; # because the first is the reference organism
      
   # for each $k-mer in the reference species see if it appears in the other species
   my ($wmer, $rwmer, $count, $rcount, $score);
   for (my $pos = 0; $pos<length($seq)-$w+1; $pos++) 
   {      
      $wmer = substr($seq, $pos, $w);
      $rwmer = getreverse($wmer);
       
      $count = 0;    # number of homologs that have the direct wmer
      $rcount = 0;   # number of homologs that do not have the direct wmer, but have the reverse compl
      foreach my $homseq (@homseqs) 
      {
         if (index($homseq, $wmer) >= 0) {
            $count++;
         } elsif (index($homseq, $rwmer) >= 0) { 
            $rcount++;
         }
      } 

      # the conservation score:
      $score = ($count + $rcount + $avg_cons*($max_hom-scalar(@homseqs)))/$max_hom;        
            
      if (($pos < length($seq)/2-$w+1) || ($pos >= length($seq)/2)) {
         print OUT sprintf("%.5f ",$score);
      }
      else {
         print OUT sprintf("%.5f ",0);
      }
   }
   for (my $pos=length($seq)-$w+1; $pos<length($seq); $pos++) {
      print OUT sprintf("%.5f ",0);
   }
   print OUT "\n";
}
close(OUT);



#######################################################################
# get the reverse complement of a seq
#######################################################################
sub getreverse
{
  my ($string) = @_;
  my $new_string = reverse($string);
  $new_string =~ tr/aAcCgGtT/tTgGcCaA/;
  return $new_string;
}

#######################################################################
# Given a fasta input file, return a reference to a hash with the
# name of each sequence as a key for each actual sequence. Thus, if you
# call '($hashref,$keysref) = getfasta_ordered($file);' you can see the sequence 
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
      #if(length($seq) > 0) {
      if(length($seqname) > 0) {
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

