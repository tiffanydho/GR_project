#######################################################################
# This Perl file contains functions used by the Perl scripts
# ... and ....
#######################################################################


$| = 1;


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
      $temp =~ s/\s//sg; 
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



#######################################################################
# Given a fasta-like input file, return a reference to a hash with the
# name of each sequence as a key for each actual sequence.
####################################################################### 
sub getfastalike_ordered
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
      $temp =~ s/\s//sg; 
      if(length($seq) > 0) {
        # $seq =~ s/\s//sg;
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
  $seqs{$seqname} = $seq;
  push(@keys, $seqname);
  return (\%seqs, \@keys);
}


#######################################################################
# Sum the elements of an array
#######################################################################
sub sum{
  my $sum = 0;
  foreach (@_) { $sum += $_; }
  return $sum;
}

#######################################################################
# Average of the elements in the array
#######################################################################
sub average{
   my $sum = 0;
   foreach (@_) { $sum += $_; }
   return $sum/scalar(@_);
}

#######################################################################
# Minimum of the elements in the array
#######################################################################
sub min{
   my $min = $_[0];
   foreach (@_) { 
      if ($min > $_) { $min = $_; }
   }
   return $min;
}

#######################################################################
# Maximum of the elements in the array
#######################################################################
sub max{
   my $max = $_[0];
   foreach (@_) { 
      if ($max < $_) { $max = $_; }
   }
   return $max;
}

#######################################################################
# Get the reverse complement of a sequence
#######################################################################
sub get_reverse
{
  my ($string) = @_;
  my $new_string = reverse($string);
  $new_string =~ tr/aAcCgGtT/tTgGcCaA/;
  return $new_string;
}


#######################################################################
# compute the index in the count table (00000000 = 0 = AAAAAAAA, 
# 0..01 = 1 = AAAAAAAC ....)
#######################################################################
sub get_index
{
   my ($wmer) = @_;
   $wmer =~ tr/aAcCgGtT/00112233/;
   my $index = 0;
   
   for (my $i=0; $i<length($wmer); $i++) {
      $index = $index * 4 + substr($wmer,$i,1);
   }
	return $index;
}


#######################################################################
# Counts the words of a particular length from a fasta file and writes
# the counts into the output array. Also takes into account reverse 
# complements.
# Input:
#    1. size of the word
#    2. fasta file to be read
#######################################################################
sub get_counts_simple
{
   my ($w, $filein) = @_;
   chomp($filein);
   
   my ($hashref,$keysref) = getfasta_ordered($filein);
   my %hash = ();
   
   foreach (@$keysref)
   {
      my $probeid = $_;
      my $probe = $$hashref{$probeid};
      
      for (my $i=0; $i<=length($probe)-$w; $i++) 
      {
         my $seq = substr($probe,$i,$w);
         if ($seq !~ /^[acgtACGT]+$/) {
            next;
         }
         
         if (exists($hash{substr($probe,$i,$w)})) {
            $hash{substr($probe,$i,$w)}++;
         }
         else {
            $hash{substr($probe,$i,$w)} = 1;
         }
         my $rev = get_reverse(substr($probe,$i,$w));
         if ($rev ne substr($probe,$i,$w))
         {
            if (exists($hash{$rev})) {
               $hash{$rev}++;
            } 
            else {
               $hash{$rev} = 1;
            }
         }
      }
   }
   
   my @sorted = sort(keys(%hash)) ;
   
   # write the results in the output array
   my @counts = ();
   my $cnt = 0;
   my $nextcnt;
   for my $j (0 .. $#sorted)
   {
      $nextcnt = 0;
      my @sep = split(//,$sorted[$j]);
      for my $k (0 .. $w-1)
      {
         my $expo = $w - 1 - $k;
         if ($sep[$k] eq 'A'|| $sep[$k] eq 'a') {
            $nextcnt = $nextcnt + (4**$expo) * 0;
         }
         elsif ($sep[$k] eq 'C'||$sep[$k] eq 'c') {
            $nextcnt = $nextcnt + (4**$expo) * 1;
         }
         elsif ($sep[$k] eq 'G'||$sep[$k] eq 'g') {
            $nextcnt = $nextcnt + (4**$expo) * 2;
         }
         else {
            $nextcnt = $nextcnt + (4**$expo) * 3;
         }
      }
      while ($nextcnt > $cnt) {
         push(@counts,0);
         $cnt++;
      }
      push(@counts, $hash{$sorted[$j]} );
      $cnt++;
   }
   while ($cnt < 4**$w) {
      push(@counts,0);
      $cnt++;
   }

   return @counts;
}


#######################################################################
# Counts the words of a particular length from a fasta file and writes
# the counts into the output array. Also takes into account reverse 
# complements.
# Input:
#    1. size of the word
#    2. fasta file to be read
#    3. probabilities in an analogous fasta format
#######################################################################

sub get_counts_info
{
   my ($w, $filein, $filein_prob) = @_;
   
   my ($hashref,$keysref) = getfasta_ordered($filein);
   my ($hashref_prob,$keysref_prob) = getfastalike_ordered($filein_prob);
   
   my %hash = ();
   my %hashprob = ();
   
   foreach (@$keysref)
   {
      my $probeid = $_;
      my $probe = $$hashref{$probeid};
      my $probe_prob = $$hashref_prob{$probeid};

      if ($probe_prob) {
         #print $probeid;
      }
      else {
         print "Error: probe $probeid not in file $filein_prob.\n";
      }
      my @probarray = split(/\s/, $probe_prob);
      
      for (my $i=0; $i<=length($probe)-$w; $i++) 
      {
         my $seq = substr($probe,$i,$w); #the seq
         if ($seq !~ /^[acgtACGT]+$/) {
            next;
         }
         my $value = $probarray[$i]; # the data

         if (exists($hash{substr($probe,$i,$w)})) {
            $hash{substr($probe,$i,$w)}++;
            $hashprob{substr($probe,$i,$w)} = $hashprob{substr($probe,$i,$w)} + $value;
         }
         else {
            $hash{substr($probe,$i,$w)} = 1;
            $hashprob{substr($probe,$i,$w)} = $value;
         }
         my $rev = get_reverse(substr($probe,$i,$w));
         if ($rev ne substr($probe,$i,$w))
         {
            if (exists($hash{$rev})) {
               $hash{$rev}++;
               $hashprob{$rev} = $hashprob{$rev} + $value;
            } 
            else {
               $hash{$rev} = 1;
               $hashprob{$rev} = $value;
            }
         }
      }
   }
   
   my @sorted = sort(keys(%hash)) ;
   
   # write the results in the output array
   my @counts = ();
   my @counts_occurrences = ();
   my $cnt = 0;
   my $nextcnt;
   for my $j (0 .. $#sorted)
   {
      $nextcnt = 0;
      my @sep = split(//,$sorted[$j]);
      for my $k (0 .. $w-1)
      {
         my $expo = $w - 1 - $k;
         if ($sep[$k] eq 'A'|| $sep[$k] eq 'a') {
            $nextcnt = $nextcnt + (4**$expo) * 0;
         }
         elsif ($sep[$k] eq 'C'||$sep[$k] eq 'c') {
            $nextcnt = $nextcnt + (4**$expo) * 1;
         }
         elsif ($sep[$k] eq 'G'||$sep[$k] eq 'g') {
            $nextcnt = $nextcnt + (4**$expo) * 2;
         }
         else {
            $nextcnt = $nextcnt + (4**$expo) * 3;
         }
      }
      while ($nextcnt > $cnt) {
         push(@counts,0);
         push(@counts_occurrences,0);
         $cnt++;
      }
      push(@counts, $hashprob{$sorted[$j]} );
      push(@counts_occurrences, $hash{$sorted[$j]} );
      $cnt++;
   }
   while ($cnt < 4**$w) {
      push(@counts,0);
      push(@counts_occurrences,0);
      $cnt++;
   }

   return (\@counts_occurrences, \@counts);
}
