
use strict;
use Cwd;
use List::Util;

if ($#ARGV<3) { die "
################################################################################
# USAGE: script -i kmer_list_or_fasta_file.txt -o output_file.txt 
#          OPTIONAL: -w motif size for PRIORITY (default 8)
#                    -p how many random bases to add on each side (default k*4)
#                    -b bkgr model order for PRIORITY (default 0)
#                    -t number of trials for PRIORITY (default 40)
#                    -a alignment file (default NULL)
#                    -c information content cutoff for refined PWM (default .15)
#                    -s path to PRIORITY software
#                       (default: /net/data/Out/raluca/PRIORITY)
#
# E.g. perl run_priority_for_kmers_v2.pl -i test_kmers.txt -o output.txt
#
# The input file can be a list of kmers (one kmer per line) or a fasta file. 
# When using a fasta file, make sure all sequences are longer than the motif 
# size  (the default for PRIORITY is 8).
#
# The output file will contain the best scoring PWM from PRIORITY (as a prob. 
# matrix), and also refined versions of the matrix that may vary in size.
# The first refined PWM is the PRIORITY PWM + k bases on each side. The second
# refined PWM is a trimmed version of the PRIORITY PWM + k bases on each side.
#
# The padding is only necessary for kmer files (default: k*4). For fasta files 
# the default is 0.
#
# The bkgr model is a 0th order model with equal probabilities. If the user 
# specifies this parameter, the model is generated from input sequences. 
#
# Recommended number of trials: 40 for kmer files, 50-80 for longer sequences.
#
# The alignment file, if set, will contain the alignment of input kmer (seqs).
# For aligned kmers, the kmers are shown in uppercase letters surrounded by ".".
# For aligned sequences, the regions matching the PRIORITY PWM are shown in
# uppercase letters; neighboring bases are shown in lowercase letters.
#
# The information content cutoff is used for trimming the refined PWMs.
################################################################################
";}

my ($kmerfile, $outputfile, $alignmentfile) = ("", "", "", "");
my ($w, $p, $bkgr, $t, $icCutoff) = (8, -1, -1, 40, 0.15);
#my $prioritydir = "/net/data/Out/raluca/PRIORITY";
my $prioritydir = "./";
while (@ARGV) 
{
   $_ = shift;
   if ($_ eq "-i") { $kmerfile = shift; next; } 
   if ($_ eq "-o") { $outputfile = shift; next; } 
   if ($_ eq "-w") { $w = shift; next; } 
   if ($_ eq "-p") { $p = shift; next; } 
   if ($_ eq "-b") { $bkgr = shift; next; } 
   if ($_ eq "-t") { $t = shift; next; } 
   if ($_ eq "-a") { $alignmentfile = shift; next; } 
   if ($_ eq "-c") { $icCutoff = shift; next; } 
   if ($_ eq "-s") { $prioritydir = shift; next; } 
}

################################################################################
# Generate timestamp for unique temporary files & get the current directory
my $time = time();
my $cwd = getcwd();
    
################################################################################
# Read the kmers in the input file
my @kmers;
my $kmer_or_fasta = "kmer"; # By default, it tries to read kmers.
my $mink = 99999;
my $maxk = -99999;
open(IN, $kmerfile)||die("Error: File \"$kmerfile\" does not exist or is not readable.\n");
while (my $line = <IN>)
{
   chomp($line);
   if (length($line)==0) { next; }
   
   if ($line =~ />/) { $kmer_or_fasta = "fasta"; last; }
   
   $line =~ /(\S+)/;
   if (length($1)<$mink) { $mink = length($1); }
   if (length($1)>$maxk) { $maxk = length($1); }
   push @kmers, $1; 
}
close(IN);

if ($kmer_or_fasta eq "fasta") # If the input file is a fasta file...
{
   my ($hashref,$keysref) = getfasta_ordered($kmerfile);
   my $seq;
   $mink = 99999;
   foreach my $probe (@$keysref)
   {
      $seq = uc($$hashref{$probe});
      push @kmers, $seq;
      if (length($seq)<$mink) { $mink = length($seq); }
   }
   print "\nFinished reading fasta file: $kmerfile\n";
}
else { print "\nFinished reading kmer file: $kmerfile\n"; }


################################################################################
# Compute how much padding to add on each side
if ($p < 0) # i.e. not set by user
{
   $p = 4 * $mink;
}
print "Padding added on each side: $p\n";


################################################################################
# Generate the sequences (kmers flanked by random sequences) and write to temp file
open(OUT, ">temp$time.fasta");
my %bases = (0, "A", 1, "C", 2, "G", 3, "T");
my $cnt = 0;
foreach my $kmer (@kmers)
{
   $cnt ++;
   print OUT ">seq$cnt\n";
   for (my $i=0; $i<$p; $i++)
   {
      print OUT $bases{int(rand(4))}; #int(rand(4)) generates an integer from 0 to 3 (inclusive)
   }
   print OUT uc($kmer);
   for (my $i=0; $i<$p; $i++)
   {
      print OUT $bases{int(rand(4))}; 
   }
   print OUT "\n";
}
close(OUT);
print "Fasta file generated: temp$time.fasta\n";


################################################################################
# Generate the background_file from the input sequences (or the default)
if ($bkgr < 0) # i.e. not set by user
{
   $bkgr = 0; # generate equiprobable distribution
   open(OUT, ">tempbkgr$time.txt");
   print OUT "0.25\n0.25\n0.25\n0.25\n";
   close(OUT);
}
else {
   `perl $prioritydir/generateMM.pl temp$time.fasta $bkgr tempbkgr$time.txt`;
}
print "Generated background file (order $bkgr): tempbkgr$time.txt\n";


################################################################################
# Generate the config file for PRIORITY
my $path = cwd(); # get the current path

open(OUT, ">tempconfig$time.txt");

print OUT "fname_path=$cwd\n";
print OUT "path_output=$cwd\n";

print OUT "tf_name=temp$time\n";

print OUT "back_file=$cwd/tempbkgr$time.txt\n";
print OUT "bkgrOrder=$bkgr\n";
print OUT "wsize=$w\n";

print OUT "trials=$t\n";
print OUT "d=0.05\n";

#print OUT "prior_dirs=\n";
print OUT "individualTF=true\n";
print OUT "tf_path=./input/TFs.txt\n\n";

print OUT "wsizeMin=3\n";
print OUT "wsizeMax=20\n";
print OUT "bkgrOrderMin=0\n";
print OUT "bkgrOrderMax=5\n";

print OUT "iter=10000\n";
print OUT "outputStep=1000\n\n";

print OUT "revflag=true\n";
print OUT "noocflag=true\n";
print OUT "phi_prior=0.5 0.5 0.5 0.5\n";
print OUT "sampling_exponents=1 1 1;1 12 1;6 1 1;6 12 1;6 12 0 \n\n";

print OUT "otherpriortype=false\n";
print OUT "putative_priortype=-1\n\n";

print OUT "pseudocounts_priortype=1\n";
print OUT "pseudocounts_putative_priortype=3\n";
print OUT "flat_prior_other_priortype=0.5\n\n";

print OUT "precision=0.0000001\n";
print OUT "logl=false\n";

close(OUT);
print "Generated config file for Priority: tempconfig$time.txt\n";


################################################################################
# Run PRIORITY
chdir($prioritydir); 
print "\nChanged current directory to ",getcwd(),"\n";
print "Running PRIORITY (for $t trials)...\n";
my $cmd = "java -jar priority2.1.0.jar -nogui -f $cwd/tempconfig$time.txt\n";
my $res = `$cmd`;
print "------------------------------------\n$res------------------------------------\n";
chdir($cwd);
print "Changed current directory to ",getcwd(),"\n";

if (length($res) == 0)
{
   unlink("temp$time.fasta");
   unlink("tempbkgr$time.txt");
   unlink("tempconfig$time.txt");
   unlink("temp$time.best.txt");
   unlink("temp$time.trials.txt");
   die("\nError: PRIORITY did not run! You may not have java installed!\n\n");
}


################################################################################
# Parse the Priority output and create the output file
if ($outputfile eq "")
{
   die "Error: output file not specified.\n";
}
open(IN, "temp$time.best.txt");
my @lines = <IN>;
close(IN);

open(OUT, ">$outputfile");
print OUT "PRIORITY motif from $kmerfile, size=$w\n";
my @A = split(/\s+/,$lines[19]); shift(@A); shift(@A); print OUT "A: ",join(" ",@A), "\n";
my @C = split(/\s+/,$lines[20]); shift(@C); shift(@C); print OUT "C: ",join(" ",@C), "\n";
my @G = split(/\s+/,$lines[21]); shift(@G); shift(@G); print OUT "G: ",join(" ",@G), "\n";
my @T = split(/\s+/,$lines[22]); shift(@T); shift(@T); print OUT "T: ",join(" ",@T), "\n";
close(OUT);

print "\nPRIORITY motif (available in $outputfile):\n\n";
print "A: ",join(" ",@A), "\n";
print "C: ",join(" ",@C), "\n";
print "G: ",join(" ",@G), "\n";
print "T: ",join(" ",@T), "\n";
print $lines[24], $lines[25];


################################################################################
# If kmer file then align all kmers to the motif (use the PRIORITY results), 
# write the alignment file and generate a new motif by replacing the random 
# characters with pseudocounts.
# If fasta file, align the sequences and use uppercase letters for the match.
my @Z;
foreach my $line (@lines)
{
   if ($line =~ /^Z:/)
   {
      @Z = split(/\s+/, $line); shift @Z; last;
   }
}

my @cntA; my @cntC; my @cntG; my @cntT; my @cnt; # to store the refined matrix
my $refnewPWM; # to store the refined and trimmed matrix

my ($hashref,$keysref) = getfasta_ordered("temp$time.fasta");
my ($seq, $newseq, $k);
my @newseqs; # will store the aligned region + 10bp each side
my $eachside = 10; #how many characters to print on each side (HARD-CODED, see below)

if ($kmer_or_fasta eq "fasta") ########## fasta
{
   # Note: ignore padding when processing fasta file
   print "\nAligned sequences from fasta file:\n\n";
   if (!($alignmentfile eq "")) { open(OUT,">$alignmentfile"); }
   for (my $i=0; $i<scalar(@Z); $i++)
   {
      $seq = $$hashref{$$keysref[$i]};
      $seq = $seq . get_reverse($seq);
      $newseq = sprintf("%10s",lc( substr($seq,$Z[$i]-1-10,10) )).
                uc( substr($seq,$Z[$i]-1,$w) ).
                lc( substr($seq,$Z[$i]-1+$w,10) );
      push @newseqs, $newseq;
      if ($Z[$i] <= length($seq)/2) # forward, print Z[i]-10 to Z[i]+w+9
      {
         print "$newseq\tF\n";
         if (!($alignmentfile eq "")) { print OUT "$newseq\tF\n"; }
      }
      else  # reverse, print Z[i]-10 to Z[i]+w+9 from forw_rev
      {
         print "$newseq\tR\n";
         if (!($alignmentfile eq "")) { print OUT "$newseq\tR\n"; }
      }
   }
   if (!($alignmentfile eq "")) { close(OUT); }
   print "\n";
   
   # Next, build a new pwm from all aligned sequences
   foreach my $seq (@newseqs)
   {
      my @temp = split(//,uc($seq));
      for (my $i=0; $i<scalar(@temp); $i++)
      {
         if ($temp[$i] eq "A") { $cntA[$i]++; }
         if ($temp[$i] eq "C") { $cntC[$i]++; }
         if ($temp[$i] eq "G") { $cntG[$i]++; }
         if ($temp[$i] eq "T") { $cntT[$i]++; }
         $cnt[$i]++;
      }
   }

   # Build the matrix from counts
   for (my $i=0; $i<scalar(@cntA); $i++)
   {
      $cntA[$i] = sprintf("%.4f",$cntA[$i]/$cnt[$i]);
      $cntC[$i] = sprintf("%.4f",$cntC[$i]/$cnt[$i]);
      $cntG[$i] = sprintf("%.4f",$cntG[$i]/$cnt[$i]);
      $cntT[$i] = 1 - $cntA[$i] - $cntC[$i] - $cntG[$i];
      if (abs($cntT[$i])<0.00001) { $cntT[$i] = 0; }
      $cntT[$i] = sprintf("%.4f",$cntT[$i]);

   }
   
   $" = " ";
      
   my ($refnewPWM) = reduce_matrix_IC([\@cntA,\@cntC,\@cntG,\@cntT], scalar(@cntA), $icCutoff);
   
   print "\nRefined PWM (available in $outputfile):\n\n";
   print "A: @cntA\n";
   print "C: @cntC\n";
   print "G: @cntG\n";
   print "T: @cntT\n";
   
   print "\nRefined PWM trimmed at $icCutoff information content (available in $outputfile):\n\n";
   print "A: @{$$refnewPWM[0]}\n";
   print "C: @{$$refnewPWM[1]}\n";
   print "G: @{$$refnewPWM[2]}\n";
   print "T: @{$$refnewPWM[3]}\n\n";
   
   open(OUT, ">>$outputfile");
   print OUT "\nRefined PWM:\n"; 
   print OUT "A: @cntA\n";
   print OUT "C: @cntC\n";
   print OUT "G: @cntG\n";
   print OUT "T: @cntT\n";
   
   print OUT "\nRefined PWM trimmed at $icCutoff information content cutoff:\n";
   print OUT "A: @{$$refnewPWM[0]}\n";
   print OUT "C: @{$$refnewPWM[1]}\n";
   print OUT "G: @{$$refnewPWM[2]}\n";
   print OUT "T: @{$$refnewPWM[3]}\n";
   close(OUT);
}
else # align kmers w/o using the PRIORITY Z array
{ 
   my @cntA; my @cntC; my @cntG; my @cntT; my @cnt; # to store the refined matrix
   my $refnewPWM; # to store the refined and trimmed matrix
   my ($seq, $newseq, $k, $kmer, $revkmer, $start, $stop, $m, $score);
   my @newseqs; # will store the aligned region + 10bp each side

   print "\nAligned sequences from kmer file :\n\n";
   if (!($alignmentfile eq "")) { open(OUT,">$alignmentfile"); }

   # Read the PRIORITY matrix and add $maxk 0.25's on each side 
   my %pwmPRIORITY;
   read_motif_from_PBM_file_into_hash(\%pwmPRIORITY, "$outputfile");
   my @temp;
   for (my $i=0; $i<$maxk; $i++) { push @temp, 0.25; }
   $pwmPRIORITY{"A"} = [@temp, @{$pwmPRIORITY{"A"}}, @temp];
   $pwmPRIORITY{"C"} = [@temp, @{$pwmPRIORITY{"C"}}, @temp];
   $pwmPRIORITY{"G"} = [@temp, @{$pwmPRIORITY{"G"}}, @temp];
   $pwmPRIORITY{"T"} = [@temp, @{$pwmPRIORITY{"T"}}, @temp];
   
   # Align each kmer to the (extended) pwm
   for (my $i=0; $i<scalar(@kmers); $i++)
   {
      $kmer = uc($kmers[$i]);
      $revkmer = get_reverse($kmer);
      $k = length($kmer);
      $m = $k/2; if ($w < $k) { $m = $w/2; }
      $start = $maxk - ($k-$m);
      $stop = $maxk + $w - $m;
      
      # get the best alignment
      my ($maxscore, $maxpos, $or) = (-99999, -1, "");
      for (my $pos = $start; $pos <= $stop; $pos++)
      {
         $score = 1;
         for (my $j=0; $j<$k; $j++)
         {
            $score = $score * $pwmPRIORITY{substr($kmer, $j, 1)}[$pos+$j];
         }
         if ($score > $maxscore) { ($maxscore, $maxpos, $or) = ($score, $pos, "F"); }
         $score = 1;
         for (my $j=0; $j<$k; $j++)
         {
            $score = $score * $pwmPRIORITY{substr($revkmer, $j, 1)}[$pos+$j];
         }
         if ($score > $maxscore) { ($maxscore, $maxpos, $or) = ($score, $pos, "R"); }
      }

      $newseq = ("." x $maxpos);
      if ($or eq "F") { $newseq .= $kmer; }
      else { $newseq .= $revkmer; }
      $newseq .= ("." x (2*$maxk+$w-$maxpos-$k));
      
      push @newseqs, $newseq;
      
      if ($or eq "F") # forward
      {
         print "$newseq\tF\n";
         if (!($alignmentfile eq "")) { print OUT "$newseq\tF\n"; }
      }
      else  # reverse
      {
         print "$newseq\tR\n";
         if (!($alignmentfile eq "")) { print OUT "$newseq\tR\n"; }
      }
   }
   if (!($alignmentfile eq "")) { close(OUT); }
   
   # Next, build a new pwm from all aligned sequences
   foreach my $seq (@newseqs)
   {
      my @temp = split(//,uc($seq));
      for (my $i=0; $i<scalar(@temp); $i++)
      {
         if ($temp[$i] eq "A") { $cntA[$i]++; $cnt[$i]++; }
         if ($temp[$i] eq "C") { $cntC[$i]++; $cnt[$i]++; }
         if ($temp[$i] eq "G") { $cntG[$i]++; $cnt[$i]++; }
         if ($temp[$i] eq "T") { $cntT[$i]++; $cnt[$i]++; }
         if ($temp[$i] eq ".") { $cntA[$i]+=0.05;  $cntC[$i]+=0.05; $cntG[$i]+=0.05; $cntT[$i]+=0.05; $cnt[$i]+=0.2;}
      }
   }
      
   # Build the matrix from counts
   for (my $i=0; $i<scalar(@cntA); $i++)
   {
      $cntA[$i] = sprintf("%.4f",$cntA[$i]/$cnt[$i]);
      $cntC[$i] = sprintf("%.4f",$cntC[$i]/$cnt[$i]);
      $cntG[$i] = sprintf("%.4f",$cntG[$i]/$cnt[$i]);
      $cntT[$i] = 1 - $cntA[$i] - $cntC[$i] - $cntG[$i];
      if (abs($cntT[$i])<0.00001) { $cntT[$i] = 0; }
      $cntT[$i] = sprintf("%.4f",$cntT[$i]);
   }
   
   $" = " ";
      
   my ($refnewPWM) = reduce_matrix_IC([\@cntA,\@cntC,\@cntG,\@cntT], scalar(@cntA), $icCutoff);
   
   print "\nRefined PWM (available in $outputfile):\n\n";
   print "A: @cntA\n";
   print "C: @cntC\n";
   print "G: @cntG\n";
   print "T: @cntT\n";
   
   print "\nRefined PWM trimmed at $icCutoff information content (available in $outputfile):\n\n";
   print "A: @{$$refnewPWM[0]}\n";
   print "C: @{$$refnewPWM[1]}\n";
   print "G: @{$$refnewPWM[2]}\n";
   print "T: @{$$refnewPWM[3]}\n\n";
   
   open(OUT, ">>$outputfile");
   print OUT "\nRefined PWM:\n"; 
   print OUT "A: @cntA\n";
   print OUT "C: @cntC\n";
   print OUT "G: @cntG\n";
   print OUT "T: @cntT\n";
   
   print OUT "\nRefined PWM trimmed at $icCutoff information content cutoff:\n";
   print OUT "A: @{$$refnewPWM[0]}\n";
   print OUT "C: @{$$refnewPWM[1]}\n";
   print OUT "G: @{$$refnewPWM[2]}\n";
   print OUT "T: @{$$refnewPWM[3]}\n";
   close(OUT);  
}

unlink("temp$time.fasta");
unlink("tempbkgr$time.txt");
unlink("tempconfig$time.txt");
unlink("temp$time.best.txt");
unlink("temp$time.trials.txt");


#######################################################################
# Trim a PWM by information content
#######################################################################
sub reduce_matrix_IC
{
	my ($ref, $W, $ic_cutoff) = @_;
	my @matrix = @$ref;
	# compute information content for each position
	my @ic = (); my $ln2 = log(2);
	for (my $i=0; $i<$W; $i++)
	{
		my $score = 0;
		for (my $j=0; $j<4; $j++) { 
			$score -= $matrix[$j][$i] * mylog($matrix[$j][$i])/$ln2;
			#print $matrix[$j][$i]," ";
		}
		push @ic, sprintf("%.4f",2-$score);
		#print "\n";
	}
   # print "IC: @ic\n";
	
	# retain the positions that have ic >= input
	my ($left, $right) = (0,$W-1);
	while ($ic[$left] < $ic_cutoff) {$left++;}
    while ($ic[$right] < $ic_cutoff) {$right--;}
   
	my @new_matrix = ();
	for (my $i=0; $i<4; $i++) {
		my @local = ();
		for (my $j=$left; $j<=$right; $j++) {
		   push @local, $matrix[$i][$j];
		}
		push @new_matrix, \@local;
	}
	return (\@new_matrix, $right-$left+1);
}


#######################################################################
# Slightly modified version of the log function
#######################################################################
sub mylog
{
	if ($_[0]==0) { return 0; }
	else { return log($_[0]); }
}

#######################################################################
# Get the reverse complement of a sequence
#######################################################################
sub get_reverse
{
  my ($string) = @_;
  my $new_string = join("",reverse(split(//,$string)));
  $new_string =~ tr/aAcCgGtT/tTgGcCaA/;
  return $new_string;
}

#######################################################################
# Read a motif from a PBM-derived PWM file (with frequencies only)
# Save the motif as a hash.
#######################################################################
sub read_motif_from_PBM_file_into_hash()
{
   my ($refhash, $filename) = @_;
   open(IN, $filename) or die ("Error: file $filename\n");
   my @letter = ("A","C","G","T");
   my $dummy = <IN>;
   for (my $i=0; $i<4; $i++) 
   {
     my $line = <IN>; chomp($line);
     my @array = split(/\s+/, $line); shift @array;
     $refhash->{$letter[$i]} = \@array;
   }
   close(IN);
}

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
