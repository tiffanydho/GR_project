      ============================================================
      *                        PRIORITY                          *
      *                     Version 2.1.0                        *
      *                      1 November 2009                         *
      *                                                          *
      *        Priority is licensed from Duke University.        *
      *    Copyright (c) 2006-2009 by Alexander J. Hartemink.    *
      *                   All rights reserved.                   *
      ============================================================
 
     
===== Installation =====

1. Unzip the zip file to the location of your choice.
2. Run priority2.1.0.jar, either by double-clicking it, or by typing 
   "java -jar priority2.1.0.jar" from the command line.  If you are 
   running on Mac OS X, you can simply double click the PRIORITY 
   application bundle instead for a more Mac-ified experience.
   By default PRIORITY Version 2.1.0 will run with a Graphical User 
   Interface (GUI). If you wish to disable the GUI and run the 
   command-line version of our application, please use the "-nogui" 
   option: "java -jar priority2.1.0.jar -nogui".


===== Contents =====

The zip file contains:
README.txt        This file.
ABOUT.txt         Contains a license overview, library license 
                  information, and author information.
license.html      The license under which PRIORITY is released.
priority2.1.0.jar A Java archive file for running PRIORITY.
PRIORITY2.1.0.app An application bundle for running PRIORITY on Mac OS X.
src/              A directory containing Java source code for PRIORITY.
config/           A directory containing configuration files for PRIORITY.
input/            A directory containing example input files for PRIORITY.
output/           A directory for PRIORITY output files.
perl/             A directory containing useful Perl scripts.


===== Example input files =====

1. ./input/TFs.txt
   Contains the names of three transcription factors, for each of which we
   provide FASTA files containing intergenic regions believed to be bound 
   by these TFs. For a complete list of TFs used in our papers (and their 
   corresponding sequence files) please see the PRIORITY web page.
2. ./input/bkgr_3rd_order.txt
   Contains a 3rd order Markov Model generated from yeast intergenic
   regions.
3. ./input/FASTA
   Contains the FASTA files corresponding to the TFs listed in TFs.txt
   Important: the sequences in the fasta file should not be wrapped (i.e.
   each sequence must occupy a single).
4. ./input/C priors
   Contains simple conservation priors for the TFs listed in TFs.txt
5. ./input/DC priors
   Contains discriminative conservation priors for the TFs in TFs.txt


===== Perl scripts =====
1. generateMM_v2.pl 
   Generates an n-th order Markov model from an input FASTA file.
2. discriminative_from_pos_and_neg.pl 
   Generates a discriminative prior from a set of positive and a set of 
   negative sequences.
3. discriminative_from_pos_and_all.pl 
   Generates a discriminative prior from a set of positive sequences and
   a set containing all sequences (positive and negative).
4. generate_fastalike_cons_simple.pl
   Generates the simple conservation prior file.
5. discriminative_INFO_from_pos_and_neg.pl
   Generate a discriminative prior file from simple priors from both the
   positive and negative sequences.

===== Parameters =====

Before running PRIORITY, be sure all the necessary parameters are set to
the correct values and all the input files contain valid data.

The user can specify the configuration file from the command line, using the parameter "-f <name-of-config-file>". By default, the file config/params is used.

You can set the following parameters using PRIORITY's GUI:
   
   1. The transcription factor(s)
         If you want to run PRIORITY just for one TF, then
         check the corresponding radio button and write the TF's 
         name in the text edit. A transcription factor name can be 
         any string, as long as there exists a FASTA sequence file 
         with the same name as the TF, and the extension ".fasta".
         To run PRIORITY for multiple TFs, just write all the names
         of the TFs in a file (one name on each line), check the 
         radio button "Select the file with TF names" and select 
         that file.

   2. The path for the sequence (FASTA) files
         Must be a directory that contains one FASTA file for each 
         selected transcription factor. Each file contains several 
         DNA sequences (in FASTA format) corresponding to DNA regions 
         believed to be bound by same TF.
         The characters accepted in the DNA sequences are "aAcCgGtT".
         All the other characters will be considered 'invalid' and
         the algorithm will not sample any motif that contains invalid
         characters.
         Important!!! The file corresponding to a TF must have the 
         same name as that TF, and the extension ".fasta".

   3. The background model
         The background model file must contain the parameters of 
         a Markov Model for DNA sequences. The model's order (k) may 
         be any integer between 0 and 5.
         For a k-th order model the file must contain exactly
            4 + 4^2 + 4^3 + ... + 4^(k+1)
         real numbers between 0 and 1, one number on each row.
         For example for a 3rd order MM the numbers represent:
            P(A) 
            P(C)
            P(G)
            P(T)
            P(A|A)
            P(C|A)
            P(G|A)
            P(T|A)
            P(A|C)
            ....
            P(T|T)
            P(A|AA)
            P(C|AA)
            .....
            P(T|TT)
            ......
            P(T|TT)
            P(A|AAA)
            P(A|AAC)
            .......
            P(A|TTT)
            ........
            P(T|TTT)
         In this case the file must contain 4+16+64+256=340 numbers.

         Important!!! Notice that every group of 4 consecutive numbers
         must add up to 1, to form a probability distribution.
         Also notice that even if the file contains a k-th order MM,
         you can choose a background model order smaller than k (in this 
         case the program will use only the first 
         4 + 4^2 + 4^3 + ... + 4^k numbers).
         To choose the background model order use the corresponding spinner
         in the bottom left part of the main window.

   4. The path for the output files
         Must be a writable directory. For each TF <tfname> the program
         will generate three output files:
         tfname.trials.txt - contains the running parameters and then for
                             each trial the best score, the corresponding
                             PSSM and the motif positions Z. When multiple
                             priors are used this file also contains C 
                             (the class/prior-type labels, and the 
                             class/prior-type counts (i.e. how many times
                             each class/prior-type appears in C)
         tfname.best.txt   - contains the running parameters, the score,
                             the PSSM and Z. When multiple priors are used
                             this file also contains C and class/prior-
                             type counts for the best scoring motif
                             (corresponding to the best trial)
         tfname.logl       - contains the scores for all the iterations in
                             all the trials. By default this file is not
                             generated (the "logl" flag in the config
                             file is set to "false"). If you wish to 
                             generate this file, please set the "logl" 
                             flag on "true" (see the "Configuration file"
                             section below).

         Important!!! The score is the logarithm of:
            P(data|PSSM,Z,background model) * P(Z|C) * P(C|gamma) * 
            P(phi) * P(gamma) / P(data|Z=0, background)
         where gamma is the prior on the classes and 
         P(data|Z=0, background) is a constant probability. For details
         please refer to the paper.


   5. The priors (prior-types)
         You can use either a single prior-type (uniform or informative) 
         or multiple prior-types (please see the "Priors" menu). If you
         choose to use a single informative prior-type you need to 
         specify the directory that contains the prior files. In this
         directory there must be a file "tfname.prior" for each selected
         TF. The ".prior" files must be in fasta-like format: 
         
              >sequence1 name
              prior[1,1] prior[1,2] prior[1,3] ...
              >sequence2 name
              prior[2,1] prior[2,2] prior[2,3] ...
              ...
              
         where prior[i,j] is the prior probability that the word starting
         at position j in sequence i is a binding site of the selected TF.
         Important!!! The sequences in a ".prior" must be in the same 
         order as in the corresponding fasta file.
         
         If you want to use multiple types of priors (e.g. for several
         classes of TFs) you need to specify for each prior-type a
         directory that contains the prior files for that prior-type
        (e.g. the class priors in ./input).
         Note that if you use more than one type of priors the algorithm
         will sample the prior-type in addition to the locations of the
         binding sites.
         
         When using multiple types of priors you can also set the
         following parameters:
   
         5.1. The "other" prior-type       
           Besides the prior-types specified using the "Add prior-type 
           directory" button, you can also choose to use the "other" 
           prior-type. This is just a uniform prior set by default to 0.5
           (please see the "Configuration file" section below for 
           information on changing this default value).
           To use the "other" prior-type, just make sure the corresponding
           checkbox is selected.
           
        5.2 The putative prior-type
           One of the prior-types selected so far (as can be seen in the
           prior-type list) can be chosen as the putative prior-type. This
           means that when computing the Dirichlet prior for the 
           prior-types the putative prior-type will get more pseudocounts
           than the other prior-types. To find out how to set the exact
           number of pseudocounts please refer to the "Configuration file"
           section below.
           To use a putative prior-type just select the corresponding
           checkbox. The prior-type selected in the list will 
           automatically become the putative prior-type. To change the
           putative prior-type (once the putative prior-type checkbox is
           selected) just double-click on the new prior-type in the 
           prior-type list.

   6. The motif length
        The motif length must a number between 3 and 20. For changing this
        range of values please refer to the "Configuration file" section
        below.
        
        IMPORTANT!!! Please make sure that the prior files are consistent
        with the motif length. For example, if you use a discriminative
        prior that was generated for words of size 8, you should set the
        motif length to 8. The algorithm will work with other motif 
        lengths, but you may get poor results.
        To prevent use of improper motif lengths, we suggest you add two
        lines at the beginning of each prior file, specifying a minimum
        and maximum motif length that should be used with the prior. For
        example, if you use a discriminative prior that was generated for
        size 8, and you want to restrict the use of the prior to a size 8,
        then the first two lines in the prior file should be:                 
           minmotiflength=8
           maxmotiflength=8

   7. The background model order
        Can be any integer between 0 and 5, but it must be consistent with
        the background model file.

  8. The number of trials
        Is the number of random restarts of the Gibbs sampler.

  9. The number of iterations
        Represents the number of iterations per trial.


===== Default values for input parameters =====         

1. The transcription factor(s)
      By default PRIORITY is run for one TF.
      The default file with TF names is "./input/TFs.txt

2. The path for the sequence files
      The default value for this parameter is "./input/FASTA"         

3. The background model
      The default background model file is "./input/bkgr_3rd_order.txt", 
      which contains a third order Markov model.
      The background model order is set to 3 by default.

4. The path for the OUTPUT files
      The default value for this parameter is "./output/"         

5. The priors (prior-types)
      By default, PRIORITY2.1.0 uses a single prior-type: DC 
      (discriminative conservation) generated for motif length 8. Please
      see papers for more details on the priors.

6. The "other" prior-type
      Is not selected by default.

7. The putative prior-type
      Is not selected by default.

8. The motif length
      The default value for this parameter is 8.

9. The background model order
      The default value for this parameter is 3.

10. The number of trials
      The default value for this parameter is 10.

11. The number of iterations
      The default value for this parameter is 10000.

We set the default values for the input parameters such that it is
easy to run the examples we provide in this archive. 
If you want to change these default values please refer to the
"Configuration file" section below.


===== Configuration file =====

The file ./config/params contains the default values for several 
parameters used by PRIORITY:

 1. prior_dirs
    - is the list of directories containing the prior files (one directory
      for each prior-type), separated by ";"

 2. otherpriortype
    - can be true or false

 3. putative_priortype
    - is -1 if no putative prior-type should be used, otherwise the index 
      of the putative prior-type in the prior-type list (the first has 
      index 0)

 4. wsizeMin, wsizeMax, wsize
    - represent the minimum, maximum, and actual motif length

 5. trials
    - the default number of trials

 6. iter
    - the default number of iterations

 7. fname_path
    - the default path for the sequence files

 8. path_output
    - the default path for the output files

 9. back_file
    - the default background file

10. bkgrOrderMin, bkgrOrderMax, bkgrOrder
    - the minimum, maximum and actual background model order

11. individualTF
    - is true if only one TF should be used, false otherwise

12. tf_path
    - the default file with TF names (in case multiple TFs are used)

13. tf_name
    - the default name of the TF (in case only one TF is used)


More advanced parameters:

14. revflag
    - specifies whether or not to look for the motif in the reverse strand
      also (can be true or false)

15. noocflag
    - specifies whether or not to allow the motif not to occur in some of
      the sequences (can be true or false)

16. outputStep
    - after each outputStep iterations, the current results will be
      displayed on the screen

17. pseudocounts_priortype, pseudocounts_putative_priortype
    - the pseudocounts used for computing the Dirichlet prior for the 
      prior-types

18. flat_prior_other_priortype
    - when using the "other" prior-type, its corresponding prior will be
      flat, with this probability everywhere

19. d
    - the scaling parameter for the priors (see paper for more details)

20. phi_prior
    - the pseudocounts for the PSSM

21. sampling_exponents
    - parameters used for sampling. Each sampling setting is described by
      three values: the exponent for the likelihood term, the exponent for
      the prior term, and whether or not to use an exponent for the 
      P(DNA site|background). Different sampling settings are separated 
      by ";".

21. precision
    - the precision used to determine whether or not the background file 
      is valid (every 4 consecutive numbers must add up to 1)

22. logl
    - specified whether or not .logl files should be generated; it is 
      false by default

===== More Information =====

For further information and support, please visit:
http://www.cs.duke.edu/~amink/software/priority/

If you have any questions, issues, or find any bugs, please do not 
hesitate to contact the developers (contact information available on 
the website).

Enjoy!
