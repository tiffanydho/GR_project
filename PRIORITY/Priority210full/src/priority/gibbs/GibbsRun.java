package priority.gibbs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import javax.swing.JOptionPane;

import priority.Parameters;
import priority.Priority;

/**
 * GibbsRun - this class implements the gibbs sampler
 * @author raluca
 */
public class GibbsRun extends Thread
{
	/* Path variables for output files: */
	String fname_output; /* the output directory for a TF (different dir for each TF) */
	String pans;         /* the complete path for the main output file (different file for each TF) */
	String pbest;        /* the best results from pans */
	String plogl;      /* the complete path for the file containing the logls (diff file for each TF) */
	
	/* Running parameters */
	double phi_prior[];      /* dirichlet prior counts for the model parameters (phi) */
	double sampling_exp[][]; /* the exponents for prior and likeliood, used for sampling */
	double back[];           /* the background model */
	String tf_names[];       /* TF names */
	
	int wsize;                  /* the window size */	
	int nc;                     /* the number of classes (prior_dirs.length + 1 for the "other" class) */	
	double comboprior[][][];    /* the priors for each class/prior, each sequence, each position (different for each TF) */
	String comboseq[];          /* the sequences for the current TF (as strings of {0,1,2,3} <-> {A,C,G,T} */
	String comboseq_names[];    /* the names of the genes corresponding to these sequences */
	
	/* Other variables */
	double cprior[];     /* dirichlet prior counts for the classes */
	double denom[][];    /* the normalization constant (for each class and each sequence) */
	double phi[][];      /* the PSSM */
	double phi_temp[][]; /* temporary phi used for printing */
	
	/* The variables we sample */
	int Z[], bestZ[], overallbestZ[] = null;
	int C[], bestC[], overallbestC[] = null;
	
	/* Variable used instead of 0 for the priors */
	double verySmall = Math.exp(-30);
	
	/* Variable used to stop the thread */
	private boolean stop_thread = false;

	/* Pointer to the main application */
	private Priority mainApp;
	
	
	/** Constructor */
	public GibbsRun(Priority mainApp) {
		this.mainApp = mainApp;
	}

	/** Stopps the current thread */
	public void stopThread() {
    	stop_thread = true;
    }
    
	/** Makes local copies of some of the variables in Parameters. */
	public void set_local_params()
	{
		wsize = Parameters.wsize;
		nc = Parameters.prior_dirs.length;
		if (Parameters.otherclass) 
			nc = nc + 1;
		phi_prior = Parameters.phi_prior;
		
		back = Parameters.back;
		tf_names = Parameters.tf_names;	
	}
	
	/** Computes COMBOPRIOR. */
	public void set_priors(int tf) throws GibbsException
	{
		/* comboprior[c][i][j] = the prior for class/priortype c, sequence i, position j */
		comboprior = new double[nc][comboseq.length][];

		/* set all the priors (except the uniform/other) */
		int c, i, j, length;
		for (c=0; c<Parameters.prior_dirs.length; c++) {
			try { 
				String filename = Parameters.prior_dirs[c] + "/" + tf_names[tf] + ".prior";
				File file = new File(filename);
				if ((!file.exists()) || (!file.isFile()) || (!file.canRead())) {
					String mess = "Error: the prior file \"" + filename + 
			    		"\" does not exist or is not readable!";
					if (Priority.useInterface) 
						javax.swing.JOptionPane.showMessageDialog(null, mess, 
								"Gibbs error", JOptionPane.ERROR_MESSAGE);					
					System.out.println(mess+"\n");
				    throw new GibbsException(mess);
			    }
				
				double temp[][] = GibbsStatic.get_positional_prior(filename);
				
				for (i=0; i<comboseq.length; i++)
				{
					if (temp[i].length != comboseq[i].length()) {
						String mess = "The size of a fasta sequence (" + comboseq[i].length() + 
						              ") does not match the size of the prior (" + temp[i].length + ").";
						throw new GibbsException(mess);
					}
				}
				comboprior[c] = temp;
				
				/* if the wsizes (min and max) are specified in the prior file
				 * check that the current wsize is in that range */
				if ((Parameters.wsizeMinFromPrior != -1) &&	(Parameters.wsizeMinFromPrior != -1))
					if ((Parameters.wsize < Parameters.wsizeMinFromPrior) ||
						(Parameters.wsize > Parameters.wsizeMaxFromPrior)) {
						String mess = "Warning: the motif length (" + Parameters.wsize +
						              ") is not in the range specified in the prior file\n  " +
						              filename + " (" + Parameters.wsizeMinFromPrior + ".." +
						              Parameters.wsizeMaxFromPrior + ").";
						if (Priority.useInterface) 
							javax.swing.JOptionPane.showMessageDialog(null, mess, 
									"Warning", JOptionPane.WARNING_MESSAGE);
						System.out.println(mess);					
					}
				
			}
			catch (GibbsException err) {
				throw err;
			}	
		}		
           
		if (Parameters.otherclass) {
			if (nc > 1) /*otherclass is NOT the only class*/
				comboprior[nc-1] = (double[][])(comboprior[0]).clone();
			for (i=0; i<comboseq.length; i++)
			{
				length = comboseq[i].length();
				comboprior[nc-1][i] = new double[length]; 
				for (j=0; j<length; j++)
					comboprior[nc-1][i][j] = Parameters.flat_prior_other_class;
			}
		}
		
		/* scale the prior to be between d and 1-d */
		for (c=0; c<nc; c++)
			for (i=0; i<comboseq.length; i++)
				for (j=0; j<comboseq[i].length(); j++) {
					if (comboprior[c][i][j] == 0)
						comboprior[c][i][j] = this.verySmall; /* e^(-30) */
					if (comboprior[c][i][j] == 1) /* to avoid values of 1 */
						comboprior[c][i][j] = comboprior[c][i][j] - this.verySmall;
					comboprior[c][i][j] = comboprior[c][i][j]
						* (1-2*Parameters.d) + Parameters.d;
				}		
		
		/* set the prior to 0 for all wmers that contain, besides 0,1,2,3 the 
		 * special character '*' (for masked sequences). Do not scale these. */
		boolean foundNotMasked;
		for (i=0; i<comboseq.length; i++) {
			foundNotMasked = false;
			for (j=0; j<comboseq[i].length()-wsize+1; j++)
				if (comboseq[i].substring(j, j+wsize).contains("*")) {
					for (c=0; c<nc; c++)  /* set the priors to 0 */
						comboprior[c][i][j] = 0;
				}
				else foundNotMasked = true;
			if (foundNotMasked == false) {
				String mess = "Error: the DNA sequence \"" + comboseq_names[i].substring(1) 
				+ "\" contains only masked wmers.";
				throw new GibbsException(mess);
			}
		}
		
		/* this is to ensure we do not sample from the middle part or the end part*/
		for (c=0; c<nc; c++)
			for (i=0; i<comboseq.length; i++) {
				if (Parameters.revflag)
					for (j=comboseq[i].length()/2-wsize+1; j<comboseq[i].length()/2; j++)
						comboprior[c][i][j]=0;
				 for (j=comboseq[i].length()-wsize+1; j<comboseq[i].length(); j++)
					comboprior[c][i][j]=0;
			}
	}
	
	
	/** The gibbs sampler */ 
	public void run()
	{
		/* for each TF: */
		for (int tf=0; tf<tf_names.length; tf++)
		{
			/**********/ if (this.stop_thread) break;
			
			/* read the sequences... */
			try { 
				String[] temp = GibbsStatic.get_DNAsequences(tf_names[tf]);
				comboseq = new String[temp.length/2];
				comboseq_names = new String[temp.length/2];
				for (int h=0; h<temp.length; h+=2) {
					comboseq_names[h/2] = temp[h];
					comboseq[h/2] = temp[h+1];
				}
			}
			catch (GibbsException err) {
				if (Priority.useInterface)
					javax.swing.JOptionPane.showMessageDialog(null, err.getMessage(), 
							"Gibbs error", JOptionPane.ERROR_MESSAGE);
				System.out.println(err.getMessage());
				continue;
			}
			
			/**********/ if (this.stop_thread) break;
			
			/* append the reverse complement */
			if (Parameters.revflag)
				for (int i=0; i<comboseq.length; i++)
					comboseq[i] = comboseq[i] + GibbsStatic.get_reverse(comboseq[i]);
			
			/* set the priors for each position in each of the sequences */
			try {
				set_priors(tf); /* modifies comboprior */
			}
			catch (GibbsException err) {
				if (Priority.useInterface)
					javax.swing.JOptionPane.showMessageDialog(null, err.getMessage(), 
							"Gibbs error", JOptionPane.ERROR_MESSAGE);
				System.out.println(err.getMessage());
				continue;
			}

			
			
			/* set the pseudocounts for gamma (the prior on the class) */
			cprior = new double[nc];
			for (int c=0; c<nc; c++) 
				cprior[c] = Parameters.pseudocounts_class;
			if (Parameters.putative_class >= 0) /* the slight advantage when the class is known */
				cprior[Parameters.putative_class] = Parameters.pseudocounts_putative_class; 
			
			
			/**********/ if (this.stop_thread) break;
			
			/* compute the normalization constant */
			denom = new double[nc][comboseq.length];
			int c,i,j;
			for (i=0; i<comboseq.length; i++)
			{
				for (c=0; c<nc; c++) {
					denom[c][i] = 1;
				
					if (comboprior[c][i].length != comboseq[i].length()) {
						String mess = "The size of a fasta sequence (" + 
						              comboseq[i].length() + 
						              ") does not match the size of the prior (" +
						              comboprior[c][i].length + ").";
						if (Priority.useInterface)
							javax.swing.JOptionPane.showMessageDialog(null, mess, 
									"Gibbs error", JOptionPane.ERROR_MESSAGE);
						System.out.println(mess+"\n");
						continue;
					}
				}
				
				for (j=0; j<comboseq[i].length(); j++)
					for (c=0; c<nc; c++) { 
						denom[c][i] = denom[c][i] + comboprior[c][i][j]/(1-comboprior[c][i][j]);
					}
			}
			
			/**********/ if (this.stop_thread) break;
			
			/* create the output files */
			fname_output = Parameters.path_output + "/" + tf_names[tf];	
			pans = fname_output + ".trials.txt";
			pbest = fname_output + ".best.txt";
			plogl = fname_output + ".logl";
			try { initial_print(tf); }
			catch (GibbsException err) {
				if (Priority.useInterface)
					javax.swing.JOptionPane.showMessageDialog(null, err.getMessage(), 
							"Gibbs error", JOptionPane.ERROR_MESSAGE);
				System.out.println(err.getMessage());
				continue;
			}
			
			
			/* initialize sampled variables */
			Z = new int[comboseq.length];
			bestZ = new int[comboseq.length];
			overallbestZ = new int[comboseq.length];
			C = new int[comboseq.length];
			bestC = new int[comboseq.length];
			overallbestC = new int[comboseq.length];
			
			phi = new double[4][wsize];
				
			for (int a=0; a<4; a++)
				for (int b=0; b<wsize; b++)
					phi[a][b] = 0;
			
			double overallbestlogl = -Double.MAX_VALUE;
			int overallbestlogl_trial = -1, overallbestlogl_iter = -1;
			
			int tr, cnt;
			System.out.print("Running TF " + tf_names[tf] + " ");
			for (tr=0; tr<Parameters.trials; tr++) //if (comboseq.length >= 10)
			{
				/**********/ if (this.stop_thread) { break; }
				System.out.print(".");

				Random rand = new Random();
				/* initialize C and Z (uniform probability) */
				for (i=0; i<comboseq.length; i++)
				{
					C[i] = rand.nextInt(nc); /* random integers between 0 and nc-1 */
					
					/* generate a weight vector that has 0 wherever the prior for the
					 * class C[i] has 0, and 1 everywhere else */
					double temp[] = new double[comboseq[i].length()];
					double tempsum = 0;
					for (j=0; j<comboseq[i].length(); j++) {
						temp[j] = comboprior[ C[i] ][i][j];
						if (temp[j] > 0)
							temp[j] = 1;
						tempsum += temp[j];
					}
					Z[i] = GibbsStatic.rand_sample(temp, tempsum);
					if (Z[i] == -1) /* debugging only */
						System.out.println("ATTN: rand_sample returned -1!!!");
				}
					
				if (tr == Parameters.trials) continue;
				
				printToScreen(tf, tr, 0, Z, C, 0, -1, -1);
				
				/**********/ if (this.stop_thread) { break; }
				
				double bestlogl = -Double.MAX_VALUE;
				int bestlogl_iter = -1;				
				double[] logl = new double[Parameters.iter];
				int index = rand.nextInt(comboseq.length); /* choose random index */
				
				
				/* run with different values for power and powerpssm */
				int power = 1;
				int powerpssm = 1;
				boolean flag_for_oldsampling = false;
				
				int n_sampling = Parameters.sampling_exp.length;
				for (int i_sampling=1; i_sampling <= n_sampling; i_sampling++) {
					if (tr < i_sampling * Parameters.trials / n_sampling) {
						powerpssm = (int)Parameters.sampling_exp[i_sampling-1][0];
						power = (int)Parameters.sampling_exp[i_sampling-1][1];
						if ((int)Parameters.sampling_exp[i_sampling-1][2] == 0)
							flag_for_oldsampling = true;
						break;
					}
				}				
				if (powerpssm != 1 || power != 1) 
					Parameters.noocflag = false;
				
				for (cnt=0; cnt<Parameters.iter; cnt++)
				{
					/**********/ if (this.stop_thread) { break; }
					if ((cnt+1) % Parameters.outputStep == 0) 
						printToScreen(tf, tr, cnt, bestZ, bestC, bestlogl, -1, -1);
					
					if((cnt+1)%Parameters.outputStep == 0)
					{
						power=Math.max(1,power-3);
						powerpssm=Math.max(1,powerpssm-1);
						if(power > 1 || powerpssm >1)
							Parameters.noocflag=false;
						else
							Parameters.noocflag=true;
						
					}
					
					/* pick the next sequence (+1 or randomly: index = rand.nextInt(comboseq.length))*/
					index = (index+1) % (comboseq.length);
					
					/* calculate the current estimate of phi */
					
					phi = GibbsStatic.calPhi(Z, comboseq, index, wsize, phi_prior);
					
					/* next we have to sample Z[index], but before that we must compute the weights:*/
					int length = comboseq[index].length();
					double[] W = new double[length+1];
					double sumW = 0;
					
					/* we initialize the weight vector for the class also (will be used later) */
					double[] Wc = new double[nc];
					double sumWc = 0;
					
					/**********/ if (this.stop_thread) { break; }
					
					for (j=0; (j<length) && (stop_thread==false); j++) 
					{	
						if (comboprior[C[index]][index][j] == 0) 
								W[j] = 0;
						else { 
							double quantity1 = GibbsStatic.calA(phi, wsize, 
									comboseq[index].substring(j,j+wsize), back, Parameters.bkgrOrder);
							if(flag_for_oldsampling) {
								quantity1=quantity1*Math.pow(GibbsStatic.calonlypssm(phi, wsize, 
									comboseq[index].substring(j,j+wsize), back, Parameters.bkgrOrder),powerpssm-1);
							}
							else {
								quantity1=Math.pow(quantity1,powerpssm);
							}
							
							double quantity2 = Math.pow(comboprior[C[index]][index][j],power) / (1-Math.pow(comboprior[C[index]][index][j],power));
							W[j] =  quantity1 * quantity2;							
							sumW += W[j];
						}
					}
										
					/**********/ if (this.stop_thread) { break; }
					
					if (GibbsStatic.max(W) <= 0)
						Z[index] = -1; /* debugging case!!! "all weights 0") Should not happen!!! */
					else {
						if (Parameters.noocflag)
							W[length] = 1; /* to sample the "no motif" case */
						else
							W[length] = 0;
						sumW += W[length];
						
						/* now we finally sample Z[index] */
						Z[index] = GibbsStatic.rand_sample(W, sumW);
						
						if(Z[index]==-1)
							System.out.println("ATTN "+index+" "+Z[index]+" "+length);
						
						/* next we compute the weights for sampling C[index] */										
						if (Z[index] == length) /* the sampler picked the appended 1 ("no motif") */ 
						{
							Z[index] = -1; 
							for (c=0; c<nc; c++) 
								Wc[c] = (1/denom[c][index]) * GibbsStatic.prob_gen_prior(c,nc,cprior,C,index);
						} 
						else /* we set the weights to sample class where a motif exists */ 
						{
							for (c=0; c<nc; c++)
							{
								if(Z[index]>length-wsize+1 || Z[index]<0)
								{
									System.out.println(index+" "+Z[index]+" "+length);
								}
								Wc[c] = comboprior[c][index][Z[index]] * 
								        GibbsStatic.prob_gen_prior(c,nc,cprior,C,index);
								
							}
						}
						sumWc = 0;
						for (c=0; c<nc; c++)
							sumWc += Wc[c];
						
					} /* end if (GibbsStatic.max(W) == 0) - the debug case */
					
					/**********/ if (this.stop_thread) { break; }
					
					/* now we sample C[index] */
					C[index] = GibbsStatic.rand_sample(Wc, sumWc);
					
					/* recompute phi */
					phi = GibbsStatic.calPhi(Z, comboseq, -1, wsize, phi_prior);
					
					/**********/ if (this.stop_thread) { break; }
					
					/* now the MAP calculation: */
					double gamma[] = new double[nc];
					double sum_gamma = 0;
					for (int cl=0; cl<nc; cl++) {
						gamma[cl] = GibbsStatic.prob_gen_prior(cl,nc,cprior,C,-1);
						sum_gamma += gamma[cl];
					}
					/* normalize it => a prob distribution */
					for (int cl=0; cl<nc; cl++)
						gamma[cl] = gamma[cl] / sum_gamma;
					
					
					/**********/ if (this.stop_thread) { break; }
					
					/* sum log(P(gamma)) and log(P(phi)) */
					logl[cnt] = GibbsStatic.logl_phi_gamma(phi, phi_prior, gamma, cprior);
					
					/**********/ if (this.stop_thread) { break; }
					
					for (i=0; i<Z.length; i++)
						if (Z[i] > -1) /* there is an occurrence of the motif */
							logl[cnt] = logl[cnt] +
							   Math.log(GibbsStatic.calA(phi, wsize, comboseq[i].substring(Z[i],
									   Z[i]+wsize), back, Parameters.bkgrOrder)) +
							   Math.log(comboprior[C[i]][i][Z[i]]) - 
							   Math.log(1-comboprior[C[i]][i][Z[i]]) - 
							   Math.log(denom[C[i]][i]) + 
							   Math.log(GibbsStatic.prob_gen_prior(C[i], nc, cprior, C, -1));
						else /* no occurrence */
							logl[cnt] = logl[cnt] -
							   Math.log(denom[C[i]][i]) +
							   Math.log(GibbsStatic.prob_gen_prior(C[i], nc, cprior, C, -1));
					
					/**********/ if (this.stop_thread) { break; }
					
					if (logl[cnt] > bestlogl) 
					{
						bestlogl = logl[cnt];
						bestlogl_iter = cnt;
						for (i=0; i<Z.length; i++) bestZ[i] = Z[i];
						for (i=0; i<C.length; i++) bestC[i] = C[i];
					}
				} /* end for (int cnt=0; cnt<iter; cnt++) */
				///**********/ if (this.stop_thread) { break; }
				
				try {
					printCurrentTrialInfo(tr, bestlogl, logl, bestZ, bestC, -1, -1);
				} catch (Exception e) { e.printStackTrace(); }
				if (overallbestlogl < bestlogl) {
					overallbestlogl = bestlogl;
					overallbestlogl_trial = tr;
					overallbestlogl_iter = bestlogl_iter; 
					for (i=0; i<bestZ.length; i++) overallbestZ[i] = bestZ[i];
					for (i=0; i<bestC.length; i++) overallbestC[i] = bestC[i];
				}
			}/* end for (tr=0; tr<this.trials; tr++) */
			System.out.print("\n");
			if (tr > 0) {
				printToScreen(tf, -1, -1, overallbestZ, overallbestC, overallbestlogl, 
						overallbestlogl_trial, overallbestlogl_iter);
				printCurrentTrialInfo(-1, overallbestlogl, null, overallbestZ, overallbestC,
						overallbestlogl_trial, overallbestlogl_iter);
				double[][] overallbest_phi = GibbsStatic.calPhi(overallbestZ, comboseq, -1, wsize, GibbsStatic.noprior);
				Parameters.setOutput(overallbestZ, overallbestC, overallbest_phi, comboseq, comboseq_names, tf_names[tf]);
			}
			
		}/* end for (int tf=0; tf<fname.length; tf++) */
		if (Priority.useInterface) {
			this.mainApp.win.activateStart();
		}
		else {
			System.out.println("Done. The results are available in " + Parameters.path_output);
		}
	   
	}/* end run() */
	
	
	
	/* *********************************************************************************
	 * ************************ Printing functions *************************************
	 * ********************************************************************************/
		
	/** Prints the information about the current trial in the output files. */
	public void printCurrentTrialInfo(int tr, double bestlogl, double logl[], int bestZ[], int bestC[],
			int overallbestlogl_trial, int overallbestlogl_iter)
	{
		/* print info for current trial */
		BufferedWriter file = null;
		try {
			if (tr <0) /* print best overall results */
				file = new BufferedWriter(new FileWriter(pbest, true));
			else {
				file = new BufferedWriter(new FileWriter(pans, true)); /* (tr+1) to start from 1, not 0*/
				file.write("\n\n" + GibbsStatic.repeatChar('-', 30) + " trial " + (tr+1) + " " +
					GibbsStatic.repeatChar('-',30));
			}
			file.write("\nBest score: " + bestlogl + "\n");

			if (tr <0) {/* print best overall results */
				file.write("obtained during trial " + (overallbestlogl_trial+1) + ", iteration " +
						(overallbestlogl_iter+1) + "\n\n");
			}
			
			phi = GibbsStatic.calPhi(bestZ, comboseq, -1, wsize, GibbsStatic.noprior);			
			
			/* print phi */
			int dec = 4;
			file.write("Phi:\n    ");
			for (int j=0; j<wsize; j++) {
				file.write(GibbsStatic.formatInt(j+1, dec+2) + " ");
			}
			file.write("\n");
			for (int i=0; i<4; i++) {
				file.write(" " + GibbsStatic.DNAchars[i] + "  ");
				for (int j=0; j<wsize; j++)
					file.write(GibbsStatic.formatDouble01(phi[i][j], dec) + " ");
				file.write("\n");
			}
			
			/* print the consensus sequence */
			/* *now write the letters with the highest frequencies */
			String letter; 
			int max;
			double phi_temp[][] = new double[phi.length][phi[0].length];
			for (int i=0; i<phi.length; i++)
				for (int j=0; j<phi[0].length; j++)
					phi_temp[i][j] = phi[i][j];
			
			file.write("\n     ");
			String line = "";
			for (int k=0; k<4; k++) {
				line = "";
				for (int j=0; j<phi_temp[0].length; j++) {
					max = 0;
					for (int i=1; i<4; i++)
						if (phi_temp[i][j] > phi_temp[max][j]) max = i;
					letter = "" + GibbsStatic.DNAchars[max];
					if (phi_temp[max][j] < 0.2)
						if (line.trim().length() == 0)
							line = line + "       ";
						else
							file.write("       ");
					else {
						if (line.trim().length() == 0) {
							if (k>0) {
								file.write("     ");
								file.write(line);
							}
							line = ".";
						}
						if (phi_temp[max][j] < 0.5)
							file.write("  " + letter.toLowerCase() + "    ");
						else
							file.write("  " + letter.toUpperCase() + "    ");
					}
					phi_temp[max][j] = -1;
				}
				if (line.trim().length() != 0) file.write("\n");				
			}

			/* print bestZ */
			file.write("\n\nZ: ");
			for (int i=0; i<bestZ.length; i++)  /*(bestZ[i]+1) to make the indices start from 1, not 0, */
				file.write((bestZ[i]+1) + " "); /* and to be consistent with the paper */ 				                                
			file.write("\n");
			
			int i, counts[] = GibbsStatic.class_counts(bestC, bestZ, nc);
			if (Parameters.multiple_priors) {
				/* print bestC and the class/priortype counts  */
				file.write("\nC: ");
				for (i=0; i<bestC.length; i++)
					file.write(this.get_class_number_for_printing(bestC[i]) + " ");
				file.write("\n");					
	
				file.write("\nPrior-type counts: ");
				for (i=0; i<Parameters.prior_dirs.length; i++)
					file.write(" " + this.get_class_number_for_printing(i) + ":" + counts[i] + " ");
				if (Parameters.otherclass)
					file.write(" " + this.get_class_number_for_printing(i) + ":" + counts[i]);
				file.write(" \n");
			}
			file.write("\n");
			
			if (Parameters.noocflag && counts[nc]!=0)
				file.write("Sequences with no occurrence of the motif: " + GibbsStatic.formatInt(counts[nc],3) + "\n\n");			
			
			if (tr < 0 && Parameters.multiple_priors) {/* the final best motif */
				file.write("Predicted prior-type: ");
				int index = 0; /* stores the index of the max */
				for (i=1; i<nc; i++)
					if (counts[index] < counts[i])
						index = i;
				//if ((index == nc-1) && Parameters.otherclass)
				//	file.write("other/uniform prior");
				//else 
				//	file.write(Parameters.prior_dirs[index]);
				file.write(""+get_class_number_for_printing(index));
			}
			file.close();
			
			if (tr >= 0 && Parameters.createLogl == true) {
				file = new BufferedWriter(new FileWriter(plogl, true));
				for (i=0; i<Parameters.iter; i++)
					file.write(logl[i] + " ");
				file.write("\n");					
				file.close();
			}			
		} catch (IOException e) {e.printStackTrace();}
		
	}
	

	/** Prints the parameters of the current run in the output files. */ 
	public void initial_print(int tf) throws GibbsException 
	{
		BufferedWriter file = null, filebest = null;
		try {
			if (Parameters.createLogl) {
				file = new BufferedWriter(new FileWriter(plogl));
				file.close();
			}
		} catch (IOException e) { throw new GibbsException("GibbsException: "+e.getMessage()); }
		try {
			file = new BufferedWriter(new FileWriter(pans));
			filebest = new BufferedWriter(new FileWriter(pbest));
			file.write("\nTranscription factor: " + tf_names[tf] + "\n\n");
			filebest.write("\nTranscription factor: " + tf_names[tf] + "\n\n");
			file.write("Trials: " + Parameters.trials + "\n");
			filebest.write("Trials: " + Parameters.trials + "\n");
			file.write("Iterations/trial: " + Parameters.iter + "\n");
			filebest.write("Iterations/trial: " + Parameters.iter + "\n");
			file.write("Background model order: " + Parameters.bkgrOrder + "\n");
			filebest.write("Background model order: " + Parameters.bkgrOrder + "\n");
			
			file.write("Prior(s):\n");
			filebest.write("Prior(s):\n");
			for (int k=0; k<Parameters.prior_dirs.length; k++) {
				file.write("\t" + this.get_class_number_for_printing(k) + " = " + Parameters.prior_dirs[k] + "\n");
				filebest.write("\t" + this.get_class_number_for_printing(k) + " = " + Parameters.prior_dirs[k] + "\n");
			}
			if (Parameters.otherclass) {
				file.write("\t" + this.get_class_number_for_printing(Parameters.prior_dirs.length) + " = other/uniform prior" + "\n");
				filebest.write("\t" + this.get_class_number_for_printing(Parameters.prior_dirs.length) + " = other/uniform prior" + "\n");
			}
			if (Parameters.putative_class != -1) {
				file.write("Putative prior-type: " + this.get_class_number_for_printing(Parameters.putative_class) + "\n");
				filebest.write("Putative prior-type: " + this.get_class_number_for_printing(Parameters.putative_class) + "\n");
			}
			else {
				file.write("\n");     /* No putative prior-type */
				filebest.write("\n"); /* No putative prior-type */
			}
			file.write("Number of sequences: " + comboseq.length + "\n");
			filebest.write("Number of sequences: " + comboseq.length + "\n");
			file.write("Motif length: " + wsize + "\n");
			filebest.write("Motif length: " + wsize + "\n\n\n");
			file.close();			
			filebest.close();
		} catch (IOException e) { e.printStackTrace(); }
	}

	
	/** Prints to screen some information about the current run (trial).
	 *  When tr = -1, it prints the information about the best run.*/
	public void printToScreen(int tf, int tr, int cnt, int bestZ[], int bestC[], double bestlogl,
			int overallbestlogl_trial, int overallbestlogl_iter)
	{
		if (Priority.useInterface == false)
			return;
		
		String str = GibbsStatic.repeatChar('-', 24);
		double[] noprior = {0,0,0,0};
		if (tr < 0) {/* the final best motif */
			this.mainApp.printOutput(GibbsStatic.repeatChar('-', 21) + 
					" TF: " + tf_names[tf] +", best result " + 
					GibbsStatic.repeatChar('-', 21) + "\n", null);	
			this.mainApp.printOutput("Best score: " + bestlogl + "\n", null);
			this.mainApp.printOutput("obtained during trial " + (overallbestlogl_trial+1) +
					", iteration "  + (overallbestlogl_iter+1) + "\n", null);
			
			this.mainApp.printOutput("PSSM: \n", null);
			phi_temp = GibbsStatic.calPhi(bestZ, comboseq, -1, wsize, noprior);		
			print_phi(phi_temp, null);
			if (Parameters.multiple_priors)
				this.mainApp.printOutput("Prior-type counts: \n", null); 
		}
		else {
			this.mainApp.printOutput(str + " TF: " + tf_names[tf] +", trial " + (tr+1) + " " + str + "\n");	
			this.mainApp.printOutput("Iteration: " + (cnt+1) + "\n");
			this.mainApp.printOutput("Best score: " + bestlogl + "\n");
			this.mainApp.printOutput("PSSM: \n");
			phi_temp = GibbsStatic.calPhi(bestZ, comboseq, -1, wsize, noprior);		
			print_phi(phi_temp);
			if (Parameters.multiple_priors)
				this.mainApp.printOutput("Prior-type counts: \n");
		}
		
		int i, counts[] = GibbsStatic.class_counts(bestC, bestZ, nc);
		
		if (tr < 0) {/* the final best motif */
			if (Parameters.multiple_priors)
			{
				for (i=0; i<Parameters.prior_dirs.length; i++)
					this.mainApp.printOutput(" " + GibbsStatic.formatInt(counts[i],3) + " ", null);
				if (Parameters.otherclass)
					this.mainApp.printOutput(" " + GibbsStatic.formatInt(counts[i],3), null);
		    }	
			this.mainApp.printOutput(" \n",  null);
			if (Parameters.noocflag && counts[nc]!=0)
				this.mainApp.printOutput("Sequences with no occurrence of the motif: " + GibbsStatic.formatInt(counts[nc],3) + "\n", null);
		}
		else {
			if (Parameters.multiple_priors)
			{
				for (i=0; i<Parameters.prior_dirs.length; i++)
					this.mainApp.printOutput(" " + GibbsStatic.formatInt(counts[i],3) + " ");
				if (Parameters.otherclass)
					this.mainApp.printOutput(" " + GibbsStatic.formatInt(counts[i],3));
			}
			this.mainApp.printOutput(" \n");
			if (Parameters.noocflag && counts[nc]!=0)
				this.mainApp.printOutput("Sequences with no occurrence of the motif: " + GibbsStatic.formatInt(counts[nc],3) + "\n");			
		}
		
		if (tr < 0 && Parameters.multiple_priors) {/* the final best motif */
			this.mainApp.printOutput("Predicted prior-type: ", null);
			int index = 0; /* stores the index of the max */
			for (i=1; i<nc; i++)
				if (counts[index] < counts[i])
					index = i;
			//if ((index == nc-1) && Parameters.otherclass)
			//	this.mainApp.printOutput("other (uniform) prior", null);
			//else 
			//	this.mainApp.printOutput(Parameters.prior_dirs[index], null);
			this.mainApp.printOutput(""+index, null);
		}
		this.mainApp.printOutput("\n\n\n", null);
	}
	
	
	/** Prints the PSSM to the screen. */
	public void print_phi(double phi[][])
	{
		if (Priority.useInterface == false)
			return;
		
		int dec = 4;
		this.mainApp.printOutput(GibbsStatic.repeatChar(' ', 4));
		for (int j=0; j<phi[0].length; j++)
			this.mainApp.printOutput(GibbsStatic.formatInt(j+1, dec+2) + " ");
		this.mainApp.printOutput("\n");
		for (int i=0; i<phi.length; i++)
		{
			this.mainApp.printOutput(" " + GibbsStatic.DNAchars[i] + "  ");
			for (int j=0; j<phi[i].length; j++)
				this.mainApp.printOutput(
						GibbsStatic.formatDouble01(phi[i][j], dec) + " ");
			this.mainApp.printOutput("\n");
		}	
	}
	
	/** Prints the PSSM to the screen. */
	public void print_phi(double phi[][], java.awt.Color color)
	{
		if (Priority.useInterface == false)
			return;

		int dec = 4;
		this.mainApp.printOutput(GibbsStatic.repeatChar(' ', 4), color);
		for (int j=0; j<phi[0].length; j++)
			this.mainApp.printOutput(GibbsStatic.formatInt(j+1, dec+2) + " ", color);
		this.mainApp.printOutput("\n", color);
		for (int i=0; i<phi.length; i++)
		{
			this.mainApp.printOutput(" " + GibbsStatic.DNAchars[i] + "  ", color);
			for (int j=0; j<phi[i].length; j++)
				this.mainApp.printOutput(
						GibbsStatic.formatDouble01(phi[i][j], dec) + " ", color);
			this.mainApp.printOutput("\n", color);
		}	
	}
	
	/** Maps the class indices (as they are used in the implementation) to
	 * the class indices used in the paper, i.e. "other" class is always 0
	 * and no other class is 0. */
	private int get_class_number_for_printing(int number) {
		if (Parameters.otherclass) 
			return (number + 1) % nc;
		else
			return number + 1;
	}
} /* end GibbsRun */
