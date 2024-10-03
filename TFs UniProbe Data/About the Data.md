
## About the dataset
(1) a “raw data” file, containing a summary of the un-normalized probe intensities, flags, and sequences extracted from the GPR files  
(2) an “all data” file, containing all the information above plus the calculated values for
expected Cy3, observed/expected Cy3, Cy3-normalized Alexa488, and spatially-detrended Alexa488  
(3) a “combinatorial” file, containing a ranked list of the normalized intensities and sequences of all combinatorial ‘all k-mer’ probes (with control spots removed)   
(4) a “regression” file, containing the coefficients determined from the R2 linear regression over Cy3 probe intensities and sequences, including the value indicating the quality of the fit.  

## About Seed-and-Wobble
(1) Calculate scores for all individual k-mers  
(2) Constructs PWMs reflecting the relaive preference of each nucleotide

“Seed-and-Wobble” describes our algorithm for PWM construction, in which
the k-mer with the greatest enrichment score (over all possible contiguous and gapped k-mers covered an equal number of times on the array) is chosen as a “seed”, and the
relative preference of each nucleotide variant is calculated by “wobbling” each position
within and outside of the seed. In the process of choosing candidate seeds, a
comprehensive look-up table containing scores for each k-mer is created.