     seqfile = ./bact-align.fa             * Path to the alignment file
     treefile = ./unroot-tree.fa           * Path to the tree file
      outfile = ./out-bam.txt              * Path to the output file

        noisy = 3          * How much rubbish on the screen
      verbose = 1          * More or less detailed report

      seqtype = 1          * Data type
        ndata = 1          * Number of data sets or loci
        icode = 0          * Genetic code
    cleandata = 0          * Remove sites with ambiguity data?

        model = 0          * Models for ω varying across lineages
	  NSsites = 0          * Models for ω varying across sites
    CodonFreq = 7          * Codon frequencies
	  estFreq = 1          * Use observed freqs or estimate freqs by ML
        clock = 0          * Clock model
    fix_omega = 0          * Estimate or fix omega
        omega = 1.0         * Initial or fixed omega
