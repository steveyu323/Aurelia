## 210507

### Group meeting summary

1. **TO-DO:** Wrap up the differential gene expression analysis (especially the overlap with sleuth/DESeq2) the findings may also help select qPCR probes
2. The word "blastema" seems abusively used in the world of regenerative science, should be much more careful of its definition when coming into reading the papers out there. Meanwhile, wnt/ $\beta$-catenin pathway probably incorporate more in the re-patterning of the cells than the initial stage of the cell-proliferation/de-differentiation.
3. **TO-DO:** Searching for more conservative regenerative markers across species, and/or cell proliferation markers. Then, map them back onto the Aurelia genome
4. **TO-DO:** Email Lea to schedule sometime next week for a more detailed discussion on the experimental setup of the initial experiments, some current thoughts:
   1. Day 3 is the current first morphological change taking place time point, but what is happening from Day1 to Day3
   2. Order L-Leucine and Insulin for use:
      1. L-Leucine methyl ester hydrochloride [https://www.sigmaaldrich.com/catalog/product/aldrich/l1002?lang=en&region=US]
      2. Insulin human [https://www.sigmaaldrich.com/catalog/product/sigma/i0908?lang=en&region=US]
   3. Real-Time bright field capturing of the first 2 days? Whether it is plausible/worthwhile to perform.



###  Differential gene expression analysis refinement

#### Design and Docs

1. Since the focus is early onset of gene program, a comparison between control VS feeding group at 27 hour should be done. However, it should also be noted that based on previous PCA mapping via both DESeq2 and raw count method (in the EdgeR) analysis, one of the 27-feeding group much closer to the control group and the other one grouped much closer to the two 51 hour samples. Based on these facts, it would be worthwhile to perform the following comparisons:
   1. Pre-amputation VS Amputation
   2. 27h feeding VS 27h control
   3. 3 putatively regenerating amputated samples  VS 5 putatively non-regenerating amputated samples
   4. One putatively regenerating amputated sample at 27h VS two 27h control
   5. Two 51h feeding VS Two 51h control
2. There are currently 4 different DEgenes pipelines having been tried out
   1. DESeq2
   2. EdgeR
   3. Salmon-sleuth
   4. Kallisto-Sleuth
3. There are currently 2 different annotation method being tried: 
   1. Blastx against mus musculus
   2. Orthofinder Ortholgroups
4. There are currently two ways for further curation
   1. Inspect the overlapping of DGE with reported regeneration programs from literatures
      1. science S2: 310 shared RRE-regulated genes during killifish and zebrafish regeneration.    
      2. science S3: 546 down-regulated genes shared between killifish and zebrafish (RNA-seq).    
      3. science S4: The regeneration response program is composed of 49 shared genes with RREs and elevated gene expression.  A gene with multiple copies in the genome is only counted as one gene.    
      4. science S7: 630 blastema marker genes
   2. gProfiler for ontology analysis

#### Computational steps towards comprehensive analysis

1. Generate DE groups for all 5 comparison groups x 4 pipelines, then add the 2 annotation columns, and the science paper curation
2. Inspect the Venn diagram for the overlapping extent of the four method for each comparison (down and up regulation separately)
3. For the most overlapped (conserved) genes, inspect the overlapping across conditions
   1. whether stimulated VS non-stimulated at 27h and 51h differs a lot and how about the wrapped-together comparison (3 VS 5)?
   2. how does the amputation DE gens overlap with the "regenerating DE genes" ?
   3. does human curation on the 2 feeding 27h samples lead to great differences in terms of DE genes list?
   4. if so what are some most most confident DE genes are what are there gProfiler/science annotation state

4. Plot some heatmaps based on the most confident genes / curated genes overlapping with science paper

   

#### Code and Results