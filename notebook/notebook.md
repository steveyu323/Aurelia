## 200105

### Literature Reading 

1) Changes in regeneration-responsive enhancers shape regenerative capacities in vertebrates

- The conserved response revealed several regeneration-responsive enhancers (RREs), including an element upstream to inhibin beta A (inhba), a known effector of vertebrate regeneration

- This element activated expression in regenerating transgenic fish, and its genomic deletion perturbed caudal fin regeneration and abrogated cardiac regeneration altogether

- The enhancer is present in mammals, shares functionally essential activator protein 1 (AP-1)–binding motifs, and responds to injury, but it cannot rescue regeneration in fish. This work suggests that changes in AP-1–enriched RREs are likely a crucial source of loss of regenerative capacities in vertebrates.

- Active enhancers and promoters are charac- terized by histone H3K27ac and H3K4me3 marks. We assayed both killifish and zebrafish genomes (~1.5 gigabases) for H3K27ac and H3K4me3 enrichment using chromatin immunoprecipitation sequencing (ChIP-seq) in samples of uninjured (0 dpa) and regener- ating (1 dpa) caudal fin.

- This shared cohort encompasses several known and essential regulators of zebrafish regeneration, includ- ing fgf20a, inhbaa, junbb, and fn1, as well as putative new regulators such as crlf1, vmp1, and tgfbr1.

- We conclude that AP-1 motifs are required for the activation of RREs in response to amputation.

- AP-1 components were not ubiquitously expressed in all cell types and were present in both mesenchymal and epithelial cells, suggest- ing that specific subunit compositions may be required to restrict expression of the enhancers to either the mesenchyme (K-IEN, Z-IEN) or epithelial cells (H-IEN).

- Additionally, the pres- ence of a predicted p53/p63–binding motif in the human enhancer, which was absent in K-IEN and Z-IEN (fig. S23A), and the abun- dance of p63 expression in basal epidermal cells (fig. S23B) suggest that interactions of AP-1 with other nuclear factors may also play a role in regulating enhancer activity.

- The Drosophila ATAC-seq dataset (Accession number: GSE102841) for the regeneration of wing imaginal discs were obtained from a previously published study (31).

- The RNA-seq data of ear pinna injury (GEO accession: GSE71761) and skin injury (GEO accession: GSE113081) in the mice M. musculus and the African spiny mice A. cahirinus were obtained from previously published studies (23, 24)

- Seurat, HOMER



### Data Grabbing

- The Drosophila ATAC-seq dataset (Accession number: GSE102841) for the regeneration of wing imaginal discs were obtained from a previously published study (31).

- The RNA-seq data of ear pinna injury (GEO accession: GSE71761) and skin injury (GEO accession: GSE113081) in the mice M. musculus and the African spiny mice A. cahirinus were obtained from previously published studies (23, 24)



## 201008

1. Email Prof.Geontoro for jellyfish sequencing result and paper
2. Paper Reading : Advances in understanding tissue regenerative capacity and mechanisms in animals
   1. Regeneration genes: Prod1 in salamander
   2. certain gene program selectively activated: zebrafish JunB phosphorylated by JNK in fin generation
   3. Little is known about the molecular underpinnings of the inverse relationship between age and regenerative capacity. However, recent studies have made a compel-ling argument that age-related increases in p16INK4a levels impede regeneration in several tissues30–32. Furthermore, ageing skeletal muscle seems to respond positively to a blood-borne factor (or factors) that is produced in young animals33, and recent publications point to Wnt signalling effectors34,35
   4. Initiation and tageting of regeneration: For example, in regener-ating zebrafish fins, these signals include Fgf20a, cer-tain Wnt ligands and activin-βA; Here, as these cells die, they synthesize and release Wnt3, a factor that rapidly stimulates local cell proliferation
3. Paper Reading: The Molecular and Cellular Choreography of Appendage Regeneration
   1. canonical WNT, retinoic acid & FGF20, Macrophage
   2. in newt, mammalian MLP intracellular localization
   3. satellite stem cell driven regeneration is a conserved mechanism
   4. dedifferentiation of muscle cell is another conserved mechanism
   5. late mode of limb budding could be ancestral and important for regenerative capacity, but not adopted by most animals
4. https://regendbase.org/experiments



## Aurelia Genome:

**Gold D.A.** & Katsuki T. (co-first authors), Li Y., Yan X., Regulski M., Ibberson D., Holstein T., Steele R.E., Jacobs D.K., Greenspan R.J. (2019) The genome of the jellyfish *Aurelia* and the evolution of animal complexity. *Nature Ecology and Evolution*. 3(1):96-104. DOI Link: [10.1038/s41559-018-0719-8](https://www.nature.com/articles/s41559-018-0719-8.epdf?author_access_token=HA9z7Wq_3-B61nnfCQusFdRgN0jAjWel9jnR3ZoTv0OoS3jcjhVPLlwHthSCDxdXSIkqCsrq8xdtwxeS31yQxr5hNE09h-d3oFfE7-Y6q_vb9yHmpHMpbEHBNC1pxuFZL8wQNKMOS7pMvv7Lqnidsw%3D%3D)  

- [https://davidadlergold.faculty.ucdavis.edu/jellyfish/]
- https://drive.google.com/drive/folders/1NC6bZ9cxWkZyofOsMPzrxIH3C7m1ySiu





### STAR Alignment

https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html

```
scp -r cyu7@login.hpc.caltech.edu:/central/groups/GoentoroLab/aohdera/03-TRIM ./Desktop/  

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install star
conda update star

gffread my.gff3 -T -o my.gtf

scp -r ./Desktop/genome cyu7@login.hpc.caltech.edu:/home/cyu7/jellyfishy/genome_data/
scp -r ./Desktop/trim_out cyu7@login.hpc.caltech.edu:/home/cyu7/jellyfishy/

submit to slurm job

!!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=757170055, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 13



Genome sequence total length = 757170055
Genome size with padding = 6784024576
Jan 14 00:58:47 ..... processing annotations GTF
Processing pGe.sjdbGTFfile=./genome_raw/aurelia_genome_annotation.gtf, found:
                35273 transcripts
                145435 exons (non-collapsed)
                93402 collapsed junctions
Total junctions: 93402
Jan 14 00:58:48 ..... finished GTF processing
Estimated genome size with padding and SJs: total=genome+SJ=6985024576 = 6784024576 + 201000000
GstrandBit=33
Number of SA indices: 1319747688
Jan 14 00:59:14 ... starting to sort Suffix Array. This may take a long time...
Number of chunks: 7;   chunks size limit: 1743195056 bytes
Jan 14 00:59:34 ... sorting Suffix Array chunks and saving them to disk...



mkdir a_63
mv *.out a_63/
mv Aligned.out.sam a_63/
mv SJ.out.tab a_63/


STAR --genomeDir ./aurelia_genome_index \
 --runThreadN 12 \
 --readFilesIn ./trim_out/SRR12953064_1_trimmed.fq

mkdir a_64
mv *.out a_64/
mv Aligned.out.sam a_64/
mv SJ.out.tab a_64/

STAR --genomeDir ./aurelia_genome_index \
 --runThreadN 12 \
 --readFilesIn ./trim_out/SRR12953065_1_trimmed.fq

mkdir a_65
mv *.out a_65/
mv Aligned.out.sam a_65/
mv SJ.out.tab a_65/

$grep -B4 -A8 "@SRR12953062.7870335 HISEQ:622:H3NMJBCXY:2:1113:16227:42625 length=100
TCTGGATGTGTGATATGATGATTTGGTTTGGGACTGCGGGCTAGAACGATCATTGGTTG" ./trim_out/SRR12953062_1_trimmed.fq

cat ./trim_out/SRR12953062_1_trimmed.fq | paste - - - - | awk -F '\t' '(length($2)!=length($4))'

cat ./trim_out/SRR12953062_1_trimmed.fq  | paste - - - - | awk '(length($2)==length($4))' | tr "\t" "\n"


samtools view -Su Aligned.out.sam | samtools sort - sorted_56




cd a_64 
samtools view -Su Aligned.out.sam | samtools sort -o sorted_64.bam
cd ..
mv a_64/sorted_64.bam ./stringtie_input/sorted_64.bam

cd a_65
samtools view -Su Aligned.out.sam | samtools sort -o sorted_65.bam
cd ..
mv a_65/sorted_65.bam ./stringtie_input/sorted_65.bam


stringtie sorted_56.bam  -e -o ../out_56/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
stringtie sorted_57.bam  -e -o ../out_57/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
stringtie sorted_58.bam  -e -o ../out_58/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
stringtie sorted_59.bam  -e -o ../out_59/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
stringtie sorted_60.bam  -e -o ../out_60/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
stringtie sorted_61.bam  -e -o ../out_61/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
stringtie sorted_62.bam  -e -o ../out_62/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
stringtie sorted_63.bam  -e -o ../out_63/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3
stringtie sorted_64.bam  -e -o ../out_64/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
stringtie sorted_65.bam  -e -o ../out_65/out.gtf -p 6 -G /central/groups/GoentoroLab/YuC/jellyfishy/genome_raw/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 
```





## 210122

### Paper Reading

1. **Effect of a leucine-enriched essential amino acids mixture on muscle recovery**
   1. leucine-enriched es- sential amino acid mixture (LEAA) consumption suppressed exercise-induced elevation of muscle damage markers in blood, which suggests that LEAA could attenuate muscle damage and aid muscle recovery.
   2. Namely, mTOR promotes muscle regeneration through kinase-independent and kinase-dependent mechanisms at the stages of nascent myofiber formation and myofiber growth, respectively
   3. chain amino acids (BCAAs) increases the anabolism and decreases the catabolism of muscle proteins
   4. Leucine, an EAA, activates mTOR signaling pathway18) and has a key role in the initiation of muscle protein synthe- sis19–26)

2. **Genome and single-cell RNA-sequencing of the earthworm Eisenia andrei identifies cellular mechanisms underlying regeneration**
   1. Using Ki-67 immunofluorescent label- ing, we found that cell proliferation initiated at 24h post- amputation, and at 48 and 72 h post-amputation the proliferating cells increased rapidly and gradually migrated to the center of cross sections (Fig. 2d and Supplementary Fig. 6). At 5 days post- amputation, the wound healing was fully accomplished and a small blastema appeared in center of the amputation plane (Supplementary Fig. 4). At 6 and 7 days post-amputation, the blastema persistently experienced growth and elongation (Sup- plementary Fig. 4). Although the newly produced body segments were not observed at 14 days post-amputation, the base of out- growth has accumulated pigments (Supplementary Fig. 4). At 18 days post-amputation, new body segments arise, and at 28 days post-amputation the obvious body segments take shape in regenerative appendages 
   2.  In total, 6,048 DEGs that changed their expression at one or more regeneration time points were identified, and these genes demonstrated a temporal order in their expression profiles (Supplementary Fig. 8). Gene enrichment analysis found that many biological processes important for development were commonly upregulated across all regeneration stages, including gene transcription (GO:0006351), Wnt signaling pathway (GO:0016055), cell surface receptor signaling pathway (GO:0007166), multicellular organism devel- opment (GO:0007275), and anatomical structure development (GO:0048856) (Supplementary Data 1). 
   3. The neighboring genes of 19 DEL2s, such as EGR1, FOSL, BMP10, HUNB and MMP17, are frequently reported to participate in regeneration
   4. Furthermore, we performed a randomization test and found five gene families standing out as showing significantly higher proportion of differentially expressed genes (P < 0.05, χ2 test), including ZNFX1, EGFR, NNP, HELZ2 and SACS
   5. To understand the large-scale gene interactions involved in regeneration, we conducted a weighted gene coexpression network analysis (WGCNA)
   6.  Previous studies reveal that FOS participates in neoblast maintenance and the wound response program in planarians36,37 and is a key factor in the cell signaling system activated immediately after cell damage
   7. Among them, the expression of the brown module was most significantly correlated with the regeneration stage (6 h) (r = 0.53, P = 0.003) (Fig. 5b and Supplementary Fig. 23). Genes enriched in this module participate in signal transduction, transcription and translation, implying an increas- ing level of cell communication and biochemical processes via the synthesis of mRNAs and proteins in response to regeneration (Supplementary Data 4). 
   8. The list of driver genes in the brown module triggered by the regeneration process includes several genes involved in cellular proliferation, differentiation and programmed cell death, such as FOS (intramodule member- ship = 0.9587) and HUNB (intramodule membership=0.934 (Fig. 5c, and Supplementary Table 11).
   9. Two other modules, red and blue, containing genes with increased expression until 12 h of regeneration (Fig. 5d, e, and Supplementary Fig. 24), were also enriched in genes involved in biosynthetic processes and the regulation of cell growth (Supplementary Data 5 and 6). Additionally, the blue module was enriched in genes involved in energy metabolism that are necessary for cell proliferation and growth. 
   10. Gene ontology enrichment analysis. Gene functional enrichments at three levels (biological process, molecular function and cellular component) were performed using a web-based gene analysis tool, g:Profiler (rev1705) (http://biit.cs.ut.ee/ gprofiler/). The p-value was adjusted by Benjamini-Hochberg FDR.
   11. Identification of coexpression networks in early regenerative processes. Analysis was carried out in R on a 64-bit LINUX platform with 65.7 GB memory. Modules/or networks were constructed using WGCNA34 (v1.67). Modules were defined as branches of the hierarchical cluster tree using the dynamic tree cut method. For each module, the expression patterns were summarized by the module eigengene (ME), defined as the right singular vector of the standardized expression patterns. MEs were also defined as the first principal component calculated using PCA, which can summarize module behavior. Pairs of modules with high ME correlations (R > 0.8) were merged. MEs for modules were plotted by using the ggplot2 library in R. These MEs were tested for correlation with phenotypes (regeneration time points) adjusted by a linear regression model. In more detail, a weighted signed network was computed based on a fit to scale-free topology, with a threshold softPower of 10 chosen (as it was the smallest that resulted in a scale-free R2 fit). A topological overlap dendrogram was used to define modules with a minimum module size of 80 genes and the deepSplit parameter set to 2. The connectivity of every gene in every module was assessed by correlation to the MEs, or kMEs. Module membership (MM) was regarded as intramodular connectivity. MM can be combined as a systematic biological method to obtain driver genes in networks, which are highly interconnected nodes within coexpression gene mod- ules. The driver genes were defined by the WGCNA connectivity algorithm. Each module network was viewed by VisANT (v5.0) (http://www.visantnet.org/ visantnet.html), which allows users to input an edge file and a node file from a WGCNA module.



BLASTx

wget https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl

update_blastdb.pl --decompress nr

export PATH=$PATH:/central/groups/GoentoroLab/YuC/jellyfishy/blast/ncbi-blast-2.11.0+/bin



blastx -db ./blastdb/nr \
 -query test.fasta \
 -taxids 10090 \
 -outfmt 7 \
 -out OUTPUT_test.tab \
 -max_target_seqs 1 \
 -show_gis \
 -num_threads 12

blastx -db ./blastdb/nr \
 -query diff_seg.fasta \
 -taxids 10090 \
 -outfmt 7 \
 -out OUTPUT_diff.tab \
 -max_target_seqs 1 \
 -show_gis \
 -num_threads 12



blastx -db ./blastdb/nr \
 -query /central/groups/GoentoroLab/YuC/jellyfishy/data/genome_raw/Aurelia.V.1.1.transcripts.fasta \
 -taxids 10090 \
 -outfmt 7 \
 -out OUTPUT_diff.tab \
 -max_target_seqs 1 \
 -show_gis \
 -num_threads 12

https://arxiv.org/abs/1611.07308

https://arxiv.org/pdf/1606.05908.pdf

https://github.com/zfjsail/gae-pytorch

https://github.com/pytorch/examples/tree/master/vae





### 210131

- Rerun blastx by dividing the transcript fasta file into 2 parts to maximize the parallelization

  ```
  scp cyu7@login.hpc.caltech.edu:/central/groups/GoentoroLab/YuC/jellyfishy/blast/OUTPUT_all_1.tab ./Desktop/
  scp cyu7@login.hpc.caltech.edu:/central/groups/GoentoroLab/YuC/jellyfishy/blast/OUTPUT_all_2.tab ./Desktop/
  
  ```

  


### 200202

#### TO-DOs after group meeting 

1. use a log scale and then normalize the samples based on library size to perform a better PCA
2. Get the data that Aki mentioned and to build modules for VAE/GAE
3. Run DESeq2 pipeline to see if there are any significant differences
4. Finalize the suimmary of g:Profiler results .









### Installation

```
$ git clone https://github.com/AntixK/PyTorch-VAE
$ cd PyTorch-VAE
$ pip install -r requirements.txt
```





The VAE/GAE Design Doc for learning conserved genetic elements through multi-organism regeneration RNA-seq data

1. Unlike scVI or other single cell based program, each gene rather than each cell should be treated as a data. Or conversely, each sample should be treated as a datapoint
2. We should collect as much data and as diverse data as possible to better generalize the model
3. State of art VAE implementation is ready to use
4. For GAE. prior knowledge on pathway/GRN should be applied to construct adjacency matrix



1. Data Collection:

   1) Jellyfish amputation data from David [with feeding/ non-feeding condition] (10 samples)

   2) 



**Deep evolutionary origin of limb and fin regeneration**

SRR2885871 RNA-seq of axolotl nonregenerating upper armtissue

SRR2885875 RNA-seq of axolotl nonregenerating upper armtissue

SRR2885873 RNA-seq of axolotl nonregenerating upper armtissue

SRR2885866  axolotl proximal blastema 

SRR2885865  axolotl proximal blastema 

AM_NRF01 non-regenerating fin 

AM_NRF02 non-regenerating fin 

AM_NRF03 non-regenerating fin 

AM_BF01 blastema fin

AM_BF02 blastema fin

AM_BF03 blastema fin



**Genome and single-cell RNA-sequencing of the earthworm *Eisenia andrei* identifies cellular mechanisms underlying**

https://www.nature.com/articles/s41467-020-16454-8#Sec27



TO-DO: The projects and and sequencing repository for different projects are messily documented and so it is slightly hard to actually compile enough samples out of it for sufficient learning models. Instead, it would be of high interests to perform a cross-species regeneration single cell experiment based on a set of significantly differentially expressed genes/Seurat's program on integration. With each cell as a node, or with a graphical representation of genes 







### 200215 

1. Finish Summary of Comparison between EdgeR and DESeq2 pipeline
2. Compile all the annotations from other studies
   1. Al
   2. nature commu:  GO programs shared during regeneration
   3. science S2: 310 shared RRE-regulated genes during killifish and zebrafish regeneration.
   4. science S3: 546 down-regulated genes shared between killifish and zebrafish (RNA-seq).
   5. science S4: The regeneration response program is composed of 49 shared genes with RREs and elevated gene expression.  A gene with multiple copies in the genome is only counted as one gene.
   6. science S7: 630 blastema marker genes
3. Heatmap for curated/naturally ranked set of genes



## 200310

### Paper Reading:

1. A novel approach to comparative RNA-Seq does not support a conserved set of genes underlying animal regeneration
   1. Orthofinder for ortholog gene mapping
   2. Download the compiled available dataset to cluster for future reference
   3. 160 deCOGs overlap with data
2. Lior's data on starvation:
   1. overlapping pathway with out DEG in different cell types
   2. Some questions related to practical conduct of future experiment design:
      1. #number of cells in one Aurelia
      2. #conditions to be labelled
      3. identification of regeneration markers



## 200402

1. Check specific pathway environment:
   1. Arginine-decomposition
2. More DEG pipeline run through via HMS's guide
3. Transcript calling using different approaches including kallisto and salmon to compare the reliability of STAR mapping gene count
4. Aligning the blastx directly to human/Using orthofinder to construct Aureia-Mouse-Human ortholog gene sets
5. Paper Reading:
   1. Gamma-delta cells regulate the intestinal response to nutrient sensing
      1. cell markers defined for epithelial cells, stromal cells, and immune cells
      2. removal of ribosomal, mitochondrial, immunoglobulinm and HLA genes to reduce unwnated batch effects
      3. ComBat for batch correction
      4. Dirchlet-multinomial regression model, tests for differences in cell composition between conditions while accounting for the proportions of all of the other cell subsets
   2. An early cell shape transition drives evolutionary expansion of the human forebrain
      1. The ability to generate brain organoids from induced pluripo- tent stem cells (iPSCs) has enabled researchers to study a vari- ety of neurodevelopmental processes that were previously inac- cessible
      2. In total, we performed RNA-seq of 42 samples of $9,000 organoids, collected from day 0, correspond- ing to pre-neurulation tissue, up until day 25, corresponding to fully committed neurogenesis
      3. 1) genes were expressed at > 10 TPM in at least one time point in all replicates of both species, 2) genes displayed a fold-change of > 1.5 TPM between any two time points in at least one species. This resulted in a list of 3,526 genes. 
      4. Log2 normalized TPMs were then z-scaled across the three time points. Next, the gene list was further filtered in order to remove genes with variable expression patterns across replicates by removing genes with a squared difference > 6 in either species. The squared difference value was obtained by calculating the squared difference of z-scores per time point between all three replicates and then taking the sum of this number for all time points. 
      5. This list of 2,905 genes was used for time course sequencing data analysis using the TCseq package (Wu and Gu, 2020) available online (https://rdrr.io/bioc/TCseq/f/ inst/doc/TCseq.pdf). Replicates were kept separate for the analysis, meaning that the input was 17,430 unique patterns of gene expression representing 2,905 genes in 2 species and 3 replicates per species. The TCSeq analysis was run using the *timeclust* func- tion with the settings *algo = ‘cm’, k = 10*, resulting in the unsupervised soft clustering of gene expression patterns into 10 clusters with similar z-scaled temporal patterns. Replicates of 563 genes were found to be in 3 different clusters in at least one species and were removed from downstream analysis resulting in 2,342 genes (Data S1). 22.5% of genes, 527 genes, were found to never be in the same cluster between species showing a robustly different expression pattern (Data S1). Genes were assigned to the cluster where 2 or all of the replicates were found per species. 59% of genes were assigned to different clusters between species. GO term enrich- ment analysis was performed on the genes present in each of the 10 clusters generated by TCseq per species, resulting in a total of 159 GO:BP terms found in both species, of which 85 were found to be moving between species (Data S1). GO term analysis on the 527 genes moving robustly between species revealed 67 were linked to enriched cell morphogenesis-related terms (‘‘cell morpho- genesis,’’ ‘‘cell part morphogenesis,’’ ‘‘cellular component morphogenesis’’). 
      6. These genes were intersected with a list of 1,639 confirmed transcription factors (Lambert et al., 2018), revealing 8 transcription factors related to cell morphogenesis

## Orthofinder

Installation of orthofinder`conda install orthofinder`

Running orthofinder on testing dataset ` orthofinder -f ExampleData/`

To construct the ortholog groups on drosophila, mouse, human, and Aurelia, we need to first download the corresponding proteome

Human proteome: 

`wget http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz`

Mouse proteome:

`wget http://ftp.ensembl.org/pub/release-103/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz`

Drosophila proteome:

`wget http://ftp.ensembl.org/pub/release-103/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.32.pep.all.fa.gz`

Aurelia proteome: (Using Daivid Gold compiled version )

`cp /central/groups/GoentoroLab/YuC/jellyfishy/data/genome_raw/Aurelia.V.1.1.proteins.fasta /central/groups/GoentoroLab/YuC/jellyfishy/proteomes/Aurelia.V.1.1.proteins.fasta`

The Aurelia proteome contains '.' in the fasta file (due to potential alignment reason?). Dots in the sequences should be removed before running through orthofinder

`sed '/^[^>]/s/\.//g' Aurelia.V.1.1.proteins.fa > Aurelia.V.1.1.proteins.cleaned.fa`

```bash
gunzip *.gz
for f in *fa ; do python /central/groups/GoentoroLab/YuC/jellyfishy/OrthoFinder/tools/primary_transcript.py $f ; done
```

```
orthofinder -f primary_transcripts/
```

Results are stored in `/central/groups/GoentoroLab/YuC/jellyfishy/proteomes/primary_transcripts/OrthoFinder/Results_Apr05` after completion of the run



## Salmon

#### Installation:

``` bash
conda create -n salmon salmon
conda deactivate salmon
```

#### To Use Salmon 

```bash
mkdir /central/groups/GoentoroLab/YuC/jellyfishy/salmon
cd /central/groups/GoentoroLab/YuC/jellyfishy/salmon
conda activate salmon

# build salmon index from transcriptome file from David Gold
salmon index -t /central/groups/GoentoroLab/YuC/jellyfishy/data/genome_raw/Aurelia.V.1.1.transcripts.fasta -i aurelia_index

# quantify the transcripts of regeneration/starvation samples
#!/bin/bash
for fn in /central/groups/GoentoroLab/YuC/jellyfishy/data/trim_out/SRR129530{56..65}_1_trimmed.fq;
do
samp=`basename ${fn}`
echo "Path ${fn}"
echo "Processing sample ${samp}"
salmon quant -i aurelia_index -l A \
				-r ${fn} \
				--numBootstraps 10 \
        -p 8 --validateMappings -o quants/${samp}_quant
done 



```

### Download the Salmon output to local environment

`scp -r cyu7@login.hpc.caltech.edu:/central/groups/GoentoroLab/YuC/jellyfishy/salmon/quants  ./Desktop/salmon`

### Use R package Wasabi for converting Salmon output to be suitable for Sleuth analysis

- #62 salmon output was corrupted, to fix the issue, re-run salmon and inspect the log file

- It seems that Aki's trimming method does not generate the error, and so for now using Aki's trimmed data as input into salmon

- ```bash
  #!/bin/bash
  for fn in /central/groups/GoentoroLab/aohdera/03-TRIM/SRR12953062_trim.fastq.gz;
  do
  samp=`basename ${fn}`
  echo "Path ${fn}"
  echo "Processing sample ${samp}"
  salmon quant -i aurelia_index -l A \
  				-r ${fn} \
  				--numBootstraps 10 \
          -p 8 --validateMappings -o quants/${samp}_quant
  done 
  ```

- `scp -r cyu7@login.hpc.caltech.edu:/central/groups/GoentoroLab/YuC/jellyfishy/salmon/quants  ./Desktop/salmon`



### Kallisto

```{bash}
~/miniconda3/bin/kallisto index -i transcripts.idx /central/groups/GoentoroLab/YuC/jellyfishy/data/genome_raw/Aurelia.V.1.1.transcripts.fasta

#!/bin/bash
for fn in /central/groups/GoentoroLab/YuC/jellyfishy/data/trim_out/SRR129530{56..65}_1_trimmed.fq;
do
samp=`basename ${fn}`
echo "Path ${fn}"
echo "Processing sample ${samp}"
~/miniconda3/bin/kallisto quant -i transcripts.idx -o ${samp}_quant -b 100 --single -l 180 -s 20 ${fn}
done 


```



`scp -r cyu7@login.hpc.caltech.edu:/central/groups/GoentoroLab/YuC/jellyfishy/kallisto/quants  ./Desktop/kallisto`



## 210421

`scp -r cyu7@login.hpc.caltech.edu:/central/groups/GoentoroLab/YuC/jellyfishy/proteomes/primary_transcripts/OrthoFinder/Results_Apr05 ./Desktop/orthofinder`

### Orthofinder Result Interpretation

- Number of species       4
- Number of genes 95194
- Number of genes in orthogroups  79797
- Number of unassigned genes      15397
- Percentage of genes in orthogroups      83.8
- Percentage of unassigned genes  16.2
- Number of orthogroups   19465
- Number of genes in species-specific orthogroups 23117
- Percentage of genes in species-specific orthogroups     24.3
- *The Orthogroups are stored in Orthogroups/Orthogroups.tsv*