# 16s_microbiome_analysis_workflow
### An example to analyze 16s microbiome sequences for taxonomy and functional prediction.
note: some codes can be simplified. There are several repeated steps that can be combined if run this workflow in order.
(i.e. use `asv_taxa.biom` as input to PICRUSt2 might generate stratified gene tables with taxonomic annotations.)
1. required package:
  * Python: `scipy`, `pands`, `statsmodels`
  * R: `vegan`, `DADA2`, `ALDEx2`, `biomformat`, `ComplexHeatmap`, `RColorBrewer`
  * Commandline tools: `PICRUSt2`, `QIIME2`, `HUMAnN2`,`biom-format`
2. workflow:
  * the following workflow is designed for paired-end sequencing
  * replace `~/.../` to actual working paths
  * Example datasets from `Young versus Aged Microbiomes in Germ-Free Mice: Increased Short-Chain Fatty Acids and Improved Cognitive Performance`
  * rename raw sequences files to `forward/reverse/barcodes.fastq.bz2`
  * the following is down in shell/commandline. (Using Terminal in Mac) make sure that `qiime2` is in the current working environment
  * decompress/recompress files for qiime2
    ```
    cd ~/.../ #raw data file path
    bzip2 -d *.bz2 #decompress
    gzip *.fastq #recompress
    ```
  * demultiplexing sequences using qiime2 (this step takes time)
    ```
    qiime tools input \
    --type EMPPairedEnd Sequences \
    --input-path ~/.../ \
    --output-path ~/.../xx.qza
    cd ~/.../ #qiime input data path (output path from last command)
    qiime demux emp-paired \
    --i-seqs paired_CA.qza \
    --m-barcodes-file barcodes.txt \
    --m-barcodes-column BarcodeSequence \
    --o-per-sample-sequences demux.qza \
    --o-error-correction-details demux-details.qza
    #extracting data for dada2_workflow.R
    mkdir ~/../qiime2_data #where you want to store extracted data
    qiime tools extract \
    --input-path demux.qza \
    --output-path ~/../qiime2_data
    ```
    * outputs include a folder of demultiplexed sequences with name `*_R1_001.fastq.gz` and `*_R2_001.fastq.gz`. This folder path will be the input for next step with DADA2.
  * `dada2` (in `R/RStudio`) denoising and taxonomy assignment;
    * we used `SILVA-v138` database for taxonomy assignment. Please check DADA2 website to download the newest version of SILVA or other database.
    * `dada2_workflow.R` see DADA2 website for full description
    * outputs include `asv.biom`, `asv.fna` and `asv_genus.tsv` files for taxonomy and functional analysis
  * taxonomy analysis. this workflow used `centered log-ratio` transformation instead of traditional `rarefaction` and `relative abundance`
    * `ALDEx2` and `vegan` (in `R/RStudio`)
    * preprocessing:
      * convert `asv.biom` to `asv.tsv` (in commandline)
      ```
      biom convert -i asv.biom -o asv.tsv --to-tsv
      ```
      * group ASVs by taxonomy: `asv_preprocessing.ipynb` (in `jupyter notebook`)
    * `taxonomy_analysis.R`
  * functional analysis using `PICRUSt2`(in `commandline`). make sure `PICRUSt2` and `HUMAnN2` is in the current environment
    ```
    cd ~/.../dada2_output #the path contains dada2 outputs
    picrust2_pipeline.py \
    -s asv_all.fna \
    -i otu_all.biom \
    -o ../picrust2_out \
    -p 1 \
    --stratified \ #give out stratified results
    --coverage \  #give out pathway coverage as in humann2
    --metagenome_contrib #contribution of each ATV
    ```
    ```
    # to normalize counts to relative abundabce
    cd ~/.../picrust2_out/ #path contains PICRUSt2 outputs
    # PICRUSt2 output is gzip compressed. either decompress or directly input .gz files
    humann2_renorm_table \
    -i EC_metagenome_out/pred_metagenome_unstrat.tsv \
    -u relab \
    -o EC_metagenome_out/ecs_relab.tsv
    humann2_renorm_table \
    -i EC_metagenome_out/pred_metagenome_strat.tsv \
    -u relab \
    -o EC_metagenome_out/ecs_relab_strat.tsv
    humann2_renorm_table \
    -i pathways_out/path_abun_unstrat.tsv \
    -u relab \
    -o pathways_out/pathabun_relab.tsv
    ```
    ```
    # PICRUSt2 function
    add_descriptions.py \
    -i EC_metagenome_out/ecs_relab.tsv \
    -o EC_metagenome_out/ecs_relab_names.tsv \
    -m EC
    add_descriptions.py \
    -i pathways_out/ecs_relab_strat.tsv \
    -o pathways_out/ecs_relab_strat_names.tsv \
    -m METACYC
    add_descriptions.py \
    -i pathways_out/pathabun_relab.tsv \
    -o pathways_out/pathabun_relab_names.tsv \
    -m METACYC
    ```
    * statistical analysis for unstatified `PICRUSt2` outputs: `gene_picrust2_stats.ipynb`, `pathway_picrust2_stats.ipynb` (in `jupyter notebook`) and `p_adjustmenet.R` in (`R/RStudio`) (optional, see notes in `gene_picrust2_stats.ipynb`); for stratified contribution of certain genes (i.e. SCFAs): `species_contribution.ipynb` (in `jupyter notebook`) and `p_adjustmenet.R` in (`R/RStudio`)
    * visualize selected enzymes by heatmap `picrust_heatmap.R` (note: the input here are the average values of scaled relative abundance of each enzyme at each time points (used `scale()` function in `R` the process is not included here)).
3. future direction: using `Bokeh` or `plotly` for interactive plotting.
4. citations:
  * `Virtanen, Pauli, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski et al. "SciPy 1.0: fundamental algorithms for scientific computing in Python." Nature methods 17, no. 3 (2020): 261-272.`
  * `Jeff Reback, Wes McKinney, jbrockmendel, Joris Van den Bossche, Tom Augspurger, Phillip Cloud, gfyoung, et al. “Pandas-dev/pandas: Pandas 1.0.3”. Zenodo, March 18, 2020. doi:10.5281/zenodo.3715232.`
  * `Seabold, Skipper, and Josef Perktold. “statsmodels: Econometric and statistical modeling with python.” Proceedings of the 9th Python in Science Conference. 2010.`
  * `Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre
  Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter
  Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2019). vegan:
  Community Ecology Package. R package version 2.5-6.
  https://CRAN.R-project.org/package=vegan`
  * `Callahan, Benjamin J., Paul J. McMurdie, Michael J. Rosen, Andrew W. Han, Amy Jo A. Johnson, and Susan P. Holmes. "DADA2: high-resolution sample inference from Illumina amplicon data." Nature methods 13, no. 7 (2016): 581.`
  * `Fernandes, AD, Macklaim, JM, Linn, TG, Reid, G, Gloor, GB (2013). “ANOVA-Like Differential Gene Expression Analysis of Single-Organism and Meta-RNA-Seq.”PLoS ONE, 2013, volume 8, issue 7, e67019. http://dx.doi.org/10.1371%2/journal.pone.0067019`
  * `Fernandes, D. A, Reid, J., Macklaim, M. J, McMurrough, T.A, Edgell, D.R., Gloor, B. G (2014). “Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis.” Microbiome, 2014, volume 2, 15. http://doi:10.1186/2049-2618-2-15`
  * `Gloor GB, Macklaim JM, Fernandes AD (2016). “Displaying Variation in Large Datasets: a Visual Summary of Effect Sizes.” Journal of Computational and Graphical Statistics. http://dx.doi.org/10.1080/10618600.2015.1131161`
  *  `Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in
  multidimensional genomic data. Bioinformatics.`
  * `Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer`
  * `Douglas, Gavin M., Vincent J. Maffei, Jesse Zaneveld, Svetlana N. Yurgel, James R. Brown, Christopher M. Taylor, Curtis Huttenhower, and Morgan GI Langille. "PICRUSt2: An improved and extensible approach for metagenome inference." BioRxiv (2019): 672295.`
  * `Franzosa, Eric A., Lauren J. McIver, Gholamali Rahnavard, Luke R. Thompson, Melanie Schirmer, George Weingart, Karen Schwarzberg Lipson et al. "Species-level functional profiling of metagenomes and metatranscriptomes." Nature methods 15, no. 11 (2018): 962-968.`
  * `Bolyen, Evan, Jai Ram Rideout, Matthew R. Dillon, Nicholas A. Bokulich, Christian C. Abnet, Gabriel A. Al-Ghalith, Harriet Alexander et al. "Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2." Nature biotechnology 37, no. 8 (2019): 852-857.`
  * `Paul J. McMurdie and Joseph N Paulson (2019). biomformat: An interface package
  for the BIOM file format. https://github.com/joey711/biomformat/,
  http://biom-format.org/.`
