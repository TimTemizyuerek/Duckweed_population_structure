     ## load packages ####
     
     library(geosphere) # to compute geographical distances
     library(ape) # tree
     library(adegenet) #PCA
     library(LEA) ## for population structure analysis
     library(RColorBrewer) ## ready made color paletes
     library(vegan) ## for mantel()
     
     # install_version("SNPfiltR", "1.0.1")
     library(snpStats) ## for linkage disequilibrium
     library(SNPfiltR)
     library(vcfR)
     library(dartR)
     
     ## custom functions #### 
     
     ## convert GPS to format that facilitates distance calculations
     convert_dmm_to_dd = function(coord_str) {
          parts = strsplit(coord_str, "°")[[1]]
          return(as.numeric(parts[1]) + as.numeric(parts[2]) / 60)
     }
     
     ## transform hamdist to rarecurve friendly dataframe
     hamdist_to_rarecurve = function(hamdist_subset) {
       
       ## set up loop
       clone_list = list()
       whiler=TRUE; n=0
       while (whiler == TRUE){
         
         ## run through matrix by cols
         n=n+1
         runner_col = hamdist_subset[,n,drop=FALSE]
         
         ## are there clones?
         if (length(which(runner_col<0.02)) > 1) {
           
           ## if yes save clones and remove from hamdist
           clone_list[[n]] = rownames(hamdist_subset)[which(runner_col<0.02)]
           hamdist_subset = hamdist_subset[-which(rownames(hamdist_subset) %in% clone_list[[n]][-1]),-which(colnames(hamdist_subset) %in% clone_list[[n]][-1]), drop=FALSE]
           
         } else {
           
           ## if not, just save individual as only representative of that clone
           clone_list[[n]] = colnames(runner_col)
         }
         
         ## stop the loop of hamdist runs out
         if (n == ncol(hamdist_subset)) {
           whiler=FALSE
         }
         
         
       }
       
       ## transform hamdist to rarecurve format
       for (n in 1:length(clone_list)) {hamdist_subset[1,n] = length(clone_list[[n]])}
       out_rarecurve = hamdist_subset[1,,drop=FALSE]
       
       return(out_rarecurve)
     }
     
     
## PIPELINE & FILTERING ####
     ## initial filter with vcftools ####
     
     ## Landoltia
     ## put together VCFtools command
     landoltia_vcftools_cmd = paste("vcftools",
                                    "--vcf", shQuote("/mnt/c/WSL_storage/landoltia_populations.snps.vcf"),
                                    "--recode",
                                    "--recode-INFO-all",
                                    "--mac 3",
                                    "--minDP 5",
                                    "--min-alleles 2",
                                    "--max-alleles 2",
                                    "--max-missing 0.5",
                                    "--out", shQuote("landoltia_post_vcftools"))
          
     # call command
     system2("wsl", args = landoltia_vcftools_cmd)
     
     ## Lemna
     ## put together VCFtools command
     lemna_vcftools_cmd = paste("vcftools",
                                "--vcf", shQuote("/mnt/c/WSL_storage/lemna_populations.snps.vcf"),
                                "--recode",
                                "--recode-INFO-all",
                                "--mac 3",
                                "--minDP 5",
                                "--min-alleles 2",
                                "--max-alleles 2",
                                "--max-missing 0.5",
                                "--out", shQuote("lemna_post_vcftools"))
                                    
     
     # call command
     system2("wsl", args = lemna_vcftools_cmd)

     ## comprehensive filter in R ####     

     ## count individuals before cleaning
     # library(openxlsx)
     # landoltia_raw = read.vcfR("C:/WSL_storage/landoltia_populations.snps.vcf")
     # lemna_raw = read.vcfR("C:/WSL_storage/lemna_populations.snps.vcf")
     # 
     # raw_collection = read.xlsx("C:/Users/timte/Desktop/Brisbane/Chapter 1/Duckweed collection/Duckweed collection 03.09.2025.xlsx")
     # landoltia_collection_raw = raw_collection[raw_collection[,1] %in% colnames(landoltia_raw@gt),]
     # lemna_collection_raw = raw_collection[raw_collection[,1] %in% colnames(lemna_raw@gt),]
     # 
     # test = raw_collection[!raw_collection[,1] %in% c(landoltia_collection_raw[,1], lemna_collection_raw[,1]),]
     # test[,1]
     
     ## load SNP data and population maps
     landoltia_postvcftools = read.vcfR("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/landoltia_post_vcftools.recode.vcf")
     lemna_postvcftools = read.vcfR("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/lemna_post_vcftools.recode.vcf")
     
     landoltia_popmap = read.table("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/landoltia_popmap.txt", sep="\t")
     names(landoltia_popmap) = c("id", "pop")
     lemna_popmap = read.table("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/lemna_popmap.txt", sep="\t")
     names(lemna_popmap) = c("id", "pop")
     
     ## filter out everything below 5 depth and under 29 quality
     landoltia_postcvftools_qlty = hard_filter(landoltia_postvcftools, depth = 5, gq = 30)
     lemna_postcvftools_qlty = hard_filter(lemna_postvcftools, depth = 5, gq = 30)
     
     ## filter for allele balance between 0.1 and 0.9
     landoltia_postcvftools_qlty_ab = filter_allele_balance(landoltia_postcvftools_qlty, min.ratio = 0.1, max.ratio = 0.9)
     lemna_postcvftools_qlty_ab = filter_allele_balance(lemna_postcvftools_qlty, min.ratio = 0.1, max.ratio = 0.9)
     
     ## filter for max depth
     landoltia_postcvftools_qlty_ab_maxdp = max_depth(landoltia_postcvftools_qlty_ab, maxdepth = 100)
     lemna_postcvftools_qlty_ab_maxdp = max_depth(lemna_postcvftools_qlty_ab, maxdepth = 100)
     
     ##remove invariant SNPs generated during the genotype filtering steps
     landoltia_postcvftools_qlty_ab_maxdp = min_mac(landoltia_postcvftools_qlty_ab_maxdp, min.mac = 1)
     lemna_postcvftools_qlty_ab_maxdp = min_mac(lemna_postcvftools_qlty_ab_maxdp, min.mac = 1)
     
     ## remove samples with more than 60% missing SNPs
     landoltia_postcvftools_qlty_ab_maxdp_smpl06 = missing_by_sample(landoltia_postcvftools_qlty_ab_maxdp, cutoff = 0.6)
     lemna_postcvftools_qlty_ab_maxdp_smpl06 = missing_by_sample(lemna_postcvftools_qlty_ab_maxdp, cutoff = 0.6)
     
     ## subset population map to match retained samples
     landoltia_popmap = landoltia_popmap[landoltia_popmap$id %in% colnames(landoltia_postcvftools_qlty_ab_maxdp_smpl06@gt),]
     lemna_popmap = lemna_popmap[lemna_popmap$id %in% colnames(lemna_postcvftools_qlty_ab_maxdp_smpl06@gt),]
     
     #remove invariant sites generated by dropping individuals
     landoltia_postcvftools_qlty_ab_maxdp_smpl06 = min_mac(landoltia_postcvftools_qlty_ab_maxdp_smpl06, min.mac = 1)
     lemna_postcvftools_qlty_ab_maxdp_smpl06 = min_mac(lemna_postcvftools_qlty_ab_maxdp_smpl06, min.mac = 1)
     
     # ## Landoltia  - look at the effect of missing data on PCA
     #      ## set missing data thresholds
     #      cust_thresholds = c(0.9,0.92,0.94,0.96,0.98)
     #      
     #      ## apply missing data thresholds and compute PCAs
     #      landoltia_pca_plot_miss = assess_missing_data_pca(landoltia_postcvftools_qlty_ab_maxdp_smpl06,
     #                                                        thresholds = cust_thresholds, popmap = landoltia_popmap)
     #      
     #      ## add color to dataframe
     #      for (n in 1:length(landoltia_pca_plot_miss)) {
     #      
     #           ## make sinister palette
     #           black_to_red = colorRampPalette(c("black", "red"))
     #           
     #           ## make and add vector
     #           landoltia_pca_plot_miss[[n]]$col_vec = ifelse(landoltia_pca_plot_miss[[n]][,"missing"] < 0.05, black_to_red(7)[1],
     #                                                         ifelse(landoltia_pca_plot_miss[[n]][,"missing"] < 0.1, black_to_red(7)[2],
     #                                                                ifelse(landoltia_pca_plot_miss[[n]][,"missing"] < 0.15, black_to_red(7)[3],
     #                                                                       ifelse(landoltia_pca_plot_miss[[n]][,"missing"] < 0.20, black_to_red(7)[4],
     #                                                                              ifelse(landoltia_pca_plot_miss[[n]][,"missing"] < 0.25, black_to_red(7)[5],
     #                                                                                     ifelse(landoltia_pca_plot_miss[[n]][,"missing"] < 0.30, black_to_red(7)[6], black_to_red(7)[7]))))))
     #           
     #      }
     #      
     #      ## plot PCAs
     #      par(mfrow=c(2,3))
     #      for (n in 1:length(landoltia_pca_plot_miss)) {
     #           
     #           ## plot PCA results with missingness
     #           plot(landoltia_pca_plot_miss[[n]][,"PC1"],landoltia_pca_plot_miss[[n]][,"PC2"], 
     #                pch=19, col=landoltia_pca_plot_miss[[n]][,"col_vec"], main=paste(cust_thresholds[n]," % SNP missing", sep=""))
     #           legend("topright", legend = c("0.05 missing", "0.1 miss", "0.15 miss", "0.2 miss", "0.25 miss", "0.3 miss", "> 0.3 miss" ),
     #                  pch=19, col=black_to_red(7), cex=0.5)
     #      }
     # 
     # ## Lemna  - look at the effect of missing data on PCA
     #      ## set missing data thresholds
     #      cust_thresholds = c(0.9,0.92,0.94,0.96,0.98)
     #      
     #      ## apply missing data thresholds and compute PCAs
     #      lemna_pca_plot_miss = assess_missing_data_pca(lemna_postcvftools_qlty_ab_maxdp_smpl06,
     #                                                        thresholds = cust_thresholds, popmap = lemna_popmap)
     #      
     #      ## add color to dataframe
     #      for (n in 1:length(lemna_pca_plot_miss)) {
     #           
     #           ## make sinister palette
     #           black_to_red = colorRampPalette(c("black", "red"))
     #           
     #           ## make and add vector
     #           lemna_pca_plot_miss[[n]]$col_vec = ifelse(lemna_pca_plot_miss[[n]][,"missing"] < 0.05, black_to_red(7)[1],
     #                                                         ifelse(lemna_pca_plot_miss[[n]][,"missing"] < 0.1, black_to_red(7)[2],
     #                                                                ifelse(lemna_pca_plot_miss[[n]][,"missing"] < 0.15, black_to_red(7)[3],
     #                                                                       ifelse(lemna_pca_plot_miss[[n]][,"missing"] < 0.20, black_to_red(7)[4],
     #                                                                              ifelse(lemna_pca_plot_miss[[n]][,"missing"] < 0.25, black_to_red(7)[5],
     #                                                                                     ifelse(lemna_pca_plot_miss[[n]][,"missing"] < 0.30, black_to_red(7)[6], black_to_red(7)[7]))))))
     #           
     #      }
     #      
     #      ## plot PCAs
     #      par(mfrow=c(2,3))
     #      for (n in 1:length(lemna_pca_plot_miss)) {
     #           
     #           ## plot PCA results with missingness
     #           plot(lemna_pca_plot_miss[[n]][,"PC1"],lemna_pca_plot_miss[[n]][,"PC2"], 
     #                pch=19, col=lemna_pca_plot_miss[[n]][,"col_vec"], main=paste(cust_thresholds[n]," % SNP missing", sep=""))
     #           legend("topright", legend = c("0.05 missing", "0.1 miss", "0.15 miss", "0.2 miss", "0.25 miss", "0.3 miss", "> 0.3 miss" ),
     #                  pch=19, col=black_to_red(7), cex=0.5)
     #      }     
      
      
     ## apply more stringent missing SNP filter
     landoltia_final = missing_by_snp(landoltia_postcvftools_qlty_ab_maxdp_smpl06, cutoff = .95)
     lemna_final = missing_by_snp(lemna_postcvftools_qlty_ab_maxdp_smpl06, cutoff = .95)
     
     ## export data in vcf format for spider to make structure files
     # write.vcf(landoltia_final, "landoltia_final.vcf.gz")
     # write.vcf(lemna_final, "lemna_final.vcf.gz")
     
## GENETIC DISTANCE AND CLONES ####
     ## LANDOLTIA compute hamming distance ####
     
     ## replace homozygote, heterozygote, homozygote with 0,1,2
     landoltia_final_genind = vcfR2genind(landoltia_final)
     landoltia_genotypes_double = t(landoltia_final_genind@tab)
     landoltia_genotypes = landoltia_genotypes_double[seq(1, nrow(landoltia_genotypes_double), by = 2), ]
     
     ## create matrix for manual hamming distance calculation
     landoltia_hamdist_raw = matrix(NA, nrow=ncol(landoltia_genotypes), ncol=ncol(landoltia_genotypes))
     colnames(landoltia_hamdist_raw) = colnames(landoltia_genotypes)
     rownames(landoltia_hamdist_raw) = colnames(landoltia_genotypes)
     
     ##  compute hamming distance
     for (n in 1:ncol(landoltia_genotypes)) {
          
          ## compute hamming distance
          runner_vec = apply(landoltia_genotypes, 2, function(x) sum(abs(landoltia_genotypes[,n] - x), na.rm=TRUE))
          
          ## store in matrix
          landoltia_hamdist_raw[,n] = runner_vec
          
     }
     
     ## hardcode meaningful order of samples
     landoltia_ordered_names = c("P1S2",
                                 "P2S2",
                                 "P4S2",
                                 "P10S1b","P10S3","P10S3b","P10S4","P10S4b","P10S6","P10S7b","P10S8b","P10S9b","P10S10b","P10S11b","P10S12b","P10S13b","P10S14b","P10S15","P10S16","P10S18","P10S22","P10S23","P10S26","P10S31",
                                 "P11S10","P11S11","P11S12",
                                 "P12S4","P12S5","P12S6",
                                 "P13S1","P13S2","P13S3",
                                 "P14S1","P14S7","P14S8","P14S9","P14S12","P14S13","P14S16","P14S17","P14S18","P14S22","P14S23","P14S24","P14S28","P14S29","P14S30","P14S37","P14S38","P14S39","P14S40","P14S42","P14S46","P14S47","P14S48",
                                 "P15S1","P15S2","P15S3",
                                 "P16S4","P16S5","P16S6",
                                 "P17S4","P17S6",
                                 "P18S4","P18S6",
                                 "P19S4","P19S5","P19S10","P19S11","P19S16","P19S17","P19S18","P19S22","P19S23","P19S24","P19S25","P19S26","P19S27","P19S34","P19S35","P19S36","P19S40","P19S41","P19S42","P19S43","P19S44","P19S45","P19S52","P19S53","P19S54",
                                 "P20S1","P20S2","P20S3",
                                 "P22S1","P22S2","P22S3",
                                 "P23S4","P23S5",
                                 "P24S1","P24S2","P24S3",
                                 "P25S1","P25S2","P25S3",
                                 "P26S1","P26S2","P26S3",
                                 "P27S1","P27S2","P27S3","P27S7","P27S8","P27S9","P27S13","P27S14","P27S15","P27S22","P27S23","P27S24","P27S26","P27S28","P27S30","P27S31","P27S32","P27S33","P27S34","P27S35","P27S36","P27S37","P27S38","P27S49","P27S51",
                                 "P28S1","P28S2","P28S3",
                                 "P32S1","P32S2","P32S3",
                                 "P34S1","P34S2","P34S3",
                                 "P36S1","P36S2","P36S3","P36S7","P36S8","P36S9","P36S16","P36S17","P36S18","P36S19","P36S20","P36S21","P36S25","P36S26","P36S27","P36S34","P36S35","P36S36","P36S40","P36S41")
     
     ## order hamming distance matrix
     landoltia_hamdist_raw = landoltia_hamdist_raw[landoltia_ordered_names, landoltia_ordered_names]
     
     ## scale by maximum distance
     landoltia_hamdist = round(landoltia_hamdist_raw/(nrow(landoltia_final@fix)*2),3)
     
     ## LANDOLTIA identify clones based on distribution ####
     
     ## visualize the distribution
     par(mfrow=c(1,2))
     hist(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)], breaks=100,
          main = "Landoltia", col="purple", xlab="")
     abline(v=0.02, col="red", lty=2)
     hist(sort(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)])[1:5000], breaks=100,
          main = "Landoltia zoom in", col="purple", xlab="")
     abline(v=0.02, col="red", lty=2)
     landoltia_clone_cutoff = 0.02
     
     ## find plants with 0.02 difference
     clonal_candidates = which(landoltia_hamdist <= landoltia_clone_cutoff, arr.ind = T)
     rownames(clonal_candidates) = NULL
     
     ## remove self-comparison
     clonal_candidates = clonal_candidates[which(clonal_candidates[,1] != clonal_candidates[,2]),]
     
     ## remove duplicates
     n=0; whiler = TRUE
     while (whiler == TRUE) {
          
          ## counter
          n=n+1
          
          ## check for duplicates
          if (n >= nrow(clonal_candidates)) {
               ## stop loop
               whiler = FALSE
          }
          else if (sum(clonal_candidates[,1] == clonal_candidates[n,2] & clonal_candidates[,2] == clonal_candidates[n,1])==1) {
               ## remove if it does exist
               clonal_candidates = clonal_candidates[-n,]
          }
     }
     
     ## make clonal groups
     clone_list = vector(mode='list'); lc = 0
     for (n in 1:length(sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2]))))))) {
          
          ## list counter
          runner_candidate = sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2])))))[n]
          
          ## start vector
          clone_vector = runner_candidate
          
          ## clone vector
          whiler = TRUE
          while(whiler==TRUE) {
               pre_clone_vector = clone_vector
               clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,1]))),]))
               clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,2]))),]))
               if (identical(pre_clone_vector,clone_vector)){whiler = FALSE}
          }
          
          ## add indices to list
          clone_list[[n]] = sort(clone_vector)
          
     }
     
     ## remove duplicates
     clone_str = sapply(clone_list, paste, collapse = "-")
     clone_list = clone_list[!duplicated(clone_str)]
     
     ## extract sample names
     landoltia_clone_list = lapply(clone_list, function(x) rownames(landoltia_hamdist)[x])
     lengths = unlist(lapply(landoltia_clone_list, length))
     
     ## save in dataframe
     landoltia_clone_df = as.data.frame(matrix(NA, ncol=length(landoltia_clone_list), nrow=max(lengths)))
     landoltia_clone_names = vector()
     for (n in 1:length(landoltia_clone_list)) {landoltia_clone_names[n] = paste0("clone",n, sep="")}
     names(landoltia_clone_df) = landoltia_clone_names
     for (n in 1:length(landoltia_clone_list)) {
          
          ## transform to vector and put in matrix
          runner = landoltia_clone_list[[n]]
          runner = c(runner, rep(NA, max(lengths)-length(runner)))
          landoltia_clone_df[,n] = runner
          
     }
     
     ## export table
     # write.csv(landoltia_clone_df, paste("Landoltia_clones_hamming0-",landoltia_clone_cutoff,".csv", sep=""))
     
     ## LANDOLTIA create clone, no clone and no clone + alive hamming distance ####
     
     ## remove clones
     landoltia_noclone = landoltia_final
     for (n in 1:ncol(landoltia_clone_df)) {
          
          ## run through landoltia_clone_df by col
          runner_clone = landoltia_clone_df[,n]; runner_clone = runner_clone[!is.na(runner_clone)]
          
          # identify and remove clones
          to_keep = which(!colnames(landoltia_noclone@gt) %in% runner_clone[2:length(runner_clone)])
          landoltia_noclone = landoltia_noclone[, to_keep]
          
     }
     landoltia_noclone_names = colnames(landoltia_noclone@gt)[-1]
     
     ## adjust hamdist
     landoltia_noclone_hamdist = landoltia_hamdist[which(names(landoltia_hamdist[1,]) %in% landoltia_noclone_names),
                                                   which(names(landoltia_hamdist[1,]) %in% landoltia_noclone_names)]     
     
     ## landoltia still in collection
     landoltia_alive = c("P2S2",
                         "P11S12",
                         "P12S4", "P12S5", "P12S6",
                         "P13S1", "P13S2",
                         "P14S7", "P14S8", "P14S9", "P14S12", "P14S13", "P14S16", "P14S18", "P14S22","P14S23", "P14S24", "P14S30", "P14S37", "P14S38", "P14S41", "P14S42", "P14S46", "P14S48", "P14S52", "P14S53",
                         "P15S3",
                         "P19S11", "P19S12", "P19S16", "P19S17", "P19S22", "P19S23", "P19S26", "P19S34", "P19S36", "P19S40", "P19S41", "P19S42", "P19S43", "P19S52", "P19S53", "P19S54",
                         "P23S4", "P23S5", "P23S6",
                         "P24S2", "P24S3",
                         "P25S1", "P25S3",
                         "P26S1", "P26S3",
                         "P27S3", "P27S7", "P27S13", "P27S26", "P27S28", "P27S29", "P27S30", "P27S35", "P27S51", 
                         "P28S3",
                         "P36S1", "P36S2", "P36S7", "P36S16", "P36S34", "P36S41")
     
     ## gendist matrix for living plants
     landoltia_noclone_alive_hamdist = landoltia_noclone_hamdist[which(names(landoltia_noclone_hamdist[1,]) %in% landoltia_alive),
                                                                 which(names(landoltia_noclone_hamdist[1,]) %in% landoltia_alive)]
     
     
     ## LEMNA compute hamming distance ####
     
     ## replace homozygote, heterozygote, homozygote with 0,1,2
     lemna_final_genind = vcfR2genind(lemna_final)
     lemna_genotypes_double = t(lemna_final_genind@tab)
     lemna_genotypes = lemna_genotypes_double[seq(1, nrow(lemna_genotypes_double), by = 2), ]
     
     ## create matrix for manual hamming distance calculation
     lemna_hamdist_raw = matrix(NA, nrow=ncol(lemna_genotypes), ncol=ncol(lemna_genotypes))
     colnames(lemna_hamdist_raw) = colnames(lemna_genotypes)
     rownames(lemna_hamdist_raw) = colnames(lemna_genotypes)
     
     ##  compute hamming distance
     for (n in 1:ncol(lemna_genotypes)) {
       
       ## compute hamming distance
       runner_vec = apply(lemna_genotypes, 2, function(x) sum(abs(lemna_genotypes[,n] - x), na.rm=TRUE))
       
       ## store in matrix
       lemna_hamdist_raw[,n] = runner_vec
       
     }
     
     ## hardcode meaningful order of samples
     lemna_ordered_names = c("P5S1",
                             "P6S1","P6S3",
                             "P7S1","P7S2","P7S3",
                             "P10S1","P10S2","P10S2b","P10S5","P10S5b","P10S6b","P10S7","P10S8","P10S9","P10S10","P10S11","P10S12","P10S13","P10S14","P10S15b","P10S19","P10S20","P10S21","P10S28","P10S30",
                             "P11S7","P11S8","P11S9",
                             "P14S10","P14S11","P14S15","P14S19","P14S20","P14S21","P14S25","P14S26","P14S27","P14S31","P14S32","P14S34","P14S35","P14S36","P14S43","P14S44","P14S45","P14S49","P14S50","P14S55","P14S56","P14S57",
                             "P16S2","P16S3",
                             "P17S1","P17S2","P17S3",
                             "P18S2",
                             "P19S1","P19S2","P19S3","P19S7","P19S8","P19S9","P19S13","P19S14","P19S15","P19S19","P19S20","P19S21","P19S28","P19S29","P19S30","P19S31","P19S32","P19S33","P19S37","P19S38","P19S39","P19S46","P19S47","P19S48","P19S49","P19S50","P19S51",
                             "P21S1","P21S2",
                             "P22S4","P22S6",
                             "P23S3",
                             "P27S4","P27S5","P27S10","P27S11","P27S12","P27S16","P27S17","P27S18","P27S19","P27S20","P27S21","P27S25","P27S27","P27S39","P27S40","P27S41","P27S42","P27S43","P27S44","P27S45","P27S46","P27S47","P27S48","P27S52",
                             "P28S4","P28S5","P28S6",
                             "P30S1","P30S2","P30S3",
                             "P31S1","P31S2","P31S3",
                             "P32S4","P32S5","P32S6",
                             "P33S1","P33S2","P33S3",
                             "P34S5","P34S6",
                             "P35S1","P35S2","P35S3",
                             "P36S4","P36S5","P36S6","P36S10","P36S11","P36S12","P36S13","P36S14","P36S15","P36S22","P36S23","P36S24","P36S28","P36S29","P36S30","P36S31","P36S32","P36S33","P36S37","P36S38","P36S39")
     
     ## order hamming distance matrix
     lemna_hamdist_raw = lemna_hamdist_raw[lemna_ordered_names, lemna_ordered_names]
     
     ## scale by maximum distance
     lemna_hamdist = round(lemna_hamdist_raw/(nrow(lemna_final@fix)*2),3)
     
     ## LEMNA identify clones based on distribution ####
     
     ## visualize the distribution
     par(mfrow=c(1,2))
     hist(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)], breaks=100,
          main = "Lemna", col="darkgreen", xlab="")
     abline(v=0.02, col="red", lty=2)
     hist(sort(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)])[1:5000], breaks=100,
          main = "Lemna zoom in", col="darkgreen", xlab="")
     abline(v=0.02, col="red", lty=2)
     lemna_clone_cutoff = 0.02
     
     ## find plants with 0.02 difference
     clonal_candidates = which(lemna_hamdist <= lemna_clone_cutoff, arr.ind = T)
     rownames(clonal_candidates) = NULL
     
     ## remove self-comparison
     clonal_candidates = clonal_candidates[which(clonal_candidates[,1] != clonal_candidates[,2]),]
     
     ## remove duplicates
     n=0; whiler = TRUE
     while (whiler == TRUE) {
       
       ## counter
       n=n+1
       
       ## check for duplicates
       if (n >= nrow(clonal_candidates)) {
         ## stop loop
         whiler = FALSE
       }
       else if (sum(clonal_candidates[,1] == clonal_candidates[n,2] & clonal_candidates[,2] == clonal_candidates[n,1])==1) {
         ## remove if it does exist
         clonal_candidates = clonal_candidates[-n,]
       }
     }
     
     ## make clonal groups
     clone_list = vector(mode='list'); lc = 0
     for (n in 1:length(sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2]))))))) {
       
       ## list counter
       runner_candidate = sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2])))))[n]
       
       ## start vector
       clone_vector = runner_candidate
       
       ## clone vector
       whiler = TRUE
       while(whiler==TRUE) {
         pre_clone_vector = clone_vector
         clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,1]))),]))
         clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,2]))),]))
         if (identical(pre_clone_vector,clone_vector)){whiler = FALSE}
       }
       
       ## add indices to list
       clone_list[[n]] = sort(clone_vector)
       
     }
     
     ## remove duplicates
     clone_str = sapply(clone_list, paste, collapse = "-")
     clone_list = clone_list[!duplicated(clone_str)]
     
     ## extract sample names
     lemna_clone_list = lapply(clone_list, function(x) rownames(lemna_hamdist)[x])
     lengths = unlist(lapply(lemna_clone_list, length))
     
     ## save in dataframe
     lemna_clone_df = as.data.frame(matrix(NA, ncol=length(lemna_clone_list), nrow=max(lengths)))
     lemna_clone_names = vector()
     for (n in 1:length(lemna_clone_list)) {lemna_clone_names[n] = paste0("clone",n, sep="")}
     names(lemna_clone_df) = lemna_clone_names
     for (n in 1:length(lemna_clone_list)) {
       
       ## transform to vector and put in matrix
       runner = lemna_clone_list[[n]]
       runner = c(runner, rep(NA, max(lengths)-length(runner)))
       lemna_clone_df[,n] = runner
       
     }
     
     ## export table
     # write.csv(lemna_clone_df, paste("lemna_clones_hamming0-",lemna_clone_cutoff,".csv", sep=""))
     
     ## LEMNA create clone, no clone and no clone + alive hamming distance ####
     
     ## remove clones
     lemna_noclone = lemna_final
     for (n in 1:ncol(lemna_clone_df)) {
       
       ## run through lemna_clone_df by col
       runner_clone = lemna_clone_df[,n]; runner_clone = runner_clone[!is.na(runner_clone)]
       
       # identify and remove clones
       to_keep = which(!colnames(lemna_noclone@gt) %in% runner_clone[2:length(runner_clone)])
       lemna_noclone = lemna_noclone[, to_keep]
       
     }
     lemna_noclone_names = colnames(lemna_noclone@gt)[-1]
     
     ## adjust hamdist
     lemna_noclone_hamdist = lemna_hamdist[which(names(lemna_hamdist[1,]) %in% lemna_noclone_names),
                                                   which(names(lemna_hamdist[1,]) %in% lemna_noclone_names)]     
     
     ## lemna still in collection
     lemna_alive = c("P7S2", "P7S3",
                     "P10S28",
                     "P11S8",
                     "P14S25", "P14S26", "P14S43", "P14S45",
                     "P19S19", "P19S20", "P19S29", "P19S30", "P19S31", "P19S38", "P19S39", "P19S49",
                     "P23S1", "P23S2", "P23S3", "P23S8",
                     "P27S4", "P27S10", "P27S42","P27S52",
                     "P29S1",
                     "P31S3",
                     "P35S2", "P35S3")
     
     ## gendist matrix for living plants
     lemna_noclone_alive_hamdist = lemna_noclone_hamdist[which(names(lemna_noclone_hamdist[1,]) %in% lemna_alive),
                                                                 which(names(lemna_noclone_hamdist[1,]) %in% lemna_alive)]
     
     
     
     
     
     
     
## DATA VISUALISATION / RESULTS ####
     ## FIGURE 1: Kinship Heatmaps ####
     
     ## LANDOLTIA
     
       ## set parameters
       col_pal = colorRampPalette(c("purple", "white"))
       colbreaks = 10
       clone_cutoff = 0.02
           
       ## set up split screen layout
       par(mar=c(2,2,2,2),oma=c(3,3,3,3))
       split.screen(rbind(c(0, 0.8, 0, 1), 
                          c(0.8, 1, 0, 1)))
       
       ## create heatmap
       screen(1)
       image(1:ncol(landoltia_hamdist), 1:nrow(landoltia_hamdist),
             landoltia_hamdist, axes = FALSE, frame=FALSE, xlab="", ylab="",
             col=c("black",col_pal(colbreaks-2)),
             breaks=c(0,seq(from=clone_cutoff, to=max(c(landoltia_hamdist,lemna_hamdist)), length.out = colbreaks-1)))
       
       
       ## generate axis positions and labels for waterbody
       landoltia_waterbody = unique(substr(colnames(landoltia_hamdist),1,3))
       landoltia_waterbody_pos = vector()
       for(n in 1:length(landoltia_waterbody)) {landoltia_waterbody_pos[n] = median(which(substr(colnames(landoltia_hamdist),1,3) == landoltia_waterbody[n]))}
       mtext(landoltia_waterbody, side = 1, at = landoltia_waterbody_pos, cex = 0.5, line=0.5 ,las=2)
       mtext(landoltia_waterbody, side = 2, at = landoltia_waterbody_pos, cex = 0.5, line=0.5 ,las=2)
       
       ## generate separating lines
       landoltia_deep_waterbody = c("P10", "P14", "P19", "P27", "P36")
       landoltia_deep_waterbody_pos = list()
       for(n in 1:length(landoltia_deep_waterbody)) {
         
         runner_min = min(which(substr(colnames(landoltia_hamdist),1,3) == landoltia_deep_waterbody[n]))-0.5
         runner_max = max(which(substr(colnames(landoltia_hamdist),1,3) == landoltia_deep_waterbody[n]))+0.5
         
         landoltia_deep_waterbody_pos[[n]] = c(runner_min, runner_max) 
         
         }
       abline(v=(unlist(landoltia_deep_waterbody_pos)), col="white")  
       abline(h=(unlist(landoltia_deep_waterbody_pos)), col="white")  
       close.screen(1)
       
       ## add legend
       screen(2)
       par(mar=c(2,1,2,4))
       plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,colbreaks),xaxt="n",yaxt="n",bty="n")
       for (n in 1:(colbreaks-1)) {rect(1,n-1,2,n, col=c("black",col_pal(colbreaks-2))[n])}
       axis(side=4, at=c(1:(colbreaks-1))-0.5, 
            labels=round(seq(from=min(c(landoltia_hamdist, lemna_hamdist)), 
                             to=max(c(landoltia_hamdist, lemna_hamdist)), 
                             length.out=colbreaks-1),2), 
            las=2, lwd=0, lwd.ticks=1)
       close.screen(2)
       close.screen(all.screens=TRUE)
       
     ## LEMNA
       
       ## set parameters
       col_pal = colorRampPalette(c("darkgreen", "white"))
       colbreaks = 10
       clone_cutoff = 0.02
       
       ## set up split screen layout
       par(mar=c(2,2,2,2),oma=c(3,3,3,3))
       split.screen(rbind(c(0, 0.8, 0, 1), 
                          c(0.8, 1, 0, 1)))
       
       ## create heatmap
       screen(1)
       image(1:ncol(lemna_hamdist), 1:nrow(lemna_hamdist),
             lemna_hamdist, axes = FALSE, frame=FALSE, xlab="", ylab="",
             col=c("black",col_pal(colbreaks-2)),
             breaks=c(0,seq(from=clone_cutoff, to=max(c(landoltia_hamdist,lemna_hamdist)), length.out = colbreaks-1)))
       
       
       ## generate axis positions and labels for waterbody
       lemna_waterbody = unique(substr(colnames(lemna_hamdist),1,3))
       lemna_waterbody_pos = vector()
       for(n in 1:length(lemna_waterbody)) {lemna_waterbody_pos[n] = median(which(substr(colnames(lemna_hamdist),1,3) == lemna_waterbody[n]))}
       mtext(lemna_waterbody, side = 1, at = lemna_waterbody_pos, cex = 0.5, line=0.5 ,las=2)
       mtext(lemna_waterbody, side = 2, at = lemna_waterbody_pos, cex = 0.5, line=0.5 ,las=2)
       
       ## generate separating lines
       lemna_deep_waterbody = c("P10", "P14", "P19", "P27", "P36")
       lemna_deep_waterbody_pos = list()
       for(n in 1:length(lemna_deep_waterbody)) {
         
         runner_min = min(which(substr(colnames(lemna_hamdist),1,3) == lemna_deep_waterbody[n]))-0.5
         runner_max = max(which(substr(colnames(lemna_hamdist),1,3) == lemna_deep_waterbody[n]))+0.5
         
         lemna_deep_waterbody_pos[[n]] = c(runner_min, runner_max) 
         
       }
       abline(v=(unlist(lemna_deep_waterbody_pos)), col="white")  
       abline(h=(unlist(lemna_deep_waterbody_pos)), col="white")  
       close.screen(1)
       
       ## add legend
       screen(2)
       par(mar=c(2,1,2,4))
       plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,colbreaks),xaxt="n",yaxt="n",bty="n")
       for (n in 1:(colbreaks-1)) {rect(1,n-1,2,n, col=c("black",col_pal(colbreaks-2))[n])}
       axis(side=4, at=c(1:(colbreaks-1))-0.5, 
            labels=round(seq(from=min(c(landoltia_hamdist, lemna_hamdist)), 
                             to=max(c(landoltia_hamdist, lemna_hamdist)), 
                             length.out=colbreaks-1),2), 
            las=2, lwd=0, lwd.ticks=1)
       close.screen(2)
       close.screen(all.screens=TRUE)
       
     ## FIGURE 2: Distance Decay Plots ####
     
     ## read and transform coordinates
     all_coordinates = read.csv("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/duckweed_coordinates.csv")
     all_coordinates$latitude = sapply(all_coordinates[,"GPS_S"], convert_dmm_to_dd)
     all_coordinates$longitute = sapply(all_coordinates[,"GPS_E"], convert_dmm_to_dd)
     
     ## assemble plot
     
     par(mfrow=c(1,2))
     
     ## LANDOLTIA
     
       ## subset to landoltia
       landoltia_coor = all_coordinates[which(colnames(landoltia_hamdist) %in% all_coordinates[,"ID"]),]
       
       ## calculate geographic distance
       landoltia_geodist = round(distm(landoltia_coor[, c("longitute", "latitude")], fun = distVincentyEllipsoid),2)
       rownames(landoltia_geodist) = landoltia_coor[,"ID"]; colnames(landoltia_geodist) = landoltia_coor[,"ID"]
       
       ## reduce matrix to triangle
       landoltia_hamdist_vector = landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)]
       landoltia_geodist_vector = landoltia_geodist[lower.tri(landoltia_geodist, diag=FALSE)]
       
       # replace 0 values, for log transform
       landoltia_geodist_vector[which(landoltia_geodist_vector==0)] = 1
       
       ## assemble plot
       plot(log(landoltia_geodist_vector), landoltia_hamdist_vector,
            pch=21, bg="purple", cex=0.8, ylim=c(-0.01,max(c(lemna_hamdist,landoltia_hamdist))),
            xlab="Geographic distance log(m)", 
            ylab="Genetic distance")
       abline(v=log(10), lty=2); text(log(10), -0.01, label="10m")
       abline(v=log(100), lty=2); text(log(100), -0.01, label="100m")
       abline(v=log(1000), lty=2); text(log(1000), -0.01, label="1km")
       abline(v=log(10000), lty=2); text(log(10000), -0.01, label="10km")
       abline(v=log(100000), lty=2); text(log(100000), -0.01, label="100km")
       
     ## LEMNA
       
       ## subset to lemna
       lemna_coor = all_coordinates[which(colnames(lemna_hamdist) %in% all_coordinates[,"ID"]),]
       
       ## calculate geographic distance
       lemna_geodist = round(distm(lemna_coor[, c("longitute", "latitude")], fun = distVincentyEllipsoid),2)
       rownames(lemna_geodist) = lemna_coor[,"ID"]; colnames(lemna_geodist) = lemna_coor[,"ID"]
       
       ## reduce matrix to triangle
       lemna_hamdist_vector = lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)]
       lemna_geodist_vector = lemna_geodist[lower.tri(lemna_geodist, diag=FALSE)]
       
       # replace 0 values, for log transform
       lemna_geodist_vector[which(lemna_geodist_vector==0)] = 1
       
       ## assemble plot
       plot(log(lemna_geodist_vector), lemna_hamdist_vector,
            pch=21, bg="darkgreen", cex=0.8, ylim=c(-0.01,max(c(lemna_hamdist, landoltia_hamdist))),
            xlab="Geographic distance log(m)", 
            ylab="Genetic distance")
       abline(v=log(10), lty=2); text(log(10), -0.01, label="10m")
       abline(v=log(100), lty=2); text(log(100), -0.01, label="100m")
       abline(v=log(1000), lty=2); text(log(1000), -0.01, label="1km")
       abline(v=log(10000), lty=2); text(log(10000), -0.01, label="10km")
       abline(v=log(100000), lty=2); text(log(100000), -0.01, label="100km")

     
     ## Mantel tests
     
     ## LANDOLTIA
     vegan::mantel(landoltia_hamdist, landoltia_geodist, permutations = 1000)
     
     ## LEMNA
     vegan::mantel(lemna_hamdist, lemna_geodist, permutations = 1000)
     
       
     ## LANDOLTIA within vs across pond diversity ####
          
     ## subset the ponds   
     landoltia_p10 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P10"),which(substr(colnames(landoltia_hamdist),1,3) == "P10")]
     landoltia_p14 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P14"),which(substr(colnames(landoltia_hamdist),1,3) == "P14")]
     landoltia_p19 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P19"),which(substr(colnames(landoltia_hamdist),1,3) == "P19")]
     landoltia_p27 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P27"),which(substr(colnames(landoltia_hamdist),1,3) == "P27")]
     landoltia_p36 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P36"),which(substr(colnames(landoltia_hamdist),1,3) == "P36")]
     
     ## draw and fill plot
     plot(NULL, xlim=c(0.5,6.5), ylim=c(0,max(landoltia_hamdist)), xlab="", ylab="Hamming distance",
          main="Landoltia: within and across pond diversity", xaxt = "n")
     axis(1, at = 1:6, las=2, labels = c("Pond10", "Pond14", "Pond19", "Pond27", "Pond36", "all data"))
     
     ## clone cutoff
     abline(h=landoltia_clone_cutoff, lty=2)
     text(0.2, landoltia_clone_cutoff+8, label="clone cutoff", cex=0.8, pos=4)
     
     ## pond 10
     points(runif(length(landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)]), min = 0.8, max = 1.2),
            landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)], pch=19, col=scales::alpha("purple",0.1))
     points(1, mean(landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(1,1), 
           c(mean(landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)])+sd(landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)]),
             mean(landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)])-sd(landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)])),
           lwd=2)
     
     ## pond 14
     points(runif(length(landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)]), min = 1.8, max = 2.2),
            landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)], pch=19, col=scales::alpha("purple",0.1))
     points(2, mean(landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(2,2), 
           c(mean(landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)])+sd(landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)]),
             mean(landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)])-sd(landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)])),
           lwd=2)
     
     ## pond 19
     points(runif(length(landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)]), min = 2.8, max = 3.2),
            landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)], pch=19, col=scales::alpha("purple",0.1))
     points(3, mean(landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(3,3), 
           c(mean(landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)])+sd(landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)]),
             mean(landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)])-sd(landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)])),
           lwd=2)
     
     ## pond 27
     points(runif(length(landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)]), min = 3.8, max = 4.2),
            landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)], pch=19, col=scales::alpha("purple",0.1))
     points(4, mean(landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(4,4), 
           c(mean(landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)])+sd(landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)]),
             mean(landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)])-sd(landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)])),
           lwd=2)
     
     ## pond 36
     points(runif(length(landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)]), min = 4.8, max = 5.2),
            landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)], pch=19, col=scales::alpha("purple",0.1))
     points(5, mean(landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(5,5), 
           c(mean(landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)])+sd(landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)]),
             mean(landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)])-sd(landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)])),
           lwd=2)
     
     # ## all data
     # points(runif(length(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)]), min = 5.8, max = 6.2),
     #        landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)], pch=19, col=scales::alpha("purple",0.1))
     # points(6, mean(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)]), pch=21, bg="black", cex=2)
     # lines(c(6,6), 
     #       c(mean(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)])+sd(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)]),
     #         mean(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)])-sd(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)])),
     #       lwd=2)
     
     ## random permutation means
     perm_mean = vector()
     for (n in 1:10000) {
       
       runner_sample = sample(colnames(landoltia_hamdist), 22)
       runner_matrix = landoltia_hamdist[runner_sample,runner_sample]
       perm_mean[n] = mean(runner_matrix[lower.tri(runner_matrix, diag=FALSE)])
       
     }
     points(rep(6,10000), perm_mean)
     
     
     sd
          
     ## LEMNA within vs across pond diversity ####
     
     ## subset the ponds   
     lemna_p10 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P10"),which(substr(colnames(lemna_hamdist),1,3) == "P10")]
     lemna_p14 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P14"),which(substr(colnames(lemna_hamdist),1,3) == "P14")]
     lemna_p19 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P19"),which(substr(colnames(lemna_hamdist),1,3) == "P19")]
     lemna_p27 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P27"),which(substr(colnames(lemna_hamdist),1,3) == "P27")]
     lemna_p36 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P36"),which(substr(colnames(lemna_hamdist),1,3) == "P36")]
     
     ## draw and fill plot
     plot(NULL, xlim=c(0.5,6.5), ylim=c(0,max(lemna_hamdist)), xlab="", ylab="Hamming distance",
          main="Lemna: within and across pond diversity (OUTLIER NOT SHOWN)", xaxt = "n")
     axis(1, at = 1:6, las=2, labels = c("Pond10", "Pond14", "Pond19", "Pond27", "Pond36", "all data"))
     
     ## clone cutoff
     abline(h=lemna_clone_cutoff, lty=2)
     text(0.2, lemna_clone_cutoff+8, label="clone cutoff", cex=0.8, pos=4)
     
     ## pond 10
     points(runif(length(lemna_p10[lower.tri(lemna_p10, diag=FALSE)]), min = 0.8, max = 1.2),
            lemna_p10[lower.tri(lemna_p10, diag=FALSE)], pch=19, col=scales::alpha("darkgreen",0.1))
     points(1, mean(lemna_p10[lower.tri(lemna_p10, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(1,1), 
           c(mean(lemna_p10[lower.tri(lemna_p10, diag=FALSE)])+sd(lemna_p10[lower.tri(lemna_p10, diag=FALSE)]),
             mean(lemna_p10[lower.tri(lemna_p10, diag=FALSE)])-sd(lemna_p10[lower.tri(lemna_p10, diag=FALSE)])),
           lwd=2)
     
     ## pond 14
     points(runif(length(lemna_p14[lower.tri(lemna_p14, diag=FALSE)]), min = 1.8, max = 2.2),
            lemna_p14[lower.tri(lemna_p14, diag=FALSE)], pch=19, col=scales::alpha("darkgreen",0.1))
     points(2, mean(lemna_p14[lower.tri(lemna_p14, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(2,2), 
           c(mean(lemna_p14[lower.tri(lemna_p14, diag=FALSE)])+sd(lemna_p14[lower.tri(lemna_p14, diag=FALSE)]),
             mean(lemna_p14[lower.tri(lemna_p14, diag=FALSE)])-sd(lemna_p14[lower.tri(lemna_p14, diag=FALSE)])),
           lwd=2)
     
     ## pond 19
     points(runif(length(lemna_p19[lower.tri(lemna_p19, diag=FALSE)]), min = 2.8, max = 3.2),
            lemna_p19[lower.tri(lemna_p19, diag=FALSE)], pch=19, col=scales::alpha("darkgreen",0.1))
     points(3, mean(lemna_p19[lower.tri(lemna_p19, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(3,3), 
           c(mean(lemna_p19[lower.tri(lemna_p19, diag=FALSE)])+sd(lemna_p19[lower.tri(lemna_p19, diag=FALSE)]),
             mean(lemna_p19[lower.tri(lemna_p19, diag=FALSE)])-sd(lemna_p19[lower.tri(lemna_p19, diag=FALSE)])),
           lwd=2)
     
     ## pond 27
     points(runif(length(lemna_p27[lower.tri(lemna_p27, diag=FALSE)]), min = 3.8, max = 4.2),
            lemna_p27[lower.tri(lemna_p27, diag=FALSE)], pch=19, col=scales::alpha("darkgreen",0.1))
     points(4, mean(lemna_p27[lower.tri(lemna_p27, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(4,4), 
           c(mean(lemna_p27[lower.tri(lemna_p27, diag=FALSE)])+sd(lemna_p27[lower.tri(lemna_p27, diag=FALSE)]),
             mean(lemna_p27[lower.tri(lemna_p27, diag=FALSE)])-sd(lemna_p27[lower.tri(lemna_p27, diag=FALSE)])),
           lwd=2)
     
     ## pond 36
     points(runif(length(lemna_p36[lower.tri(lemna_p36, diag=FALSE)]), min = 4.8, max = 5.2),
            lemna_p36[lower.tri(lemna_p36, diag=FALSE)], pch=19, col=scales::alpha("darkgreen",0.1))
     points(5, mean(lemna_p36[lower.tri(lemna_p36, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(5,5), 
           c(mean(lemna_p36[lower.tri(lemna_p36, diag=FALSE)])+sd(lemna_p36[lower.tri(lemna_p36, diag=FALSE)]),
             mean(lemna_p36[lower.tri(lemna_p36, diag=FALSE)])-sd(lemna_p36[lower.tri(lemna_p36, diag=FALSE)])),
           lwd=2)
     
     ## all data
     points(runif(length(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)]), min = 5.8, max = 6.2),
            lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)], pch=19, col=scales::alpha("darkgreen",0.1))
     points(6, mean(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)]), pch=21, bg="black", cex=2)
     lines(c(6,6), 
           c(mean(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)])+sd(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)]),
             mean(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)])-sd(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)])),
           lwd=2)
     
     ## landoltia plot PCA ####

     ## transform file
     landoltia_final_pca = vcfR2genlight(landoltia_final)
     
     par(mfrow=c(1,3))
     ## run and plot PCA
     no_missing_data = tab(landoltia_final_pca, NA.method = "mean")
     
     dudipca_result = dudi.pca(no_missing_data, scannf = FALSE, nf = 3)
     plot(dudipca_result$li[,1], dudipca_result$li[,2], 
          xlab = "PC1", ylab = "PC2",
          main = "Landoltia ade4::dudi.pca()", pch = 19, col = "purple")
     text(dudipca_result$li[,1], dudipca_result$li[,2], 
          labels = rownames(dudipca_result$li), pos = 3, cex = 0.4, col="gray50")
     
     glpca_result = glPca(landoltia_final_pca, nf=3)
     plot(glpca_result$scores[,"PC1"], glpca_result$scores[,"PC2"],
          xlab = "PC1", ylab = "PC2",
          main = "Landoltia adegenet::glPca", pch = 19, col = "purple")
     text(glpca_result$scores[,"PC1"], glpca_result$scores[,"PC2"],
          labels = rownames(glpca_result$scores), pos = 3, cex = 0.4, col="gray50")
     
     lea_pca = LEA::pca(landoltia_final_geno)
     plot(lea_pca$projections, main="landoltia LEA::pca", pch=19, col="purple")
     
     ## lemna plot PCA ####

     ## transform file
     lemna_final_pca = vcfR2genlight(lemna_final)
     
     par(mfrow=c(1,2))
     ## run and plot PCA
     no_missin_data = tab(lemna_final_pca, NA.method = "mean")
     dudipca_result = dudi.pca(no_missin_data, scannf = FALSE, nf = 3)
     plot(dudipca_result$li[,1], dudipca_result$li[,2], 
          xlab = "PC1", ylab = "PC2",
          main = "Lemna ade4::dudi.pca()", pch = 19, col = "darkgreen")
     text(dudipca_result$li[,1], dudipca_result$li[,2], 
          labels = rownames(dudipca_result$li), pos = 3, cex = 0.4, col="gray50")
     
     glpca_result = glPca(lemna_final_pca, nf=3)
     plot(glpca_result$scores[,"PC1"], glpca_result$scores[,"PC2"],
          xlab = "PC1", ylab = "PC2",
          main = "Lemna adegenet::glPca", pch = 19, col = "darkgreen")
     text(glpca_result$scores[,"PC1"], glpca_result$scores[,"PC2"],
          labels = rownames(glpca_result$scores), pos = 3, cex = 0.4, col="gray50")
     
     
     ## landoltia population structure with LEA ####
     
     ## transform data
     landoltia_final_geno = vcf2geno("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/landoltia_final.vcf")
     
     ## run analysis
     landoltia_project = NULL
     landoltia_project = snmf(landoltia_final_geno,
                              K = 1:10,
                              entropy = TRUE,
                              repetitions = 20,
                              project = "new")
                   
     
     ## plot entropy to identify best number of clusters
     plot(landoltia_project, main="Landoltia")
     
     ## STRUCTURE plot
     best = which.min(cross.entropy(landoltia_project, K = 9))
     my.colors <- rainbow(9)
     barchart(landoltia_project, K = 9, run = best,
              border = NA, space = 0,
              col = my.colors,
              xlab = "Individuals",
              ylab = "Ancestry proportions",
              main = "Ancestry matrix") -> bp
     axis(1, at = 1:length(bp$order),
          labels = bp$order, las=1,
          cex.axis = .4)
     
     
     ## lemna population structure with LEA ####
     
     ## transform data
     lemna_final_geno = vcf2geno("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/lemna_final.vcf")
     
     ## run analysis
     lemna_project = NULL
     lemna_project = snmf(lemna_final_geno,
                    K = 1:10,
                    entropy = TRUE,
                    repetitions = 10,
                    project = "new")
     
     ## plot entropy to identify best number of clusters
     plot(lemna_project, main="Lemna")
     
     ## STRUCTURE plot
     best = which.min(cross.entropy(lemna_project, K = 7))
     my.colors <- brewer.pal(7, "Set3")
     barchart(lemna_project, K = 7, run = best,
              border = NA, space = 0,
              col = my.colors,
              xlab = "Individuals",
              ylab = "Ancestry proportions",
              main = "Ancestry matrix") -> bp
     axis(1, at = 1:length(bp$order),
          labels = bp$order, las=1,
          cex.axis = .4)
     
     ## trees ####
     
     tree = bionj(landoltia_hamdist)
     tree$edge.length[tree$edge.length < 0] <- 0
     chronotree = chronos(tree,control = chronos.control(iter.max = 1000), model="relaxed")
     plot.phylo(chronotree, cex=0.4, type = "fan", tip.color=ifelse(grepl("P10", tree$tip.label), "red", 
                                                                    ifelse(grepl("P14", tree$tip.label), "blue", 
                                                                           ifelse(grepl("P19", tree$tip.label), "darkgreen",
                                                                                  ifelse(grepl("P27", tree$tip.label), "purple",
                                                                                         ifelse(grepl("P36", tree$tip.label), "orange", "black")))))
                
     )
     
     tree = bionj(lemna_hamdist)
     tree$edge.length[tree$edge.length < 0] <- 0
     chronotree = chronos(tree,control = chronos.control(iter.max = 1000), model="relaxed")
     plot.phylo(chronotree, cex=0.4, type = "fan", tip.color=ifelse(grepl("P10", tree$tip.label), "red", 
                                                                    ifelse(grepl("P14", tree$tip.label), "blue", 
                                                                           ifelse(grepl("P19", tree$tip.label), "darkgreen",
                                                                                  ifelse(grepl("P27", tree$tip.label), "purple",
                                                                                         ifelse(grepl("P36", tree$tip.label), "orange", "black")))))
                
     )
       
     ## rarefaction curves ####
     
     landoltia_pond10 = landoltia_hamdist[c(which(substr(colnames(landoltia_hamdist),1,3) == "P10")),c(which(substr(colnames(landoltia_hamdist),1,3) == "P10"))]      
     landoltia_pond14 = landoltia_hamdist[c(which(substr(colnames(landoltia_hamdist),1,3) == "P14")),c(which(substr(colnames(landoltia_hamdist),1,3) == "P14"))]      
     landoltia_pond19 = landoltia_hamdist[c(which(substr(colnames(landoltia_hamdist),1,3) == "P19")),c(which(substr(colnames(landoltia_hamdist),1,3) == "P19"))]      
     landoltia_pond27 = landoltia_hamdist[c(which(substr(colnames(landoltia_hamdist),1,3) == "P27")),c(which(substr(colnames(landoltia_hamdist),1,3) == "P27"))]      
     landoltia_pond36 = landoltia_hamdist[c(which(substr(colnames(landoltia_hamdist),1,3) == "P36")),c(which(substr(colnames(landoltia_hamdist),1,3) == "P36"))]      
     
     landoltia_rarecurve_all = hamdist_to_rarecurve(landoltia_hamdist)
     landoltia_rarecurve_P10 = hamdist_to_rarecurve(landoltia_pond10)
     landoltia_rarecurve_P14 = hamdist_to_rarecurve(landoltia_pond14)
     landoltia_rarecurve_P19 = hamdist_to_rarecurve(landoltia_pond19)
     landoltia_rarecurve_P27 = hamdist_to_rarecurve(landoltia_pond27)
     landoltia_rarecurve_P36 = hamdist_to_rarecurve(landoltia_pond36)
     
     landoltia_rareplot_all = rarecurve(landoltia_rarecurve_all)
     landoltia_rareplot_P10 = rarecurve(landoltia_rarecurve_P10)
     landoltia_rareplot_P14 = rarecurve(landoltia_rarecurve_P14)
     landoltia_rareplot_P19 = rarecurve(landoltia_rarecurve_P19)
     landoltia_rareplot_P27 = rarecurve(landoltia_rarecurve_P27)
     landoltia_rareplot_P36 = rarecurve(landoltia_rarecurve_P36)
     
     plot(1:length(landoltia_rareplot_all[[1]]),landoltia_rareplot_all[[1]], type="l")
     lines(1:length(landoltia_rareplot_P10[[1]]),landoltia_rareplot_P10[[1]])
     lines(1:length(landoltia_rareplot_P14[[1]]),landoltia_rareplot_P14[[1]])
     lines(1:length(landoltia_rareplot_P19[[1]]),landoltia_rareplot_P19[[1]])
     lines(1:length(landoltia_rareplot_P27[[1]]),landoltia_rareplot_P27[[1]])
     lines(1:length(landoltia_rareplot_P36[[1]]),landoltia_rareplot_P36[[1]])
     
     
     lemna_pond10 = lemna_hamdist[c(which(substr(colnames(lemna_hamdist),1,3) == "P10")),c(which(substr(colnames(lemna_hamdist),1,3) == "P10"))]      
     lemna_pond14 = lemna_hamdist[c(which(substr(colnames(lemna_hamdist),1,3) == "P14")),c(which(substr(colnames(lemna_hamdist),1,3) == "P14"))]      
     lemna_pond19 = lemna_hamdist[c(which(substr(colnames(lemna_hamdist),1,3) == "P19")),c(which(substr(colnames(lemna_hamdist),1,3) == "P19"))]      
     lemna_pond27 = lemna_hamdist[c(which(substr(colnames(lemna_hamdist),1,3) == "P27")),c(which(substr(colnames(lemna_hamdist),1,3) == "P27"))]      
     lemna_pond36 = lemna_hamdist[c(which(substr(colnames(lemna_hamdist),1,3) == "P36")),c(which(substr(colnames(lemna_hamdist),1,3) == "P36"))]      
    
     lemna_rarecurve_all = hamdist_to_rarecurve(lemna_hamdist)
     lemna_rarecurve_P10 = hamdist_to_rarecurve(lemna_pond10)
     lemna_rarecurve_P14 = hamdist_to_rarecurve(lemna_pond14)
     lemna_rarecurve_P19 = hamdist_to_rarecurve(lemna_pond19)
     lemna_rarecurve_P27 = hamdist_to_rarecurve(lemna_pond27)
     lemna_rarecurve_P36 = hamdist_to_rarecurve(lemna_pond36)
     
     lemna_rareplot_all = rarecurve(lemna_rarecurve_all)
     lemna_rareplot_P10 = rarecurve(lemna_rarecurve_P10)
     lemna_rareplot_P14 = rarecurve(lemna_rarecurve_P14)
     lemna_rareplot_P19 = rarecurve(lemna_rarecurve_P19)
     lemna_rareplot_P27 = rarecurve(lemna_rarecurve_P27)
     lemna_rareplot_P36 = rarecurve(lemna_rarecurve_P36)
     
     plot(1:length(lemna_rareplot_all[[1]]),lemna_rareplot_all[[1]], type="l")
     plot(1:length(lemna_rareplot_P10[[1]]),lemna_rarecurve_all_P10[[1]], type="l", ylim=c(0,5), xlim=c(0,22))
     lines(1:length(lemna_rareplot_P10[[1]]),lemna_rareplot_P10[[1]])
     lines(1:length(lemna_rareplot_P14[[1]]),lemna_rareplot_P14[[1]])
     lines(1:length(lemna_rareplot_P19[[1]]),lemna_rareplot_P19[[1]])
     lines(1:length(lemna_rareplot_P27[[1]]),lemna_rareplot_P27[[1]])
     lines(1:length(lemna_rareplot_P36[[1]]),lemna_rareplot_P36[[1]])
     
     
     
     sd
     ## mantel test and the like ####
     
     library(ecodist)
     MRM(as.dist(landoltia_hamdist) ~ as.dist(landoltia_geodist), nperm = 999)
     
     mantel(landoltia_hamdist,landoltia_geodist)
     
## scraps ####
     ## landoltia compute geographical distance heatmap ####
     
     
     #add ordered names
     landoltia_geodist_matrix = landoltia_geodist_matrix[landoltia_ordered_names, landoltia_ordered_names]
     
     ## create heatmap
     landoltia_geodist_matrix = as.matrix(landoltia_geodist_matrix)
     heatmap(landoltia_geodist_matrix,
             main = "Geographic Distance Heatmap", Colv=NA, Rowv=NA,
             col = colorRampPalette(c("firebrick1", "dodgerblue"))(100))
     
     ## lemna compute geographical distance heatmap ####
     
     ## read file
     all_coordinates = read.csv("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/duckweed_coordinates.csv")
     
     ## only keep the samples with genetic data
     lemna_coordinates = all_coordinates[which(all_coordinates[,"ID"] %in% dimnames(lemna_gendist_matrix)[[2]]),]; rownames(lemna_coordinates) <- NULL
     
     ## convert my GPS to something one can caluculate distances with
     convert_dmm_to_dd = function(coord_str) {
          parts = strsplit(coord_str, "°")[[1]]
          return(as.numeric(parts[1]) + as.numeric(parts[2]) / 60)
     }
     lemna_coordinates$latitude = sapply(lemna_coordinates[,"GPS_S"], convert_dmm_to_dd)
     lemna_coordinates$longitute = sapply(lemna_coordinates[,"GPS_E"], convert_dmm_to_dd)
     rownames(lemna_coordinates) = lemna_coordinates[,"ID"]
     
     ## calculate geographic distance
     lemna_geodist_matrix = distm(lemna_coordinates[, c("longitute", "latitude")], fun = distHaversine)
     rownames(lemna_geodist_matrix) = lemna_coordinates[,"ID"]; colnames(lemna_geodist_matrix) = lemna_coordinates[,"ID"]
     
     #add ordered names
     lemna_geodist_matrix = lemna_geodist_matrix[lemna_ordered_names, lemna_ordered_names]
     
     ## create heatmap
     lemna_geodist_matrix = as.matrix(lemna_geodist_matrix)
     heatmap(lemna_geodist_matrix,
             main = "Geographic Distance Heatmap", Colv=NA, Rowv=NA,
             col = colorRampPalette(c("firebrick1", "dodgerblue"))(100))
     
     
     
     
     ## plot depth of stack for each individual ####
     
     ## set directory
     dir = "C:/Users/timte/Desktop/Brisbane/Chapter 1/Duckweed GBS/Second run early 2025/Landoltia_denovo_out/"
     
     # ## make list of zipped files
     # ustacks_zipped = list.files(dir, pattern = "tags.tsv.gz")
     # ustacks_zipped = ustacks_output[-1]
     # 
     # ## extract files
     # for (n in 1:length(ustacks_zipped)) {R.utils::gunzip(paste(dir,ustacks_output[n], sep=""))}
     
     ## make list of unzipped files
     ustacks_output = list.files(dir, pattern = "tags.tsv")
     ustacks_output = ustacks_output[-1]
     
     ## run through files and create plots
     
     for (n in 1:length(ustacks_output)) {
          
          ## load lines from file
          file_lines = readLines(paste(dir, ustacks_output[n], sep=""))
          
          ## split string into list (abomination)
          file_lines = strsplit(file_lines, "\t")
          
          ## extract first number, which is the stacks ID
          stacks_IDs = unlist(lapply(file_lines, function(x) x[[1]]))
          stacks_IDs = as.numeric(stacks_IDs[3:length(stacks_IDs)-1])
          
          ## remove two of each value (Chat GPT)
          tab <- table(stacks_IDs)
          keep_counts <- pmax(tab - 2, 0)
          stacks_IDs <- as.numeric(unlist(mapply(rep, as.numeric(names(keep_counts)), keep_counts)))
          
          ## removes scientific numbers for plotting
          options(scipen=999)
          
          ## generate histogram
          png(paste("Landoltia ",substr(ustacks_output[n], 1, (nchar(ustacks_output[n])-9)), ".png",sep=""), width=700, heigh=700, res=100)
          plot(unique(stacks_IDs), log(sort(table(stacks_IDs))), 
               ylab="within stack copy number (log)", xlab="stacks",
               main=paste("Landoltia ",substr(ustacks_output[n], 1, (nchar(ustacks_output[n])-9)), ".png",sep=""),
               pch=20, col="purple")
          abline(h=log(10), lty=2)
          text(1, log(10), labels="10", pos=4)
          abline(h=log(100), lty=2)
          text(1, log(100), labels="100", pos=4)
          abline(h=log(1000), lty=2)
          text(1, log(1000), labels="1000", pos=4)
          dev.off()
          
     }
     ## plot overlap in stacks (can't work, too much data) ####
     
     dat = read.table("C:/Users/timte/Desktop/Brisbane/Chapter 1/Duckweed GBS/Second run early 2025/Landoltia_denovo_out/R_populations.haplotypes.tsv",
                      header=TRUE, nrows=100000, sep="\t")
     
     ## transform into 1 and 0
     dat[3:185][dat[3:185]!="-"]=1
     dat[3:185][dat[3:185]=="-"]=0
     
     ## transform to numeric and sort plants in decreasing order
     dat[] = lapply(dat, as.numeric)
     dat = dat[,order(colSums(dat), decreasing = TRUE)]
     
     ## sort loci in decreasing order
     dat[] = dat[order(rowSums(dat[,3:183]), decreasing = TRUE),]
     
     ## traditional plotting -> visually collapses plot, doesn't show everything
     cor = which(dat[,3:183] == 1, arr.ind = TRUE)
     png("Landoltia_stack_overlap_new_method.png", width=10000, heigh=7000)
     plot(NA, xlim = c(1, ncol(dat)-3), ylim = c(1, nrow(dat)),
          xlab = "", ylab = "", xaxt = "n", yaxt = "n")
     points(cor[, "col"], cor[, "row"], pch = 15, cex = 0.5)
     dev.off()
     
     ## plot density plot
     plot(density(dat[,"Cnt"]),xlab="number of samples in which each stack is found",
          ylab="Density")
     
     
     
     
     
     
     
     ## Landoltia identify perfect clones based in landoltia_gendist ####
     
     ## transform into matrix of genotypes
     landoltia_genotypes = extract.gt(landoltia_final, element = "GT", as.numeric = TRUE)
     
     ## compute distance, transform to matrix and transform to number of SNP difference
     landoltia_gendist_dist = dist(t(landoltia_genotypes), method = "euclidean")
     landoltia_gendist_matrix = as.matrix(landoltia_gendist_dist)
     landoltia_gendist_matrix = landoltia_gendist_matrix^2
     
     ## find plants that have only 3 SNPs difference
     clonal_candidates = which(landoltia_gendist_matrix == 0, arr.ind = T)
     rownames(clonal_candidates) = NULL
     
     ## remove self-comparison
     clonal_candidates = clonal_candidates[which(clonal_candidates[,1] != clonal_candidates[,2]),]
     
     ## remove duplicates
     n=0; whiler = TRUE
     while (whiler == TRUE) {
          
          ## counter
          n=n+1
          
          ## check for duplicates
          if (sum(clonal_candidates[,1] == clonal_candidates[n,2] & clonal_candidates[,2] == clonal_candidates[n,1])==1) {
               
               ## remove if it does exist
               clonal_candidates = clonal_candidates[-n,]
          }
          
          ## stop loop
          if (n == nrow(clonal_candidates)) {
               whiler = FALSE
          }
          
     }
     
     ## make clonal groups
     clone_list = vector(mode='list'); lc = 0
     for (n in 1:length(sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2]))))))) {
          
          ## list counter
          runner_candidate = sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2])))))[n]
          
          ## start vector
          clone_vector = runner_candidate
          
          ## clone vector
          whiler = TRUE
          while(whiler==TRUE) {
               pre_clone_vector = clone_vector
               clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,1]))),]))
               clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,2]))),]))
               if (identical(pre_clone_vector,clone_vector)){whiler = FALSE}
          }
          
          ## add indices to list
          clone_list[[n]] = sort(clone_vector)
          
     }
     
     ## remove duplicates
     clone_str = sapply(clone_list, paste, collapse = "-")
     clone_list = clone_list[!duplicated(clone_str)]
     
     ## extract sample names
     landoltia_clone_list = lapply(clone_list, function(x) rownames(landoltia_gendist_matrix)[x])
     lengths = unlist(lapply(landoltia_clone_list, length))
     
     ## print as table
     landoltia_out = as.data.frame(matrix(NA, ncol=length(landoltia_clone_list), nrow=max(lengths)))
     
     ## create names
     landoltia_clone_names = vector()
     for (n in 1:length(landoltia_clone_list)) {landoltia_clone_names[n] = paste0("clone",n, sep="")}
     names(landoltia_out) = landoltia_clone_names
     
     ## fill dataframe
     for (n in 1:length(landoltia_clone_list)) {
          
          ## transform to vector and put in matrix
          runner = landoltia_clone_list[[n]]
          runner = c(runner, rep(NA, max(lengths)-length(runner)))
          landoltia_out[,n] = runner
          
     }
     
     ## export table
     write.csv(landoltia_out, "Landoltia_clones.csv")
     
     ## Lemna identify perfect clones baes on lemna_gendist ####
     
     ## transform into matrix of genotypes
     lemna_genotypes = extract.gt(lemna_final, element = "GT", as.numeric = TRUE)
     
     ## compute distance, transform to matrix and transform to number of SNP difference
     lemna_gendist_dist = dist(t(lemna_genotypes), method = "euclidean")
     lemna_gendist_matrix = as.matrix(lemna_gendist_dist)
     lemna_gendist_matrix = lemna_gendist_matrix^2
     
     ## find plants that have only 3 SNPs difference
     clonal_candidates = which(lemna_gendist_matrix == 0, arr.ind = T)
     rownames(clonal_candidates) = NULL
     
     ## remove self-comparison
     clonal_candidates = clonal_candidates[which(clonal_candidates[,1] != clonal_candidates[,2]),]
     
     ## remove duplicates
     n=0; whiler = TRUE
     while (whiler == TRUE) {
          
          ## counter
          n=n+1
          
          ## check for duplicates
          if (sum(clonal_candidates[,1] == clonal_candidates[n,2] & clonal_candidates[,2] == clonal_candidates[n,1])==1) {
               
               ## remove if it does exist
               clonal_candidates = clonal_candidates[-n,]
          }
          
          ## stop loop
          if (n == nrow(clonal_candidates)) {
               whiler = FALSE
          }
          
     }
     
     ## make clonal groups
     clone_list = vector(mode='list'); lc = 0
     for (n in 1:length(sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2]))))))) {
          
          ## list counter
          runner_candidate = sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2])))))[n]
          
          ## start vector
          clone_vector = runner_candidate
          
          ## clone vector
          whiler = TRUE
          while(whiler==TRUE) {
               pre_clone_vector = clone_vector
               clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,1]))),]))
               clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,2]))),]))
               if (identical(pre_clone_vector,clone_vector)){whiler = FALSE}
          }
          
          ## add indices to list
          clone_list[[n]] = sort(clone_vector)
          
     }
     
     ## remove duplicates
     clone_str = sapply(clone_list, paste, collapse = "-")
     clone_list = clone_list[!duplicated(clone_str)]
     
     ## extract sample names
     lemna_clone_list = lapply(clone_list, function(x) rownames(lemna_gendist_matrix)[x])
     lengths = unlist(lapply(lemna_clone_list, length))
     
     ## print as table
     lemna_out = as.data.frame(matrix(NA, ncol=length(lemna_clone_list), nrow=max(lengths)))
     
     ## create names
     lemna_clone_names = vector()
     for (n in 1:length(lemna_clone_list)) {lemna_clone_names[n] = paste0("clone",n, sep="")}
     names(lemna_out) = lemna_clone_names
     
     ## fill dataframe
     for (n in 1:length(lemna_clone_list)) {
          
          ## transform to vector and put in matrix
          runner = lemna_clone_list[[n]]
          runner = c(runner, rep(NA, max(lengths)-length(runner)))
          lemna_out[,n] = runner
          
     }
     
     ## export table
     write.csv(lemna_out, "Lemna_clones.csv")
     
     ## landoltia clone heatmap gendist ####
     
     ## hardcode meaningful order of samples
     landoltia_ordered_names = c("P1S2",
                                 "P2S2",
                                 "P4S2",
                                 "P10S1b","P10S3","P10S3b","P10S4","P10S4b","P10S6","P10S7b","P10S8b","P10S9b","P10S10b","P10S11b","P10S12b","P10S13b","P10S14b","P10S15","P10S16","P10S18","P10S22","P10S23","P10S26","P10S31",
                                 "P11S10","P11S11","P11S12",
                                 "P12S4","P12S5","P12S6",
                                 "P13S1","P13S2","P13S3",
                                 "P14S1","P14S7","P14S8","P14S9","P14S12","P14S13","P14S16","P14S17","P14S18","P14S22","P14S23","P14S24","P14S28","P14S29","P14S30","P14S37","P14S38","P14S39","P14S40","P14S42","P14S46","P14S47","P14S48",
                                 "P15S1","P15S2","P15S3",
                                 "P16S4","P16S5","P16S6",
                                 "P17S4","P17S6",
                                 "P18S4","P18S6",
                                 "P19S4","P19S5","P19S10","P19S11","P19S16","P19S17","P19S18","P19S22","P19S23","P19S24","P19S25","P19S26","P19S27","P19S34","P19S35","P19S36","P19S40","P19S41","P19S42","P19S43","P19S44","P19S45","P19S52","P19S53","P19S54",
                                 "P20S1","P20S2","P20S3",
                                 "P22S1","P22S2","P22S3",
                                 "P23S4","P23S5",
                                 "P24S1","P24S2","P24S3",
                                 "P25S1","P25S2","P25S3",
                                 "P26S1","P26S2","P26S3",
                                 "P27S1","P27S2","P27S3","P27S7","P27S8","P27S9","P27S13","P27S14","P27S15","P27S22","P27S23","P27S24","P27S26","P27S28","P27S30","P27S31","P27S32","P27S33","P27S34","P27S35","P27S36","P27S37","P27S38","P27S49","P27S51",
                                 "P28S1","P28S2","P28S3",
                                 "P32S1","P32S2","P32S3",
                                 "P34S1","P34S2","P34S3",
                                 "P36S1","P36S2","P36S3","P36S7","P36S8","P36S9","P36S16","P36S17","P36S18","P36S19","P36S20","P36S21","P36S25","P36S26","P36S27","P36S34","P36S35","P36S36","P36S40","P36S41")
     
     ## transform into matrix of genotypes
     landoltia_genotypes = extract.gt(landoltia_final, element = "GT", as.numeric = TRUE)
     
     ## compute distance (in dist object)
     landoltia_gendist_dist = dist(t(landoltia_genotypes), method = "euclidean")
     
     ## transform to matrix and order samples
     landoltia_gendist_matrix = as.matrix(landoltia_gendist_dist)
     landoltia_gendist_matrix = landoltia_gendist_matrix[landoltia_ordered_names, landoltia_ordered_names]
     
     ## set up split screen layout
     split.screen(rbind(c(0, 0.8, 0, 1), 
                        c(0.8, 1, 0, 1)))
     
     ## create palette
     purple_white = colorRampPalette(c("purple", "thistle1"))
     
     ## create heatmap
     screen(1)
     image(landoltia_gendist_matrix, axes = FALSE, frame=FALSE, 
           col=ifelse(landoltia_gendist_matrix == 0,"purple", "thistle1"),
           main="Landoltia clones")
     ## add axis
     axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(landoltia_gendist_matrix)), 
          labels = landoltia_ordered_names, las =2, cex.axis=0.8)
     axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(landoltia_gendist_matrix)), 
          labels = landoltia_ordered_names, las =2, cex.axis=0.8)
     close.screen(1)
     
     ## add legend
     screen(2)
     par(mar=c(2,1,2,4))
     plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,2),xaxt="n",yaxt="n",bty="n")
     for (n in 1:2) {rect(1,n-1,2,n, col=purple_white(2)[n])}
     axis(side=4, at=c(1:2)-0.5, labels=c("clone", "non-clone"), las=2, lwd=0, lwd.ticks=1)
     close.screen(2)
     close.screen(all.screens=TRUE)
     
     
     ## lemna clone heatmap gendist ####
     
     ## hardcode meaningful order of samples
     lemna_ordered_names = c("P5S1",
                             "P6S1","P6S3",
                             "P7S1","P7S2","P7S3",
                             "P10S1","P10S2","P10S2b","P10S5","P10S5b","P10S6b","P10S7","P10S8","P10S9","P10S10","P10S11","P10S12","P10S13","P10S14","P10S15b","P10S19","P10S20","P10S21","P10S28","P10S30",
                             "P11S7","P11S8","P11S9",
                             "P14S10","P14S11","P14S15","P14S19","P14S20","P14S21","P14S25","P14S26","P14S27","P14S31","P14S32","P14S34","P14S35","P14S36","P14S43","P14S44","P14S45","P14S49","P14S50","P14S55","P14S56","P14S57",
                             "P16S2","P16S3",
                             "P17S1","P17S2","P17S3",
                             "P18S2",
                             "P19S1","P19S2","P19S3","P19S7","P19S8","P19S9","P19S13","P19S14","P19S15","P19S19","P19S20","P19S21","P19S28","P19S29","P19S30","P19S31","P19S32","P19S33","P19S37","P19S38","P19S39","P19S46","P19S47","P19S48","P19S49","P19S50","P19S51",
                             "P21S1","P21S2",
                             "P22S4","P22S6",
                             "P23S3",
                             "P27S4","P27S5","P27S10","P27S11","P27S12","P27S16","P27S17","P27S18","P27S19","P27S20","P27S21","P27S25","P27S27","P27S39","P27S40","P27S41","P27S42","P27S43","P27S44","P27S45","P27S46","P27S47","P27S48","P27S52",
                             "P28S4","P28S5","P28S6",
                             "P30S1","P30S2","P30S3",
                             "P31S1","P31S2","P31S3",
                             "P32S4","P32S5","P32S6",
                             "P33S1","P33S2","P33S3",
                             "P34S5","P34S6",
                             "P35S1","P35S2","P35S3",
                             "P36S4","P36S5","P36S6","P36S10","P36S11","P36S12","P36S13","P36S14","P36S15","P36S22","P36S23","P36S24","P36S28","P36S29","P36S30","P36S31","P36S32","P36S33","P36S37","P36S38","P36S39")
     
     ## transform into matrix of genotypes
     lemna_genotypes = extract.gt(lemna_final, element = "GT", as.numeric = TRUE)
     
     ## compute distance (in dist object)
     lemna_gendist_dist = dist(t(lemna_genotypes), method = "euclidean")
     
     ## transform to matrix and order samples
     lemna_gendist_matrix = as.matrix(lemna_gendist_dist)
     lemna_gendist_matrix = lemna_gendist_matrix[lemna_ordered_names, lemna_ordered_names]
     
     ## set up split screen layout
     split.screen(rbind(c(0, 0.8, 0, 1), 
                        c(0.8, 1, 0, 1)))
     
     ## create palette
     purple_white = colorRampPalette(c("darkgreen", "darkseagreen1"))
     
     ## create heatmap
     screen(1)
     image(lemna_gendist_matrix, axes = FALSE, frame=FALSE, 
           col=ifelse(lemna_gendist_matrix == 0,"darkgreen", "darkseagreen1"),
           main="Lemna clones")
     ## add axis
     axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(lemna_gendist_matrix)), 
          labels = lemna_ordered_names, las =2, cex.axis=0.8)
     axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(lemna_gendist_matrix)), 
          labels = lemna_ordered_names, las =2, cex.axis=0.8)
     close.screen(1)
     
     ## add legend
     screen(2)
     par(mar=c(2,1,2,4))
     plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,2),xaxt="n",yaxt="n",bty="n")
     for (n in 1:2) {rect(1,n-1,2,n, col=purple_white(2)[n])}
     axis(side=4, at=c(1:2)-0.5, labels=c("clone", "non-clone"), las=2, lwd=0, lwd.ticks=1)
     close.screen(2)
     close.screen(all.screens=TRUE)
     
     ## landoltia kinship heatmap Eucledian ####
     
     ## hardcode meaningful order of samples
     landoltia_ordered_names = c("P1S2",
                                 "P2S2",
                                 "P4S2",
                                 "P10S1b","P10S3","P10S3b","P10S4","P10S4b","P10S6","P10S7b","P10S8b","P10S9b","P10S10b","P10S11b","P10S12b","P10S13b","P10S14b","P10S15","P10S16","P10S18","P10S22","P10S23","P10S26","P10S31",
                                 "P11S10","P11S11","P11S12",
                                 "P12S4","P12S5","P12S6",
                                 "P13S1","P13S2","P13S3",
                                 "P14S1","P14S7","P14S8","P14S9","P14S12","P14S13","P14S16","P14S17","P14S18","P14S22","P14S23","P14S24","P14S28","P14S29","P14S30","P14S37","P14S38","P14S39","P14S40","P14S42","P14S46","P14S47","P14S48",
                                 "P15S1","P15S2","P15S3",
                                 "P16S4","P16S5","P16S6",
                                 "P17S4","P17S6",
                                 "P18S4","P18S6",
                                 "P19S4","P19S5","P19S10","P19S11","P19S16","P19S17","P19S18","P19S22","P19S23","P19S24","P19S25","P19S26","P19S27","P19S34","P19S35","P19S36","P19S40","P19S41","P19S42","P19S43","P19S44","P19S45","P19S52","P19S53","P19S54",
                                 "P20S1","P20S2","P20S3",
                                 "P22S1","P22S2","P22S3",
                                 "P23S4","P23S5",
                                 "P24S1","P24S2","P24S3",
                                 "P25S1","P25S2","P25S3",
                                 "P26S1","P26S2","P26S3",
                                 "P27S1","P27S2","P27S3","P27S7","P27S8","P27S9","P27S13","P27S14","P27S15","P27S22","P27S23","P27S24","P27S26","P27S28","P27S30","P27S31","P27S32","P27S33","P27S34","P27S35","P27S36","P27S37","P27S38","P27S49","P27S51",
                                 "P28S1","P28S2","P28S3",
                                 "P32S1","P32S2","P32S3",
                                 "P34S1","P34S2","P34S3",
                                 "P36S1","P36S2","P36S3","P36S7","P36S8","P36S9","P36S16","P36S17","P36S18","P36S19","P36S20","P36S21","P36S25","P36S26","P36S27","P36S34","P36S35","P36S36","P36S40","P36S41")
     
     ## transform into matrix of genotypes
     landoltia_genotypes = extract.gt(landoltia_final, element = "GT", as.numeric = TRUE)
     
     ## compute distance (in dist object)
     landoltia_gendist_dist = dist(t(landoltia_genotypes), method = "euclidean")
     
     ## transform to matrix and order samples
     landoltia_gendist_matrix = as.matrix(landoltia_gendist_dist)
     landoltia_gendist_matrix = landoltia_gendist_matrix[landoltia_ordered_names, landoltia_ordered_names]
     
     ## create custom heatmap
     ## make gradient
     purple_white = colorRampPalette(c("purple", "white"))
     landoltia_number_breaks = 10
     
     ## set up split screen layout
     #par(mar=c(1,1,1,1),oma=c(3,3,3,3))
     split.screen(rbind(c(0, 0.8, 0, 1), 
                        c(0.8, 1, 0, 1)))
     
     ## create heatmap
     screen(1)
     image(landoltia_gendist_matrix, axes = FALSE, frame=FALSE, 
           col=c("black",purple_white(landoltia_number_breaks-2)),
           breaks=c(0,seq(from=min(landoltia_gendist_matrix)+0.01, to=max(landoltia_gendist_matrix), length.out = landoltia_number_breaks-1)),
           main="Landoltia Euclidean SNP Distance")
     ## add axis
     axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(landoltia_gendist_matrix)), 
          labels = landoltia_ordered_names, las =2, cex.axis=0.8)
     axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(landoltia_gendist_matrix)), 
          labels = landoltia_ordered_names, las =2, cex.axis=0.8)
     close.screen(1)
     
     ## add legend
     screen(2)
     par(mar=c(2,1,2,4))
     plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,landoltia_number_breaks),xaxt="n",yaxt="n",bty="n")
     for (n in 1:(landoltia_number_breaks-1)) {rect(1,n-1,2,n, col=c("black",purple_white(landoltia_number_breaks-2))[n])}
     axis(side=4, at=c(1:(landoltia_number_breaks-1))-0.5, labels=round(seq(from=min(landoltia_gendist_matrix), to=max(landoltia_gendist_matrix), length.out=landoltia_number_breaks-1),0), las=2, lwd=0, lwd.ticks=1)
     close.screen(2)
     close.screen(all.screens=TRUE)
     
     ## lemna kinship heatmap Eucledian ####
     
     lemna_ordered_names = c("P5S1",
                             "P6S1","P6S3",
                             "P7S1","P7S2","P7S3",
                             "P10S1","P10S2","P10S2b","P10S5","P10S5b","P10S6b","P10S7","P10S8","P10S9","P10S10","P10S11","P10S12","P10S13","P10S14","P10S15b","P10S19","P10S20","P10S21","P10S28","P10S30",
                             "P11S7","P11S8","P11S9",
                             "P14S10","P14S11","P14S15","P14S19","P14S20","P14S21","P14S25","P14S26","P14S27","P14S31","P14S32","P14S34","P14S35","P14S36","P14S43","P14S44","P14S45","P14S49","P14S50","P14S55","P14S56","P14S57",
                             "P16S2","P16S3",
                             "P17S1","P17S2","P17S3",
                             "P18S2",
                             "P19S1","P19S2","P19S3","P19S7","P19S8","P19S9","P19S13","P19S14","P19S15","P19S19","P19S20","P19S21","P19S28","P19S29","P19S30","P19S31","P19S32","P19S33","P19S37","P19S38","P19S39","P19S46","P19S47","P19S48","P19S49","P19S50","P19S51",
                             "P21S1","P21S2",
                             "P22S4","P22S6",
                             "P23S3",
                             "P27S4","P27S5","P27S10","P27S11","P27S12","P27S16","P27S17","P27S18","P27S19","P27S20","P27S21","P27S25","P27S27","P27S39","P27S40","P27S41","P27S42","P27S43","P27S44","P27S45","P27S46","P27S47","P27S48","P27S52",
                             "P28S4","P28S5","P28S6",
                             "P30S1","P30S2","P30S3",
                             "P31S1","P31S2","P31S3",
                             "P32S4","P32S5","P32S6",
                             "P33S1","P33S2","P33S3",
                             "P34S5","P34S6",
                             "P35S1","P35S2","P35S3",
                             "P36S4","P36S5","P36S6","P36S10","P36S11","P36S12","P36S13","P36S14","P36S15","P36S22","P36S23","P36S24","P36S28","P36S29","P36S30","P36S31","P36S32","P36S33","P36S37","P36S38","P36S39")
     
     
     ## transform into matrix of genotypes
     lemna_genotypes = extract.gt(lemna_final, element = "GT", as.numeric = TRUE)
     
     ## compute distance (in dist object)
     lemna_gendist_dist = dist(t(lemna_genotypes), method = "euclidean")
     
     ## transform to matrix and order samples
     lemna_gendist_matrix = as.matrix(lemna_gendist_dist)
     lemna_gendist_matrix = lemna_gendist_matrix[lemna_ordered_names, lemna_ordered_names]
     
     ## create custom heatmap
     ## make gradient
     green_white = colorRampPalette(c("darkgreen", "white"))
     lemna_number_breaks = 10
     
     ## set up split screen layout
     #par(mar=c(1,1,1,1),oma=c(3,3,3,3))
     split.screen(rbind(c(0, 0.8, 0, 1), 
                        c(0.8, 1, 0, 1)))
     
     ## create heatmap
     screen(1)
     image(lemna_gendist_matrix, axes = FALSE, frame=FALSE, 
           col=c("black",green_white(lemna_number_breaks-2)),
           breaks=c(0,seq(from=min(lemna_gendist_matrix)+0.01, to=max(lemna_gendist_matrix), length.out = lemna_number_breaks-1)),
           main="Lemna Euclidean SNP Distance")
     ## add axis
     axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(lemna_gendist_matrix)), 
          labels = lemna_ordered_names, las =2, cex.axis=0.8)
     axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(lemna_gendist_matrix)), 
          labels = lemna_ordered_names, las =2, cex.axis=0.8)
     close.screen(1)
     
     ## add legend
     screen(2)
     par(mar=c(2,1,2,4))
     plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,lemna_number_breaks),xaxt="n",yaxt="n",bty="n")
     for (n in 1:(lemna_number_breaks-1)) {rect(1,n-1,2,n, col=c("black",green_white(lemna_number_breaks-2))[n])}
     axis(side=4, at=c(1:(lemna_number_breaks-1))-0.5, labels=round(seq(from=min(lemna_gendist_matrix), to=max(lemna_gendist_matrix), length.out=lemna_number_breaks-1),0), las=2, lwd=0, lwd.ticks=1)
     close.screen(2)
     close.screen(all.screens=TRUE)
     
     
     ## LANDOLTIA clone heatmap ####
     
     ## set up split screen layout
     split.screen(rbind(c(0, 0.8, 0, 1), 
                        c(0.8, 1, 0, 1)))
     
     ## create heatmap
     screen(1)
     image(landoltia_hamdist, axes = FALSE, frame=FALSE, 
           col=ifelse(landoltia_hamdist <= 35,"purple", "thistle1"),
           main="Landoltia clones")
     ## add axis
     axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(landoltia_hamdist)), 
          labels = landoltia_ordered_names, las =2, cex.axis=0.8)
     axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(landoltia_hamdist)), 
          labels = landoltia_ordered_names, las =2, cex.axis=0.8)
     close.screen(1)
     
     ## add legend
     screen(2)
     par(mar=c(2,1,2,4))
     plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,2),xaxt="n",yaxt="n",bty="n")
     for (n in 1:2) {rect(1,n-1,2,n, col=c("purple", "thistle1")[n])}
     axis(side=4, at=c(1:2)-0.5, labels=c("clone", "non-clone"), las=2, lwd=0, lwd.ticks=1)
     close.screen(2)
     close.screen(all.screens=TRUE)
     
     ## LEMNA clone heatmap ####
     
     
     ## transform into matrix of genotypes
     lemna_genotypes = extract.gt(lemna_final, element = "GT", as.numeric = TRUE)
     
     
     
     ## set up split screen layout
     split.screen(rbind(c(0, 0.8, 0, 1), 
                        c(0.8, 1, 0, 1)))
     
     ## create heatmap
     screen(1)
     image(lemna_hamdist, axes = FALSE, frame=FALSE, 
           col=ifelse(lemna_hamdist <= 35,"darkgreen", "lightgreen"),
           main="lemna clones")
     ## add axis
     axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(lemna_hamdist)), 
          labels = lemna_ordered_names, las =2, cex.axis=0.8)
     axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(lemna_hamdist)), 
          labels = lemna_ordered_names, las =2, cex.axis=0.8)
     close.screen(1)
     
     ## add legend
     screen(2)
     par(mar=c(2,1,2,4))
     plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,2),xaxt="n",yaxt="n",bty="n")
     for (n in 1:2) {rect(1,n-1,2,n, col=c("darkgreen", "lightgreen")[n])}
     axis(side=4, at=c(1:2)-0.5, labels=c("clone", "non-clone"), las=2, lwd=0, lwd.ticks=1)
     close.screen(2)
     close.screen(all.screens=TRUE)
     
     ## landoltia living plants heatmap ####
     
     
     
     ## gendist matrix for living plants
     landoltia_alive_gendist_matrix = landoltia_gendist_matrix[which(names(landoltia_gendist_matrix[1,]) %in% landoltia_alive),which(names(landoltia_gendist_matrix[1,]) %in% landoltia_alive)]
     landoltia_alive_and_geneotyped = intersect(landoltia_alive, colnames(landoltia_alive_gendist_matrix))
     landoltia_alive_gendist_matrix = landoltia_alive_gendist_matrix[landoltia_alive_and_geneotyped, landoltia_alive_and_geneotyped]
     
     ## create custom heatmap
     ## make gradient
     green_white = colorRampPalette(c("purple", "white"))
     landoltia_number_breaks = 15
     
     ## set up split screen layout
     #par(mar=c(1,1,1,1),oma=c(3,3,3,3))
     split.screen(rbind(c(0, 0.8, 0, 1), 
                        c(0.8, 1, 0, 1)))
     
     ## create heatmap
     screen(1)
     image(landoltia_alive_gendist_matrix, axes = FALSE, frame=FALSE, 
           col=green_white(landoltia_number_breaks-1),breaks=seq(from=min(landoltia_alive_gendist_matrix), to=max(landoltia_alive_gendist_matrix), length.out = landoltia_number_breaks),
           main="Landoltia Euclidean SNP Distance")
     ## add axis
     axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(landoltia_alive_gendist_matrix)), 
          labels = colnames(landoltia_alive_gendist_matrix), las =2, cex.axis=0.8)
     axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(landoltia_alive_gendist_matrix)), 
          labels = rownames(landoltia_alive_gendist_matrix), las =2, cex.axis=0.8)
     close.screen(1)
     
     ## add legend
     screen(2)
     par(mar=c(2,1,2,4))
     plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,landoltia_number_breaks),xaxt="n",yaxt="n",bty="n")
     for (n in 1:(landoltia_number_breaks-1)) {rect(1,n-1,2,n, col=green_white(landoltia_number_breaks-1)[n])}
     axis(side=4, at=c(1:(landoltia_number_breaks-1))-0.5, labels=round(seq(from=min(landoltia_alive_gendist_matrix), to=max(landoltia_alive_gendist_matrix), length.out=landoltia_number_breaks-1),0), las=2, lwd=0, lwd.ticks=1)
     close.screen(2)
     close.screen(all.screens=TRUE)
     
     ## lemna living plants heatmap ####
     
     
     
     ## gendist matrix for living plants
     lemna_alive_gendist_matrix = lemna_gendist_matrix[which(names(lemna_gendist_matrix[1,]) %in% lemna_alive),which(names(lemna_gendist_matrix[1,]) %in% lemna_alive)]
     lemna_alive_and_geneotyped = intersect(lemna_alive, colnames(lemna_alive_gendist_matrix))
     lemna_alive_gendist_matrix = lemna_alive_gendist_matrix[lemna_alive_and_geneotyped, lemna_alive_and_geneotyped]
     
     ## create custom heatmap
     ## make gradient
     green_white = colorRampPalette(c("darkgreen", "white"))
     lemna_number_breaks = 15
     
     ## set up split screen layout
     #par(mar=c(1,1,1,1),oma=c(3,3,3,3))
     split.screen(rbind(c(0, 0.8, 0, 1), 
                        c(0.8, 1, 0, 1)))
     
     ## create heatmap
     screen(1)
     image(lemna_alive_gendist_matrix, axes = FALSE, frame=FALSE, 
           col=green_white(lemna_number_breaks-1),breaks=seq(from=min(lemna_alive_gendist_matrix), to=max(lemna_alive_gendist_matrix), length.out = lemna_number_breaks),
           main="Lemna Euclidean SNP Distance")
     ## add axis
     axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(lemna_alive_gendist_matrix)), 
          labels = colnames(lemna_alive_gendist_matrix), las =2, cex.axis=0.8)
     axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(lemna_alive_gendist_matrix)), 
          labels = rownames(lemna_alive_gendist_matrix), las =2, cex.axis=0.8)
     close.screen(1)
     
     ## add legend
     screen(2)
     par(mar=c(2,1,2,4))
     plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,lemna_number_breaks),xaxt="n",yaxt="n",bty="n")
     for (n in 1:(lemna_number_breaks-1)) {rect(1,n-1,2,n, col=green_white(lemna_number_breaks-1)[n])}
     axis(side=4, at=c(1:(lemna_number_breaks-1))-0.5, labels=round(seq(from=min(lemna_alive_gendist_matrix), to=max(lemna_alive_gendist_matrix), length.out=lemna_number_breaks-1),0), las=2, lwd=0, lwd.ticks=1)
     close.screen(2)
     close.screen(all.screens=TRUE)
     
     
     ## LEMNA compute hamming distance ####
     
     ## transform into matrix of genotypes
     lemna_genotypes = extract.gt(lemna_final, element = "GT", return.alleles = TRUE)
     
     ## create matrix for manual hamming distance calculation
     lemna_hamdist = matrix(NA, nrow=ncol(lemna_genotypes), ncol=ncol(lemna_genotypes))
     colnames(lemna_hamdist) = colnames(lemna_genotypes)
     rownames(lemna_hamdist) = colnames(lemna_genotypes)
     
     ##  compute hamming distance
     for (n in 1:ncol(lemna_genotypes)) {
       
       ## compute hamming distance
       runner_vec = apply(lemna_genotypes, 2, function(x) sum(lemna_genotypes[,n] != x, na.rm=TRUE))
       
       ## store in matrix
       lemna_hamdist[,n] = runner_vec
       
     }
     
     ## hardcode meaningful order of samples
     lemna_ordered_names = c("P5S1",
                             "P6S1","P6S3",
                             "P7S1","P7S2","P7S3",
                             "P10S1","P10S2","P10S2b","P10S5","P10S5b","P10S6b","P10S7","P10S8","P10S9","P10S10","P10S11","P10S12","P10S13","P10S14","P10S15b","P10S19","P10S20","P10S21","P10S28","P10S30",
                             "P11S7","P11S8","P11S9",
                             "P14S10","P14S11","P14S15","P14S19","P14S20","P14S21","P14S25","P14S26","P14S27","P14S31","P14S32","P14S34","P14S35","P14S36","P14S43","P14S44","P14S45","P14S49","P14S50","P14S55","P14S56","P14S57",
                             "P16S2","P16S3",
                             "P17S1","P17S2","P17S3",
                             "P18S2",
                             "P19S1","P19S2","P19S3","P19S7","P19S8","P19S9","P19S13","P19S14","P19S15","P19S19","P19S20","P19S21","P19S28","P19S29","P19S30","P19S31","P19S32","P19S33","P19S37","P19S38","P19S39","P19S46","P19S47","P19S48","P19S49","P19S50","P19S51",
                             "P21S1","P21S2",
                             "P22S4","P22S6",
                             "P23S3",
                             "P27S4","P27S5","P27S10","P27S11","P27S12","P27S16","P27S17","P27S18","P27S19","P27S20","P27S21","P27S25","P27S27","P27S39","P27S40","P27S41","P27S42","P27S43","P27S44","P27S45","P27S46","P27S47","P27S48","P27S52",
                             "P28S4","P28S5","P28S6",
                             "P30S1","P30S2","P30S3",
                             "P31S1","P31S2","P31S3",
                             "P32S4","P32S5","P32S6",
                             "P33S1","P33S2","P33S3",
                             "P34S5","P34S6",
                             "P35S1","P35S2","P35S3",
                             "P36S4","P36S5","P36S6","P36S10","P36S11","P36S12","P36S13","P36S14","P36S15","P36S22","P36S23","P36S24","P36S28","P36S29","P36S30","P36S31","P36S32","P36S33","P36S37","P36S38","P36S39")
     
     ## order hamming distance matrix 
     lemna_hamdist = lemna_hamdist[lemna_ordered_names, lemna_ordered_names]
     
     ## LEMNA identify clones based on distribution ####
     
     ## visualize the distribution
     par(mfrow=c(1,2))
     hist(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)], breaks=100,
          main = "lemna hamming distance distribution", col="darkgreen", xlab="")
     abline(v=75, col="red", lty=2)
     hist(sort(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)])[1:5000], breaks=100,
          main = "lemna zoom in", col="darkgreen", xlab="")
     abline(v=75, col="red", lty=2)
     lemna_clone_cutoff = 75
     
     ## find plants with 75 SNP difference
     clonal_candidates = which(lemna_hamdist <= lemna_clone_cutoff, arr.ind = T)
     rownames(clonal_candidates) = NULL
     
     ## remove self-comparison
     clonal_candidates = clonal_candidates[which(clonal_candidates[,1] != clonal_candidates[,2]),]
     
     ## remove duplicates
     n=0; whiler = TRUE
     while (whiler == TRUE) {
       
       ## counter
       n=n+1
       
       ## check for duplicates
       if (n >= nrow(clonal_candidates)) {
         ## stop loop
         whiler = FALSE
       }
       else if (sum(clonal_candidates[,1] == clonal_candidates[n,2] & clonal_candidates[,2] == clonal_candidates[n,1])==1) {
         ## remove if it does exist
         clonal_candidates = clonal_candidates[-n,]
       }
     }
     
     ## make clonal groups
     clone_list = vector(mode='list'); lc = 0
     for (n in 1:length(sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2]))))))) {
       
       ## list counter
       runner_candidate = sort(unique(c(unique(clonal_candidates[,1], unique(clonal_candidates[,2])))))[n]
       
       ## start vector
       clone_vector = runner_candidate
       
       ## clone vector
       whiler = TRUE
       while(whiler==TRUE) {
         pre_clone_vector = clone_vector
         clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,1]))),]))
         clone_vector = unique(c(clone_vector,clonal_candidates[unlist(sapply(clone_vector, function(x) which(x==clonal_candidates[,2]))),]))
         if (identical(pre_clone_vector,clone_vector)){whiler = FALSE}
       }
       
       ## add indices to list
       clone_list[[n]] = sort(clone_vector)
       
     }
     
     ## remove duplicates
     clone_str = sapply(clone_list, paste, collapse = "-")
     clone_list = clone_list[!duplicated(clone_str)]
     
     ## extract sample names
     lemna_clone_list = lapply(clone_list, function(x) rownames(lemna_hamdist)[x])
     lengths = unlist(lapply(lemna_clone_list, length))
     
     ## save in dataframe
     lemna_clone_df = as.data.frame(matrix(NA, ncol=length(lemna_clone_list), nrow=max(lengths)))
     lemna_clone_names = vector()
     for (n in 1:length(lemna_clone_list)) {lemna_clone_names[n] = paste0("clone",n, sep="")}
     names(lemna_clone_df) = lemna_clone_names
     for (n in 1:length(lemna_clone_list)) {
       
       ## transform to vector and put in matrix
       runner = lemna_clone_list[[n]]
       runner = c(runner, rep(NA, max(lengths)-length(runner)))
       lemna_clone_df[,n] = runner
       
     }
     
     ## export table
     # write.csv(lemna_clone_df, paste("lemna_clones_hamming0-",lemna_clone_cutoff,".csv", sep=""))
     
     ## LEMNA create clone, no clone and no clone + alive hamming distance ####
     
     ## remove clones
     lemna_noclone = lemna_final
     for (n in 1:ncol(lemna_clone_df)) {
       
       ## run through lemna_clone_df by col
       runner_clone = lemna_clone_df[,n]; runner_clone = runner_clone[!is.na(runner_clone)]
       
       # identify and remove clones
       to_keep = which(!colnames(lemna_noclone@gt) %in% runner_clone[2:length(runner_clone)])
       lemna_noclone = lemna_noclone[, to_keep]
       
     }
     lemna_noclone_names = colnames(lemna_noclone@gt)[-1]
     
     ## adjust hamdist
     lemna_noclone_hamdist = lemna_hamdist[which(names(lemna_hamdist[1,]) %in% lemna_noclone_names),
                                           which(names(lemna_hamdist[1,]) %in% lemna_noclone_names)]     
     
     ## lemna still in collection
     lemna_alive = c("P7S2", "P7S3",
                     "P10S28",
                     "P11S8",
                     "P14S25", "P14S26", "P14S43", "P14S45",
                     "P19S19", "P19S20", "P19S29", "P19S30", "P19S31", "P19S38", "P19S39", "P19S49",
                     "P23S1", "P23S2", "P23S3", "P23S8",
                     "P27S4", "P27S10", "P27S42","P27S52",
                     "P29S1",
                     "P31S3",
                     "P35S2", "P35S3")
     
     ## gendist matrix for living plants
     lemna_noclone_alive_hamdist = lemna_noclone_hamdist[which(names(lemna_noclone_hamdist[1,]) %in% lemna_alive),
                                                         which(names(lemna_noclone_hamdist[1,]) %in% lemna_alive)]
     
     
     ## LANDOLTIA kinship heatmaps ####
     
     ## landoltia with clones
     kinship_heatmap(dist_ma = landoltia_hamdist, species = "landoltia", colbreaks = 10, clone_cutoff = landoltia_clone_cutoff)
     
     ## landoltia without clones
     kinship_heatmap(dist_ma = landoltia_noclone_hamdist, species = "landoltia", colbreaks = 10, clone_cutoff = landoltia_clone_cutoff)
     
     ## landoltia without clones and alive
     kinship_heatmap(dist_ma = landoltia_noclone_alive_hamdist, species = "landoltia", colbreaks = 10, clone_cutoff = landoltia_clone_cutoff)
     
     
     
     ## LEMNA kinship heatmaps ####
     
     ## landoltia with clones
     kinship_heatmap(dist_ma = lemna_hamdist, species = "lemna", colbreaks = 10, clone_cutoff = lemna_clone_cutoff)
     
     ## landoltia without clones
     kinship_heatmap(dist_ma = lemna_noclone_hamdist, species = "lemna", colbreaks = 10, clone_cutoff = lemna_clone_cutoff)
     
     ## landoltia without clones and alive
     kinship_heatmap(dist_ma = lemna_noclone_alive_hamdist, species = "lemna", colbreaks = 10, clone_cutoff = lemna_clone_cutoff)
     
     
     ## create custom heatmap ####
     ## create custom heatmap
     kinship_heatmap = function(dist_ma, species, colbreaks, clone_cutoff) {
       
       ## make colour gradient
       if (species == "landoltia") {
         col_pal = colorRampPalette(c("purple", "white"))
       } else if (species == "lemna") {
         col_pal = colorRampPalette(c("darkgreen", "white"))
       }
       
       ## set up split screen layout
       #par(mar=c(1,1,1,1),oma=c(3,3,3,3))
       split.screen(rbind(c(0, 0.8, 0, 1), 
                          c(0.8, 1, 0, 1)))
       
       ## create heatmap
       screen(1)
       image(dist_ma, axes = FALSE, frame=FALSE, 
             col=c("black",col_pal(colbreaks-2)),
             breaks=c(0,seq(from=clone_cutoff, to=max(dist_ma), length.out = colbreaks-1)))
       ## add axis
       axis(side=1, at=seq(from=-0.00310559, to=1.00310559, length.out = ncol(dist_ma)), 
            labels = colnames(dist_ma), las =2, cex.axis=0.8)
       axis(side=2, at=seq(from=-0.00310559, to=1.00310559, length.out = nrow(dist_ma)), 
            labels = colnames(dist_ma), las =2, cex.axis=0.8)
       close.screen(1)
       
       ## add legend
       screen(2)
       par(mar=c(2,1,2,4))
       plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(0,colbreaks),xaxt="n",yaxt="n",bty="n")
       for (n in 1:(colbreaks-1)) {rect(1,n-1,2,n, col=c("black",col_pal(colbreaks-2))[n])}
       axis(side=4, at=c(1:(colbreaks-1))-0.5, labels=round(seq(from=min(dist_ma), to=max(dist_ma), length.out=colbreaks-1),2), las=2, lwd=0, lwd.ticks=1)
       close.screen(2)
       close.screen(all.screens=TRUE)
     }
     
     
     ## LEMNA distance decay plots ####
     
     ## read and transform coordinates
     all_coordinates = read.csv("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/duckweed_coordinates.csv")
     all_coordinates$latitude = sapply(all_coordinates[,"GPS_S"], convert_dmm_to_dd)
     all_coordinates$longitute = sapply(all_coordinates[,"GPS_E"], convert_dmm_to_dd)
     
     ## assemble plot
     
     par(mfrow=c(1,2))
     
     ## full dataset
     
     ## subset to lemna
     lemna_coor = all_coordinates[which(colnames(lemna_hamdist) %in% all_coordinates[,"ID"]),]
     
     ## calculate geographic distance
     lemna_geodist = round(distm(lemna_coor[, c("longitute", "latitude")], fun = distVincentyEllipsoid),2)
     rownames(lemna_geodist) = lemna_coor[,"ID"]; colnames(lemna_geodist) = lemna_coor[,"ID"]
     
     ## transform to vector
     lemna_geodist_vector = as.vector(lemna_geodist)
     lemna_hamdist_vector = as.vector(lemna_hamdist)
     
     ## replace 0 values, for log transform
     lemna_geodist_vector[which(lemna_geodist_vector==0)] = 1
     
     ## assemble plot
     plot(log(lemna_geodist_vector), lemna_hamdist_vector,
          pch=21, bg="darkgreen", cex=0.8, ylim=c(-0.01,max(lemna_hamdist)),
          xlab="pairwise log geographic distance (m)", 
          ylab="pairwise genetic distance (Hamming)",
          main="lemna distance-decay plot")
     abline(lm(lemna_hamdist_vector~log(lemna_geodist_vector)), lwd=5, col="red")
     abline(v=log(10), lty=2); text(log(10), -0.01, label="10m")
     abline(v=log(100), lty=2); text(log(100), -0.01, label="100m")
     abline(v=log(1000), lty=2); text(log(1000), -0.01, label="1km")
     abline(v=log(10000), lty=2); text(log(10000), -0.01, label="10km")
     abline(v=log(100000), lty=2); text(log(100000), -0.01, label="100km")
     
     ## no clone dataset
     
     ## subset to lemna
     lemna_noclone_coor = all_coordinates[which(colnames(lemna_noclone_hamdist) %in% all_coordinates[,"ID"]),]
     
     ## calculate geographic distance
     lemna_noclone_geodist = round(distm(lemna_noclone_coor[, c("longitute", "latitude")], fun = distVincentyEllipsoid),2)
     rownames(lemna_noclone_geodist) = lemna_noclone_coor[,"ID"]; colnames(lemna_noclone_geodist) = lemna_noclone_coor[,"ID"]
     
     ## transform to vector
     lemna_noclone_geodist_vector = as.vector(lemna_noclone_geodist)
     lemna_noclone_hamdist_vector = as.vector(lemna_noclone_hamdist)
     
     ## replace 0 values, for log transform
     lemna_noclone_geodist_vector[which(lemna_noclone_geodist_vector==0)] = 1
     
     ## assemble plot
     plot(log(lemna_noclone_geodist_vector), lemna_noclone_hamdist_vector,
          pch=21, bg="darkgreen", cex=0.8, ylim=c(-0.01,max(lemna_noclone_hamdist)),
          xlab="pairwise log geographic distance (m)", 
          ylab="pairwise genetic distance (Hamming)",
          main="lemna no clone distance-decay plot")
     abline(v=log(10), lty=2); text(log(10), -0.01, label="10m")
     abline(v=log(100), lty=2); text(log(100), -0.01, label="100m")
     abline(v=log(1000), lty=2); text(log(1000), -0.01, label="1km")
     abline(v=log(10000), lty=2); text(log(10000), -0.01, label="10km")
     abline(v=log(100000), lty=2); text(log(100000), -0.01, label="100km")
     
     ## LANDOLTIA NO CLONE VERSION
     
     ## no clone dataset
     
     ## subset to landoltia
     landoltia_noclone_coor = all_coordinates[which(colnames(landoltia_noclone_hamdist) %in% all_coordinates[,"ID"]),]
     
     ## calculate geographic distance
     landoltia_noclone_geodist = round(distm(landoltia_noclone_coor[, c("longitute", "latitude")], fun = distVincentyEllipsoid),2)
     rownames(landoltia_noclone_geodist) = landoltia_noclone_coor[,"ID"]; colnames(landoltia_noclone_geodist) = landoltia_noclone_coor[,"ID"]
     
     ## transform to vector
     landoltia_noclone_geodist_vector = as.vector(landoltia_noclone_geodist)
     landoltia_noclone_hamdist_vector = as.vector(landoltia_noclone_hamdist)
     
     ## replace 0 values, for log transform
     landoltia_noclone_geodist_vector[which(landoltia_noclone_geodist_vector==0)] = 1
     
     ## assemble plot
     plot(log(landoltia_noclone_geodist_vector), landoltia_noclone_hamdist_vector,
          pch=21, bg="purple", cex=0.8, ylim=c(-0.01,max(landoltia_noclone_hamdist)),
          xlab="pairwise log geographic distance (m)", 
          ylab="pairwise genetic distance (Hamming)",
          main="Landoltia no clone distance-decay plot")
     abline(v=log(10), lty=2); text(log(10), -0.01, label="10m")
     abline(v=log(100), lty=2); text(log(100), -0.01, label="100m")
     abline(v=log(1000), lty=2); text(log(1000), -0.01, label="1km")
     abline(v=log(10000), lty=2); text(log(10000), -0.01, label="10km")
     abline(v=log(100000), lty=2); text(log(100000), -0.01, label="100km")
     
     
     
     