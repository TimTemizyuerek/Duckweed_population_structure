     ## load packages ####
     
     library(geosphere) # to compute geographical distances
     library(ape) # tree
     library(adegenet) #PCA
     library(LEA) ## for population structure analysis
     library(RColorBrewer) ## ready made color paletes
     library(vegan) ## for mantel()
     library(data.table) ## for set()
     library(plotrix) ## for floating.pie()
     library(sf) ## for reading shapefiles read_sf()
     library(stringr) ## for str output string handling
     
     # install_version("SNPfiltR", "1.0.1")
     library(snpStats) ## for linkage disequilibrium
     library(SNPfiltR)
     library(vcfR)
     library(dartR)
     library(SNPRelate) ## PCA
     library(ASRgenomics) ## PCA
     
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
       rownames(out_rarecurve) = "abundance"
       
       return(out_rarecurve)
     }
     
     ## helper function to calculate boot strap confidence intervals for diversity estimates. 
     diversity_stat <- function(data, indices) {
       resampled <- data[indices]
       counts_boot <- tabulate(resampled, nbins = length(unique(data)))
       vegan::diversity(counts_boot, index = "simpson")
     }
     
## FILTERING and DISTANCE CALCULATION ####
     ## initial filter with vcftools ####
     
     # ## Landoltia
     # ## put together VCFtools command
     # landoltia_vcftools_cmd = paste("vcftools",
     #                                "--vcf", shQuote("/mnt/c/WSL_storage/landoltia_populations.snps.vcf"),
     #                                "--recode",
     #                                "--recode-INFO-all",
     #                                "--mac 3",
     #                                "--minDP 5",
     #                                "--min-alleles 2",
     #                                "--max-alleles 2",
     #                                "--max-missing 0.5",
     #                                "--out", shQuote("landoltia_post_vcftools"))
     #      
     # # call command
     # system2("wsl", args = landoltia_vcftools_cmd)
     # 
     # ## Lemna
     # ## put together VCFtools command
     # lemna_vcftools_cmd = paste("vcftools",
     #                            "--vcf", shQuote("/mnt/c/WSL_storage/lemna_populations.snps.vcf"),
     #                            "--recode",
     #                            "--recode-INFO-all",
     #                            "--mac 3",
     #                            "--minDP 5",
     #                            "--min-alleles 2",
     #                            "--max-alleles 2",
     #                            "--max-missing 0.5",
     #                            "--out", shQuote("lemna_post_vcftools"))
     #                                
     # 
     # # call command
     # system2("wsl", args = lemna_vcftools_cmd)

     ## comprehensive filtering in R ####     

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
      
     ## remove linked loci in the same stack
     landoltia_postcvftools_qlty_ab_maxdp_smpl06_thin = distance_thin(landoltia_postcvftools_qlty_ab_maxdp_smpl06, min.distance = 400)
     lemna_postcvftools_qlty_ab_maxdp_smpl06_thin = distance_thin(lemna_postcvftools_qlty_ab_maxdp_smpl06, min.distance = 400)
     
     ## apply more stringent missing SNP filter
     landoltia_final = missing_by_snp(landoltia_postcvftools_qlty_ab_maxdp_smpl06_thin, cutoff = .95)
     lemna_final = missing_by_snp(lemna_postcvftools_qlty_ab_maxdp_smpl06_thin, cutoff = .95)
     
     ## export data in vcf format for dartR's gl.read.vcf
     # write.vcf(landoltia_final, "landoltia_final.vcf.gz")
     # write.vcf(lemna_final, "lemna_final.vcf.gz")
      
     ## LANDOLTIA compute hamming distance ####
     
     ## make genind object
     landoltia_final_genind = vcfR2genind(landoltia_final)
     
     ## COMPUTE HAMMING DISTANCE
     
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
     
     ## rename P1S2 to P2S1 to reflect that place 1 & 2 are actually connected
     colnames(landoltia_hamdist_raw)[1] = "P2S1"; rownames(landoltia_hamdist_raw)[1] = "P2S1"
     
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
     for (n in 1:length(landoltia_clone_list)) {landoltia_clone_names[n] = paste0("la_clone",n, sep="")}
     names(landoltia_clone_df) = landoltia_clone_names
     for (n in 1:length(landoltia_clone_list)) {
       
       ## transform to vector and put in matrix
       runner = landoltia_clone_list[[n]]
       runner = c(runner, rep(NA, max(lengths)-length(runner)))
       landoltia_clone_df[,n] = runner
       
     }
     
     ## add single individual clones
     landoltia_multi_clones = as.vector(na.omit(as.vector(as.matrix(landoltia_clone_df))))
     landoltia_single_clones = colnames(landoltia_hamdist)[which(!colnames(landoltia_hamdist) %in% landoltia_multi_clones)]
     for (n in 1:length(landoltia_single_clones)) {
       
       runner_colname = paste("la_clone", ncol(landoltia_clone_df)+1, sep="")
       
       landoltia_clone_df[[runner_colname]] = c(landoltia_single_clones[n],rep(NA, nrow(landoltia_clone_df)-1))
       
     }
     
     ## export table
     # write.csv(landoltia_clone_df, paste("Landoltia_clones_hamming0-",landoltia_clone_cutoff,".csv", sep=""))
     
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
     for (n in 1:length(lemna_clone_list)) {lemna_clone_names[n] = paste0("le_clone",n, sep="")}
     names(lemna_clone_df) = lemna_clone_names
     for (n in 1:length(lemna_clone_list)) {
       
       ## transform to vector and put in matrix
       runner = lemna_clone_list[[n]]
       runner = c(runner, rep(NA, max(lengths)-length(runner)))
       lemna_clone_df[,n] = runner
       
     }
     
     ## add single individual clones
     lemna_multi_clones = as.vector(na.omit(as.vector(as.matrix(lemna_clone_df))))
     lemna_single_clones = colnames(lemna_hamdist)[which(!colnames(lemna_hamdist) %in% lemna_multi_clones)]
     for (n in 1:length(lemna_single_clones)) {
       
       runner_colname = paste("le_clone", ncol(lemna_clone_df)+1, sep="")
       
       lemna_clone_df[[runner_colname]] = c(lemna_single_clones[n],rep(NA, nrow(lemna_clone_df)-1))
       
     }
     
     ## export table
     # write.csv(lemna_clone_df, paste("lemna_clones_hamming0-",lemna_clone_cutoff,".csv", sep=""))
     
     ## remove clones for STRUCTURE ####
     
     ## identify clones to keep
     landoltia_clones_to_remove = as.vector(na.omit(unlist(landoltia_clone_df[2:nrow(landoltia_clone_df),], use.names=FALSE)))
     lemna_clones_to_remove = as.vector(na.omit(unlist(lemna_clone_df[2:nrow(lemna_clone_df),], use.names=FALSE)))

     landoltia_clones_to_keep <- setdiff(colnames(landoltia_final@gt)[-1], landoltia_clones_to_remove)
     landoltia_clones_to_keep = c("FORMAT", landoltia_clones_to_keep)

     lemna_clones_to_keep <- setdiff(colnames(lemna_final@gt)[-1], lemna_clones_to_remove)
     lemna_clones_to_keep = c("FORMAT", lemna_clones_to_keep)

     ## remove from final vcfR object
     landoltia_final_no_clone = landoltia_final[, landoltia_clones_to_keep]
     lemna_final_no_clone = lemna_final[, lemna_clones_to_keep]

     ## create population file for PGDspider
     landoltia_names_vec = colnames(landoltia_final_no_clone@gt)[-1]
     landoltia_pop_vec = substr(landoltia_names_vec,1,3)
     landoltia_pop_file = data.frame(INDIVIDUAL = landoltia_names_vec, POP = landoltia_pop_vec)
     # write.table(landoltia_pop_file, file = "landoltia_population_for_PGDSpider.txt",
     #             sep = "\t", row.names = FALSE, quote = FALSE)

     lemna_names_vec = colnames(lemna_final_no_clone@gt)[-1]
     lemna_pop_vec = c("P30", "P32", "P11", "P14", "P36", "P27", "P36", "P30", "P23", "P05",
                       "P07", "P14", "P10", "P18", "P19", "P28", "P10", "P22", "P11", "P16")
     lemna_pop_file = data.frame(INDIVIDUAL = lemna_names_vec, POP = lemna_pop_vec)
     # write.table(lemna_pop_file, file = "lemna_population_for_PGDSpider.txt",
     #             sep = "\t", row.names = FALSE, quote = FALSE)

     ## export data in vcf format for spider to make structure files
     # write.vcf(landoltia_final_no_clone, "landoltia_final_no_clone.vcf.gz")
     # write.vcf(lemna_final_no_clone, "lemna_final_no_clone.vcf.gz")

     
## DATA VISUALISATION / RESULTS ####
     ## colour set-up ####
     
     lemna_col = "#4DAF4A"
     landoltia_col = "#984EA3"
     both_col = "#999999"
     
     P10_col = "#E41A1C"
     P14_col = "#377EB8"
     P19_col = "#A65628" 
     P27_col = "#F781BF"
     P36_col = "#FF7F00"
     rest_col = "#FFFF33"
     
     ## Numbers for M&M ####
     
     ## population level diversity
     diversity(hamdist_to_rarecurve(landoltia_hamdist))
     diversity(hamdist_to_rarecurve(lemna_hamdist))
     
     mean(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)])
     mean(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)])
     
     sd(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)])
     sd(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)])
     
     ## Observed heterozygosity, expected heretozygostiy and inbreeding coefficient (FIS)
     
     ## LANDOLTIA
         
         ## extract SNP information 
         landoltia_genotype_matrix = extract.gt(landoltia_final, element = "GT", as.numeric=FALSE)
         
         ## recode to numeric
         landoltia_genotype_matrix_num = matrix(NA, nrow = nrow(landoltia_genotype_matrix), ncol = ncol(landoltia_genotype_matrix))
         rownames(landoltia_genotype_matrix_num) = rownames(landoltia_genotype_matrix)
         colnames(landoltia_genotype_matrix_num) = colnames(landoltia_genotype_matrix)
         landoltia_genotype_matrix_num[landoltia_genotype_matrix %in% c("0/0")] = 0
         landoltia_genotype_matrix_num[landoltia_genotype_matrix %in% c("0/1","1/0")] = 1
         landoltia_genotype_matrix_num[landoltia_genotype_matrix %in% c("1/1")] = 2
         landoltia_genotype_matrix_num[landoltia_genotype_matrix %in% c("./.")] = NA
         
         ## compute observed heterozygosity
         landoltia_obs_he_per_snp = rowSums(landoltia_genotype_matrix_num == 1, na.rm = TRUE) / rowSums(!is.na(landoltia_genotype_matrix_num))
         
         ## compute expected heterozygosity
         landoltia_exp_he_per_snp = apply(landoltia_genotype_matrix_num, 1, function(snp) {
           ## remove NA
           snp = snp[!is.na(snp)]
           ## count alleles
           n_ref = sum(2*(snp == 0) + 1*(snp == 1))
           n_alt = sum(2*(snp == 2) + 1*(snp == 1))
           ## compute allele frequencies
           p = n_ref / (n_ref + n_alt)
           q = n_alt / (n_ref + n_alt)
           ## compute expected hetereozygosity
           1 - (p^2 + q^2)
         })
         
         ## compute FIS (inbreeding coefficient)
         landoltia_mean_FIS_per_snp = mean(1 - (landoltia_obs_he_per_snp / landoltia_exp_he_per_snp))
    
     ## LEMNA
     
         ## extract SNP information 
         lemna_genotype_matrix = extract.gt(lemna_final, element = "GT", as.numeric=FALSE)
         
         ## recode to numeric
         lemna_genotype_matrix_num = matrix(NA, nrow = nrow(lemna_genotype_matrix), ncol = ncol(lemna_genotype_matrix))
         rownames(lemna_genotype_matrix_num) = rownames(lemna_genotype_matrix)
         colnames(lemna_genotype_matrix_num) = colnames(lemna_genotype_matrix)
         lemna_genotype_matrix_num[lemna_genotype_matrix %in% c("0/0")] = 0
         lemna_genotype_matrix_num[lemna_genotype_matrix %in% c("0/1","1/0")] = 1
         lemna_genotype_matrix_num[lemna_genotype_matrix %in% c("1/1")] = 2
         lemna_genotype_matrix_num[lemna_genotype_matrix %in% c("./.")] = NA
         
         ## compute observed heterozygosity
         lemna_obs_he_per_snp = rowSums(lemna_genotype_matrix_num == 1, na.rm = TRUE) / rowSums(!is.na(lemna_genotype_matrix_num))
         
         ## compute expected heterozygosity
         lemna_exp_he_per_snp = apply(lemna_genotype_matrix_num, 1, function(snp) {
           ## remove NA
           snp = snp[!is.na(snp)]
           ## count alleles
           n_ref = sum(2*(snp == 0) + 1*(snp == 1))
           n_alt = sum(2*(snp == 2) + 1*(snp == 1))
           ## compute allele frequencies
           p = n_ref / (n_ref + n_alt)
           q = n_alt / (n_ref + n_alt)
           ## compute expected hetereozygosity
           1 - (p^2 + q^2)
         })
         
         ## compute FIS (inbreeding coefficient)
         lemna_mean_FIS_per_snp = mean(1 - (lemna_obs_he_per_snp / lemna_exp_he_per_snp))
         
     ## waterbody df
     shared_waterbodies = data.frame(P1 = "landoltia",
                                     P2 = "landoltia",
                                     P4 = "landoltia",
                                     P5 = "lemna",
                                     P6 = "both",
                                     P7 = "lemna",
                                     P10 = "both",
                                     P11 = "both",
                                     P12 = "landoltia",
                                     P13 = "landoltia",
                                     P14 = "both",
                                     P15 = "both",
                                     P16 = "both",
                                     P17 = "both",
                                     P18 = "both",
                                     P19 = "both",
                                     P20 = "landoltia",
                                     P21 = "lemna",
                                     P22 = "both",
                                     P23 = "both",
                                     P24 = "landoltia",
                                     P25 = "landoltia",
                                     P26 = "landoltia",
                                     P27 = "both",
                                     P28 = "both",
                                     P30 = "lemna",
                                     P31 = "lemna",
                                     P32 = "both",
                                     P33 = "lemna",
                                     P34 = "both",
                                     P35 = "lemna",
                                     P36 = "both")
     
     length(which(substr(colnames(landoltia_hamdist),1,3) %in% "P10"))
     length(which(substr(colnames(landoltia_hamdist),1,3) %in% "P14"))
     length(which(substr(colnames(landoltia_hamdist),1,3) %in% "P19"))
     length(which(substr(colnames(landoltia_hamdist),1,3) %in% "P27"))
     length(which(substr(colnames(landoltia_hamdist),1,3) %in% "P36"))
     
     length(which(substr(colnames(lemna_hamdist),1,3) %in% "P10"))
     length(which(substr(colnames(lemna_hamdist),1,3) %in% "P14"))
     length(which(substr(colnames(lemna_hamdist),1,3) %in% "P19"))
     length(which(substr(colnames(lemna_hamdist),1,3) %in% "P27"))
     length(which(substr(colnames(lemna_hamdist),1,3) %in% "P36"))
     
     library(readxl)
     sampling_data = read_xlsx("C:/Users/timte/Desktop/Brisbane/Chapter 1/Duckweed collection 30.12.2024 - for material methods.xlsx")
     sampling_data = as.data.frame(sampling_data)
     
     range(as.numeric(sampling_data[,"DNA_fronds"]), na.rm=TRUE)
     
     sort(unique(sampling_data[,"ID"]))
     nrow(sampling_data)
     
     ## MAP ####
     
     ## number of waterbodies
     length(unique(c(substr(colnames(landoltia_hamdist),1,3), substr(colnames(lemna_hamdist),1,3))))
     
     ## read and transform coordinates
     all_coordinates = read.csv("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/duckweed_coordinates.csv")
     all_coordinates$latitude = sapply(all_coordinates[,"GPS_S"], convert_dmm_to_dd)
     all_coordinates$longitude = sapply(all_coordinates[,"GPS_E"], convert_dmm_to_dd)
     
     ## read micro sites data
     micro_sites = read.csv("C:/Users/timte/Desktop/Brisbane/Chapter 1/micro_sites.csv", sep=";")
     
     ## extract summary dataframe
     micro_summary = data.frame("micro_site_ID" = character(), "lemna" = integer(), "landoltia" = integer())
     for (n in unique(micro_sites[,"micro_site_ID"])) {
       micro_runner = micro_sites[which(micro_sites[,"micro_site_ID"] == n),, drop=FALSE]
       filler_row = data.frame("micro_site_ID" = n,
                               "lemna" = sum(micro_runner[,"species"] == "lemna"), 
                               "landoltia" = sum(micro_runner[,"species"] == "landoltia")) 
       micro_summary = rbind(micro_summary, filler_row)
     }
     micro_summary[,"sum"] = apply(micro_summary[,2:3], 1, sum)
     
     ## combine clone datasets
     landoltia_clone_df[19:27,] = matrix(NA, ncol=ncol(landoltia_clone_df), nrow=9)
     total_clone_df = cbind(lemna_clone_df,landoltia_clone_df)
     
     ## add clone col to micro_sites 
     clone_vec = rep(NA, nrow(micro_sites))
     for (n in colnames(total_clone_df)) {
       
       runner_col = total_clone_df[,n]
       clone_vec[which(micro_sites[,"samples"] %in% runner_col)] = n
       
       
     }
     micro_sites[,"clone"] = clone_vec
     
     ## add coordinates to each sample
     coor_df = data.frame(matrix(ncol=3, nrow=0))
     for (n in unique(micro_sites[,"samples"])) {
       
       coor_df = rbind(coor_df, all_coordinates[which(all_coordinates[,"ID"] == n),c("ID","latitude", "longitude")])
       
     }
     names(coor_df)[1] = "samples"
     micro_sites = merge(micro_sites, coor_df, by="samples")
     
     ## load shapefiles
     sf_use_s2(FALSE)
     brisbane_waterbodies = read_sf("C:/Users/timte/Desktop/Brisbane/Chapter 1/Brisbane wetland areas/Wetland_areas.shp")
     brisbane_waterbodies = st_make_valid(brisbane_waterbodies)
     brisbane_waterbodies = st_crop(brisbane_waterbodies, c(xmin = min(all_coordinates[,"longitude"])-0.01,
                                                            xmax = max(all_coordinates[,"longitude"])+0.05,
                                                            ymin = -max(all_coordinates[,"latitude"])-0.05,
                                                            ymax = -min(all_coordinates[,"latitude"])+0.05))
                                            
     brisbane_coastline = read_sf("C:/Users/timte/Desktop/Brisbane/Chapter 1/Brisbane coastline/Coastline.shp")
     brisbane_coastline = st_crop(brisbane_coastline, c(xmin = min(all_coordinates[,"longitude"])-0.01,
                                                        xmax = max(all_coordinates[,"longitude"])+0.05,
                                                        ymin = -max(all_coordinates[,"latitude"])-0.05,
                                                        ymax = -min(all_coordinates[,"latitude"])+0.05))
     
     ## align shapefiles
     brisbane_waterbodies = st_transform(brisbane_waterbodies, st_crs(brisbane_coastline))
     
     ## buffer around coast
     sf_use_s2(TRUE)
     coast_buffer = st_buffer(brisbane_coastline, dist = 5)
     
     ## separate land and sea
     is_saltwater = lengths(st_intersects(brisbane_waterbodies, coast_buffer)) > 0
     brisbane_saltwater = brisbane_waterbodies[is_saltwater, ]
     brisbane_freshwater = brisbane_waterbodies[!is_saltwater, ]
     
     ## set up plotting area
     
         split.screen(rbind(c(0,1,0,1),
                            c(0.09,0.35,0.58,0.84)))  
         
     ## plot background map
     
         ## these are criminal adjustments, but I seem to need them ... 
         xmin = min(all_coordinates[,"longitude"])-0.01
         xmax = max(all_coordinates[,"longitude"])+0.05
         ymin = -max(all_coordinates[,"latitude"])-0.05
         ymax = -min(all_coordinates[,"latitude"])+0.05
         
         ## assemble plot
         screen(1)
         par(mar=c(2,2.2,0.1,0.1))
         plot(st_geometry(brisbane_freshwater), col = "dodgerblue", border = NA)
         plot(st_geometry(brisbane_saltwater), col = "white", border = NA, add = TRUE)
         plot(st_geometry(brisbane_coastline), col = "black", lwd = 1, add = TRUE)
         rect(xmin, ymin, xmax, ymax,
              col=scales::alpha("white", 0.25), border="black", lwd=1)
         ## compute scalebar
         #distm(rbind(c(152.36,-27.74),c(152.614, -27.74)), fun = distVincentyEllipsoid)
         lines(x=c(152.08,152.324), y=c(-27.735,-27.735), lwd=2)
         text(152.202, -27.72, labels="25km")
         ## north pointer
         points(152.065,-27.65,pch=17, cex=1.8)
         lines(c(152.065,152.065), c(-27.65,-27.7),lwd=2)
         text(152.065, -27.71, labels="N")
         legend(153.15, -27.138, legend=c("Lemna", "Landoltia"),
                fill=c(lemna_col, landoltia_col), )
         
         ## manual x-axis
         xticks = round(seq(from=xmin, to=xmax-0.02, length.out = 5),2)
         segments(xticks, ymin, xticks, ymin-(ymax-ymin)*0.01)
         text(xticks, ymin - (ymax-ymin)*0.05, labels = xticks)
         
         ## manual y-axis 
         yticks = round(seq(from=ymin, to=ymax, length.out = 5),2)
         segments(xmin, yticks, xmin-(xmax-xmin)*0.005, yticks)
         mtext(side = 2, at = yticks, text = yticks, las = 1, line = -0.5)
         
         ## extract waterbody data
         waterbody_map = data.frame(matrix(0, ncol=6,nrow = 0))
         m=0; for (n in unique(substr(micro_sites[,"micro_site_ID"],1,3))){
           
           m=m+1
           
           ## extract average coordinates
           runner_waterbody_site = micro_sites[which(substr(micro_sites[,"micro_site_ID"],1,3) == n),]
           runner_coor = all_coordinates[which(all_coordinates[,"ID"] %in% runner_waterbody_site[,"samples"]),]
           
           waterbody_map[m,1] = n
           waterbody_map[m,2] = mean(runner_coor[,"latitude"])
           waterbody_map[m,3] = mean(runner_coor[,"longitude"])
           waterbody_map[m,4] = length(which(substr(runner_coor[,"species"],1,2) == "Le"))
           waterbody_map[m,5] = length(which(substr(runner_coor[,"species"],1,2) == "La"))
           waterbody_map[m,6] = waterbody_map[m,4] + waterbody_map[m,5]
           
         }
         colnames(waterbody_map) = c("micro","long", "lat", "lemna", "landoltia", "total")    
         waterbody_map[,"long"] = as.numeric(waterbody_map[,"long"]);waterbody_map[,"lat"] = as.numeric(waterbody_map[,"lat"])
         waterbody_map[,"lemna"] = as.numeric(waterbody_map[,"lemna"]);waterbody_map[,"landoltia"] = as.numeric(waterbody_map[,"landoltia"])
         waterbody_map[,"total"] = as.numeric(waterbody_map[,"total"])
         
         ## transform total for scaling in plot
         waterbody_map$scaled_total = 0.01 + (0.04*(waterbody_map[,"total"]-1))/51
         
         ## reorder so that the large piecharts are plotted first
         waterbody_map = waterbody_map[c(1,5,10,19,28,
                                         2,3,4,
                                         6,7,8,9,
                                         11,12,13,14,15,16,17,18,
                                         20,21,22,23,24,25,26,27,
                                         29,30,31,32),]
         ## add piecharts
         for (n in 1:nrow(waterbody_map)) {
           floating.pie(xpos=waterbody_map[,"lat"][n], ypos=-waterbody_map[,"long"][n], 
                        x=c(waterbody_map[,4][n], waterbody_map[,5][n]), radius=waterbody_map[,"scaled_total"][n],
                        col=c(landoltia_col, lemna_col),
                        edges=1000)
         }
         
         ## plot letter to refernce inlet
         text(waterbody_map[1,]["lat"],-waterbody_map[1,]["long"], labels="P10")
         close.screen(1)
         
     ## plot P10 inlet
     
         ## find micro sites
         P10_microsite_plotter = data.frame(matrix(0, ncol=6,nrow = 0))
         P10_microsites = micro_sites[which(substr(micro_sites[,"samples"],1,3) == "P10"),]
         m=0; for (n in unique(P10_microsites[,"micro_site_ID"])){
           
           m=m+1
           
           ## extract average coordinates
           runner_waterbody_site = P10_microsites[which(P10_microsites[,"micro_site_ID"] == n),]
           runner_coor = all_coordinates[which(all_coordinates[,"ID"] %in% runner_waterbody_site[,"samples"]),]
           
           P10_microsite_plotter[m,1] = n
           P10_microsite_plotter[m,2] = mean(runner_coor[,"latitude"])
           P10_microsite_plotter[m,3] = mean(runner_coor[,"longitude"])
           P10_microsite_plotter[m,4] = length(which(substr(runner_coor[,"species"],1,2) == "Le"))
           P10_microsite_plotter[m,5] = length(which(substr(runner_coor[,"species"],1,2) == "La"))
           P10_microsite_plotter[m,6] = P10_microsite_plotter[m,4] + P10_microsite_plotter[m,5]
           
         }
         colnames(P10_microsite_plotter) = c("micro","long", "lat", "lemna", "landoltia", "total")    
         P10_microsite_plotter[,"long"] = as.numeric(P10_microsite_plotter[,"long"]);P10_microsite_plotter[,"lat"] = as.numeric(P10_microsite_plotter[,"lat"])
         P10_microsite_plotter[,"lemna"] = as.numeric(P10_microsite_plotter[,"lemna"]);P10_microsite_plotter[,"landoltia"] = as.numeric(P10_microsite_plotter[,"landoltia"])
         P10_microsite_plotter[,"total"] = as.numeric(P10_microsite_plotter[,"total"])
         
         ## extract P10 samples
         P10_coordinates = all_coordinates[which(substr(all_coordinates[,"ID"],1,3) == "P10"),]
         
         ## cut shapefile
         P10_shapefile = st_crop(brisbane_waterbodies, c(xmin = min(P10_coordinates[,"longitude"]),
                                                         xmax = max(P10_coordinates[,"longitude"]),
                                                         ymin = -max(P10_coordinates[,"latitude"]),
                                                         ymax = -min(P10_coordinates[,"latitude"])))
         
         # extend a little to the right for plotting piecharts
         bb <- st_bbox(P10_shapefile)
         xlim <- c(bb["xmin"]-0.00008, bb["xmax"]+0.00008)
         ylim <- c(bb["ymin"]-0.00008, bb["ymax"]+0.00008)
         
         ## assemble plot
         screen(2)
         par(mar=c(0,0,0,0))
         plot(st_geometry(P10_shapefile), col = "dodgerblue", border = NA,
              xlim = xlim, ylim = ylim)
         rect(bb["xmin"]-0.00006, bb["ymin"]-0.0001, bb["xmax"]+0.00008, bb["ymax"]+0.00006,border="black", lwd=1)
         text(153.0644,-27.54035, labels="P10")
         
         ## compute scalebar
         distm(rbind(c(153.063554,-27.74),c(153.0633, -27.74)), fun = distVincentyEllipsoid)
         lines(x=c(153.0634,153.0637),y=c(-27.54105, -27.54105), lwd=2)
         text(153.06355, -27.54099, labels="25m")
         
         ## add piecharts
         for (n in 1:nrow(P10_microsite_plotter)) {
           floating.pie(xpos=P10_microsite_plotter[,"lat"][n], ypos=-P10_microsite_plotter[,"long"][n], 
                        x=c(P10_microsite_plotter[,4][n], P10_microsite_plotter[,5][n]), col=c(landoltia_col, lemna_col),
                        edges=1000, radius = 0.00003)
     }
         close.screen(2)
         close.screen(all.screens = TRUE)
         
     ## MICRO SITE: competitive environment ####
     
     # make categories for stacked barplot
     micro_summary[,"group"] = ifelse(micro_summary[,"lemna"] > 0 & micro_summary[,"landoltia"] == 0, "lemna",
                                      ifelse(micro_summary[,"landoltia"] > 0 & micro_summary[,"lemna"] == 0, "landoltia",
                                             ifelse(micro_summary[,"lemna"] > 0 & micro_summary[,"landoltia"] > 0, "both", NA)))
    
     ## get abundances per species
     lemna_counts = table(factor(micro_summary[,"sum"][micro_summary[,"group"] == "lemna"], levels = 1:6))
     landoltia_counts = table(factor(micro_summary[,"sum"][micro_summary[,"group"] == "landoltia"], levels = 1:6))
     both_counts = table(factor(micro_summary[,"sum"][micro_summary[,"group"] == "both"], levels = 1:6))
     
     ## assembly plotter dataframe
     stacked_barplotter = rbind(lemna_counts, landoltia_counts, both_counts)
     
     ## two species within micro site diversity
    
     ## calculate average within site distances for competitive environment plot
     micro_distance = data.frame("micro_site_ID" = character(), 
                                 "lemna_avg_distance" = integer(), "lemna_sd_distance" = integer(),
                                 "landoltia_avg_distanc" = integer(), "landoltia_sd_distance" = integer())
     for (n in unique(micro_sites[,"micro_site_ID"])) {
       
       micro_runner = micro_sites[which(micro_sites[,"micro_site_ID"] == n),, drop=FALSE]
       
       ## extract subset from lemna_hamdist from micro site
       runner_lemna_dist = lemna_hamdist[which(colnames(lemna_hamdist) %in% micro_runner[which(micro_runner[,"species"] == "lemna"),][,"samples"]),which(colnames(lemna_hamdist) %in% micro_runner[which(micro_runner[,"species"] == "lemna"),][,"samples"])]
       
       ## extract subset from landoltia_hamdist from micro site
       runner_landoltia_dist = landoltia_hamdist[which(colnames(landoltia_hamdist) %in% micro_runner[which(micro_runner[,"species"] == "landoltia"),][,"samples"]),which(colnames(landoltia_hamdist) %in% micro_runner[which(micro_runner[,"species"] == "landoltia"),][,"samples"])]
       
       ## compute mean
       filler_row = data.frame("micro_site_ID" = n,
                               "lemna_avg_distance" = if(length(runner_lemna_dist[lower.tri(runner_lemna_dist, diag=FALSE)]) == 0) NA else mean(runner_lemna_dist[lower.tri(runner_lemna_dist, diag=FALSE)]),
                               "lemna_sd_distance" = sd(runner_lemna_dist[lower.tri(runner_lemna_dist, diag=FALSE)]),
                               "landoltia_avg_distance" = if(length(runner_landoltia_dist[lower.tri(runner_landoltia_dist, diag=FALSE)]) == 0) NA else mean(runner_landoltia_dist[lower.tri(runner_landoltia_dist, diag=FALSE)]),
                               "landoltia_sd_distance" = sd(runner_landoltia_dist[lower.tri(runner_landoltia_dist, diag=FALSE)])) 
       micro_distance = rbind(micro_distance, filler_row)
     }
     
     ## merge to micro_summary and retain only those with 2 or more individuals
     micro_scatter_plotter_full = merge(micro_summary, micro_distance, by="micro_site_ID")
     micro_scatter_plotter = micro_scatter_plotter_full[which(micro_scatter_plotter_full[,"lemna"] >= 2 & micro_scatter_plotter_full[,"landoltia"] >= 2),]
         
     ## permutations for histogram
     landoltia_perm_mean = vector(); lemna_perm_mean = vector()
     for (n in 1:100000) {
       sim_size = sample(2:3,1)
       runner_landoltia_sample = sample(colnames(landoltia_hamdist), sim_size)
       runner_lemna_sample = sample(colnames(lemna_hamdist), sim_size)
       runner_landoltia_matrix = landoltia_hamdist[runner_landoltia_sample,runner_landoltia_sample]
       runner_lemna_matrix = lemna_hamdist[runner_lemna_sample,runner_lemna_sample]
       landoltia_perm_mean[n] = mean(runner_landoltia_matrix[lower.tri(runner_landoltia_matrix, diag=FALSE)])
       lemna_perm_mean[n] = mean(runner_lemna_matrix[lower.tri(runner_lemna_matrix, diag=FALSE)])
     }
     
     ## set up screens
     
         split.screen(rbind(c(0, 0.8, 0, 0.8),     ## scatterplot 
                            c(0.82, 1, 0, 0.8),    ## y-axis (landoltia)
                            c(0, 0.8, 0.82, 1)))   ## x-axis (lemna)
         
          
     ## scatterplot
         
         screen(1)
         par(mar=c(4.1,4.1,0,0))
         plot(NULL, xlab=expression("genetic distance within micro site (" * italic("L. aequinoctialis" * ")")),
              ylab=expression("genetic distance within micro site (" * italic("L. punctata" * ")")),
              xlim=c(0,max(lemna_perm_mean)), ylim=c(0,max(landoltia_perm_mean)))
         
         ## x=y under the points
         abline(a = 0, b = 1, lty=2, lwd=2, col="gray50")
         
         ## micro points
         points(micro_scatter_plotter[,"lemna_avg_distance"], micro_scatter_plotter[,"landoltia_avg_distance"],
                pch= 21,
                bg=ifelse(substr(micro_scatter_plotter[,"micro_site_ID"],1,3) == "P10", P10_col,
                          ifelse(substr(micro_scatter_plotter[,"micro_site_ID"],1,3) == "P14", P14_col,
                                 ifelse(substr(micro_scatter_plotter[,"micro_site_ID"],1,3) == "P19", P19_col,
                                        ifelse(substr(micro_scatter_plotter[,"micro_site_ID"],1,3) == "P27", P27_col,
                                               ifelse(substr(micro_scatter_plotter[,"micro_site_ID"],1,3) == "P36", P36_col, rest_col))))),
                cex=2)
         
         
         ## legend on top of the points/line                                              
         legend("topright", legend=c("P14", "P19", "P27", "P36", "other"),
                pt.bg=c(P14_col, P19_col, P27_col, P36_col, rest_col),
                pch=21, bty="n", pt.cex = 2)
         close.screen(1)
     
     ## add landoltia histogram
         
         screen(2)
         par(mgp=c(1,0,0))
         landoltia_breaks = pretty(range(landoltia_perm_mean), n=20)
         landoltia_micro_hist = hist(micro_scatter_plotter[,"landoltia_avg_distance"], breaks=landoltia_breaks, plot=FALSE)
         landoltia_global_hist = hist(landoltia_perm_mean, breaks=landoltia_breaks, plot=FALSE)
         ## assemble plot
         par(mar=c(4,0,0,0.5))
         plot(0, 0, type = "n", xlim = c(0, max(landoltia_micro_hist$density, landoltia_global_hist$density)), ylim = c(0,max(landoltia_perm_mean)), axes = FALSE, xlab = "Density", ylab = "")
         ## plot global
         for (n in seq_along(landoltia_global_hist$density)) {rect(0, landoltia_global_hist$breaks[n], landoltia_global_hist$density[n], landoltia_global_hist$breaks[n+1], col = both_col, border = "black")}
         ## plot micro
         for (n in seq_along(landoltia_micro_hist$density)) {rect(0, landoltia_micro_hist$breaks[n], landoltia_micro_hist$density[n], landoltia_micro_hist$breaks[n+1], col = scales::alpha(landoltia_col,0.75), border = "black")}
         box()
         legend("topright", col=c(landoltia_col, both_col),legend=c("micro", "global"), pch=15, bty="n")
         close.screen(2)
         
     ## add lemna histogram
         
         screen(3)
         par(mgp=c(1,0,0))
         lemna_breaks = pretty(range(lemna_perm_mean), n=20)
         lemna_micro_hist = hist(micro_scatter_plotter[,"lemna_avg_distance"], breaks=lemna_breaks, plot=FALSE)
         lemna_global_hist = hist(lemna_perm_mean, breaks=lemna_breaks, plot=FALSE)
         ## assemble plot
         par(mar=c(0,4,0.5,0))
         plot(0, 0, type = "n", ylim = c(0, max(lemna_micro_hist$density, lemna_global_hist$density)), xlim = c(0,max(lemna_perm_mean)), axes = FALSE, xlab = "", ylab = "Density")
         ## plot global
         for (n in seq_along(lemna_global_hist$density)) {rect(lemna_global_hist$breaks[n],0, lemna_global_hist$breaks[n+1], lemna_global_hist$density[n], ,col = both_col, border = "black")}
         ## add outlier point
         points(0.49,2, cex=0.8, pch=21, bg=both_col)
         ## plot micro
         for (n in seq_along(lemna_micro_hist$density)) {rect(lemna_micro_hist$breaks[n],0, lemna_micro_hist$breaks[n+1], lemna_micro_hist$density[n], ,col = scales::alpha(lemna_col,0.75), border = "black")}
         box()
         legend("topright", col=c(lemna_col, both_col),legend=c("micro", "global"), pch=15, bty="n")
         close.screen(3)
         
         ## wrap it up
         close.screen(all.screens = TRUE)
        
     ## t test to comparing genetic distances at micro sites
     t.test(micro_scatter_plotter[,"lemna_avg_distance"], micro_scatter_plotter[,"landoltia_avg_distance"], paired=FALSE)
     
     ## WATERBODY: within vs outside distance & diversity ####
     
     ## number of iterations
     iter = 100000
     boot_iter = 100000
     
     ## set screens
     split.screen(rbind(c(0, 0.55, 0.5, 1), 
                        c(0, 0.55, 0, 0.5),
                        c(0.55, 1, 0.5, 1), 
                        c(0.55, 1, 0, 0.5)))
     
     ## LANDOLTIA
     
         ## subset the ponds   
         landoltia_p10 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P10"),which(substr(colnames(landoltia_hamdist),1,3) == "P10")]
         landoltia_p14 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P14"),which(substr(colnames(landoltia_hamdist),1,3) == "P14")]
         landoltia_p19 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P19"),which(substr(colnames(landoltia_hamdist),1,3) == "P19")]
         landoltia_p27 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P27"),which(substr(colnames(landoltia_hamdist),1,3) == "P27")]
         landoltia_p36 = landoltia_hamdist[which(substr(colnames(landoltia_hamdist),1,3) == "P36"),which(substr(colnames(landoltia_hamdist),1,3) == "P36")]
         
         ## calculate permutations
         landoltia_perm_mean = vector()
         landoltia_perm_clone = vector()
         landoltia_perm_diversity = vector()
         for (n in 1:iter) {
           
           sim_size = sample(20:27,1)
           runner_sample = sample(colnames(landoltia_hamdist), sim_size)
           runner_matrix = landoltia_hamdist[runner_sample,runner_sample]
           landoltia_perm_mean[n] = mean(runner_matrix[lower.tri(runner_matrix, diag=FALSE)])
           landoltia_perm_clone[n] = ncol(hamdist_to_rarecurve(runner_matrix))
           landoltia_perm_diversity[n] = diversity(hamdist_to_rarecurve(runner_matrix), index = "simpson")
         }
         
         ## mean plot
         par(mar=c(0.3,4.5,0.3,0.3))
         screen(1)
         ## base plot
         plot(NULL, xlim=c(0.5,7.5), ylim=c(0,0.22), xlab="", main="", xaxt = "n", ylab="Genetic distance", las=2)
         ## mean + confidence interval of population
         boot_means <- replicate(boot_iter, mean(sample(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 7, y0 = ci[1], x1 = 7, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(7, mean(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)]), pch=21, bg=landoltia_col, cex=1.5)
         abline(v=6.5)
         ## means + confidence intervals of the deep sampled ponds
         boot_means <- replicate(boot_iter, mean(sample(landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 1, y0 = ci[1], x1 = 1, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(1, mean(landoltia_p10[lower.tri(landoltia_p10, diag=FALSE)]), pch=21, bg=landoltia_col, cex=1.5)
         boot_means <- replicate(boot_iter, mean(sample(landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 2, y0 = ci[1], x1 = 2, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(2, mean(landoltia_p14[lower.tri(landoltia_p14, diag=FALSE)]), pch=21, bg=landoltia_col, cex=1.5)
         boot_means <- replicate(boot_iter, mean(sample(landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 3, y0 = ci[1], x1 = 3, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(3, mean(landoltia_p19[lower.tri(landoltia_p19, diag=FALSE)]), pch=21, bg=landoltia_col, cex=1.5)
         boot_means <- replicate(boot_iter, mean(sample(landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 4, y0 = ci[1], x1 = 4, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(4, mean(landoltia_p27[lower.tri(landoltia_p27, diag=FALSE)]), pch=21, bg=landoltia_col, cex=1.5)
         boot_means <- replicate(boot_iter, mean(sample(landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 5, y0 = ci[1], x1 = 5, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(5, mean(landoltia_p36[lower.tri(landoltia_p36, diag=FALSE)]), pch=21, bg=landoltia_col, cex=1.5)
         ## permutation results
         clip(0,6,0,0.3)
         abline(h=quantile(landoltia_perm_mean, probs = c(0.025, 0.975)), lty=2, lwd=1)
         clip(0,6.5,0,0.3)
         points(rep(6,iter), landoltia_perm_mean, pch=21, bg=landoltia_col, cex=1.5)
         close.screen(1)
         
         ## diversity plot
         par(mar=c(3,4.5,0.3,0.3))
         screen(2)
         ## base plot
         plot(NULL, xlim=c(0.5,7.5), ylim=c(0,1.05), xlab="", main="", xaxt = "n", ylab="Diversity (Simpson)", las=2)
         axis(1, at = 1:7, las=2, labels = c("P10", "P14", "P19", "P27", "P36", "perm", "pop"))
         ## plot diversity and CI for population
         data_individuals = rep(seq_along(hamdist_to_rarecurve(landoltia_hamdist)), hamdist_to_rarecurve(landoltia_hamdist))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 7, y0 = ci[1], x1 = 7, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(7, diversity(hamdist_to_rarecurve(landoltia_hamdist), index = "simpson"), pch=21, bg=landoltia_col, cex=1.5)
         abline(v=6.5)
         ## plot diversity and CI for deep sampled ponds
         data_individuals = rep(seq_along(hamdist_to_rarecurve(landoltia_p10)), hamdist_to_rarecurve(landoltia_p10))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 1, y0 = ci[1], x1 = 1, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(1, diversity(hamdist_to_rarecurve(landoltia_p10), index = "simpson"), pch=21, bg=landoltia_col, cex=1.5)
         data_individuals = rep(seq_along(hamdist_to_rarecurve(landoltia_p14)), hamdist_to_rarecurve(landoltia_p14))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 2, y0 = ci[1], x1 = 2, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(2, diversity(hamdist_to_rarecurve(landoltia_p14), index = "simpson"), pch=21, bg=landoltia_col, cex=1.5)
         data_individuals = rep(seq_along(hamdist_to_rarecurve(landoltia_p19)), hamdist_to_rarecurve(landoltia_p19))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 3, y0 = ci[1], x1 = 3, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(3, diversity(hamdist_to_rarecurve(landoltia_p19), index = "simpson"), pch=21, bg=landoltia_col, cex=1.5)
         data_individuals = rep(seq_along(hamdist_to_rarecurve(landoltia_p27)), hamdist_to_rarecurve(landoltia_p27))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 4, y0 = ci[1], x1 = 4, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(4, diversity(hamdist_to_rarecurve(landoltia_p27), index = "simpson"), pch=21, bg=landoltia_col, cex=1.5)
         data_individuals = rep(seq_along(hamdist_to_rarecurve(landoltia_p36)), hamdist_to_rarecurve(landoltia_p36))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 5, y0 = ci[1], x1 = 5, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(5, diversity(hamdist_to_rarecurve(landoltia_p36), index = "simpson"), pch=21, bg=landoltia_col, cex=1.5)
         ## permutations
         clip(0,6,0,4)
         abline(h=quantile(landoltia_perm_diversity, probs = c(0.001, 0.999)), lty=2, lwd=1)
         clip(0,6.5,0,4)
         points(rep(6,iter), landoltia_perm_diversity, pch=21, bg=landoltia_col, cex=1.5)
         close.screen(2)
         
     ## LEMNA
     
         ## subset the ponds   
         lemna_p10 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P10"),which(substr(colnames(lemna_hamdist),1,3) == "P10")]
         lemna_p14 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P14"),which(substr(colnames(lemna_hamdist),1,3) == "P14")]
         lemna_p19 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P19"),which(substr(colnames(lemna_hamdist),1,3) == "P19")]
         lemna_p27 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P27"),which(substr(colnames(lemna_hamdist),1,3) == "P27")]
         lemna_p36 = lemna_hamdist[which(substr(colnames(lemna_hamdist),1,3) == "P36"),which(substr(colnames(lemna_hamdist),1,3) == "P36")]
         
         ## calculate permutations
         lemna_perm_mean = vector()
         lemna_perm_clone = vector()
         lemna_perm_diversity = vector()
         for (n in 1:iter) {
           
           sim_size = sample(20:27,1)
           runner_sample = sample(colnames(lemna_hamdist), sim_size)
           runner_matrix = lemna_hamdist[runner_sample,runner_sample]
           lemna_perm_mean[n] = mean(runner_matrix[lower.tri(runner_matrix, diag=FALSE)])
           lemna_perm_clone[n] = ncol(hamdist_to_rarecurve(runner_matrix))
           lemna_perm_diversity[n] = diversity(hamdist_to_rarecurve(runner_matrix), index="simpson")
         }
         
         ## mean plot
         par(mar=c(0.3,0.3,0.3,0.3))
         screen(3)
         ## base plot
         plot(NULL, xlim=c(0.5,7.5), ylim=c(0,0.22), xlab="", main="", xaxt = "n", yaxt="n", ylab="Genetic distance", las=2)
         ## mean + confidence interval of population
         boot_means <- replicate(boot_iter, mean(sample(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 7, y0 = ci[1], x1 = 7, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(7, mean(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)]), pch=21, bg=lemna_col, cex=1.5)
         abline(v=6.5)
         ## means + confidence intervals of the deep sampled ponds
         boot_means <- replicate(boot_iter, mean(sample(lemna_p10[lower.tri(lemna_p10, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 1, y0 = ci[1], x1 = 1, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(1, mean(lemna_p10[lower.tri(lemna_p10, diag=FALSE)]), pch=21, bg=lemna_col, cex=1.5)
         boot_means <- replicate(boot_iter, mean(sample(lemna_p14[lower.tri(lemna_p14, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 2, y0 = ci[1], x1 = 2, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(2, mean(lemna_p14[lower.tri(lemna_p14, diag=FALSE)]), pch=21, bg=lemna_col, cex=1.5)
         boot_means <- replicate(boot_iter, mean(sample(lemna_p19[lower.tri(lemna_p19, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 3, y0 = ci[1], x1 = 3, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(3, mean(lemna_p19[lower.tri(lemna_p19, diag=FALSE)]), pch=21, bg=lemna_col, cex=1.5)
         boot_means <- replicate(boot_iter, mean(sample(lemna_p27[lower.tri(lemna_p27, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 4, y0 = ci[1], x1 = 4, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(4, mean(lemna_p27[lower.tri(lemna_p27, diag=FALSE)]), pch=21, bg=lemna_col, cex=1.5)
         boot_means <- replicate(boot_iter, mean(sample(lemna_p36[lower.tri(lemna_p36, diag=FALSE)], replace = TRUE)))
         ci = quantile(boot_means, c(0.025, 0.975))
         arrows(x0 = 5, y0 = ci[1], x1 = 5, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(5, mean(lemna_p36[lower.tri(lemna_p36, diag=FALSE)]), pch=21, bg=lemna_col, cex=1.5)
         ## permutation results
         clip(0,6,0,0.3)
         abline(h=quantile(lemna_perm_mean, probs = c(0.025, 0.975)), lty=2, lwd=1)
         clip(0,6.5,0,0.3)
         points(rep(6,iter), lemna_perm_mean, pch=21, bg=lemna_col, cex=1.5)
         close.screen(3)
         
         ## diversity plot
         par(mar=c(3,0.3,0.3,0.3))
         screen(4)
         ## base plot
         plot(NULL, xlim=c(0.5,7.5), ylim=c(0,1.05), xlab="", main="", xaxt = "n", yaxt = "n", ylab="Diversity (Shannon)", las=2)
         axis(1, at = 1:7, las=2, labels = c("P10", "P14", "P19", "P27", "P36", "perm", "pop"))
         ## plot diversity and CI for population
         data_individuals = rep(seq_along(hamdist_to_rarecurve(lemna_hamdist)), hamdist_to_rarecurve(lemna_hamdist))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 7, y0 = ci[1], x1 = 7, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(7, diversity(hamdist_to_rarecurve(lemna_hamdist), index = "simpson"), pch=21, bg=lemna_col, cex=1.5)
         abline(v=6.5)
         ## plot diversity and CI for deep sampled ponds
         data_individuals = rep(seq_along(hamdist_to_rarecurve(lemna_p10)), hamdist_to_rarecurve(lemna_p10))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 1, y0 = ci[1], x1 = 1, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(1, diversity(hamdist_to_rarecurve(lemna_p10), index = "simpson"), pch=21, bg=lemna_col, cex=1.5)
         data_individuals = rep(seq_along(hamdist_to_rarecurve(lemna_p14)), hamdist_to_rarecurve(lemna_p14))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 2, y0 = ci[1], x1 = 2, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(2, diversity(hamdist_to_rarecurve(lemna_p14), index = "simpson"), pch=21, bg=lemna_col, cex=1.5)
         data_individuals = rep(seq_along(hamdist_to_rarecurve(lemna_p19)), hamdist_to_rarecurve(lemna_p19))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 3, y0 = ci[1], x1 = 3, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(3, diversity(hamdist_to_rarecurve(lemna_p19), index = "simpson"), pch=21, bg=lemna_col, cex=1.5)
         data_individuals = rep(seq_along(hamdist_to_rarecurve(lemna_p27)), hamdist_to_rarecurve(lemna_p27))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 4, y0 = ci[1], x1 = 4, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(4, diversity(hamdist_to_rarecurve(lemna_p27), index = "simpson"), pch=21, bg=lemna_col, cex=1.5)
         data_individuals = rep(seq_along(hamdist_to_rarecurve(lemna_p36)), hamdist_to_rarecurve(lemna_p36))
         boot_out = boot(data = data_individuals, statistic = diversity_stat, R = boot_iter)
         ci = quantile(boot_out$t, c(0.025, 0.975))
         arrows(x0 = 5, y0 = ci[1], x1 = 5, y1 = ci[2], angle = 90, code = 3, length = 0.05)
         points(5, diversity(hamdist_to_rarecurve(lemna_p36), index = "simpson"), pch=21, bg=lemna_col, cex=1.5)
         ## permutations
         clip(0,6,0,4)
         abline(h=quantile(lemna_perm_diversity, probs = c(0.001, 0.999)), lty=2, lwd=1)
         clip(0,6.5,0,4)
         points(rep(6,iter), lemna_perm_diversity, pch=21, bg=lemna_col, cex=1.5)
         close.screen(4)
         
         
     ## POPULATION: Kinship Heatmaps ####
     
     ## LANDOLTIA
     
       ## set parameters
       col_pal = colorRampPalette(c(landoltia_col, "white"))
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
       col_pal = colorRampPalette(c(lemna_col, "white"))
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
       
     ## POPULATION: Distance Decay Plots ####
     
     ## read and transform coordinates
     all_coordinates = read.csv("C:/Users/timte/Desktop/Brisbane/Chapter 1/Second run early 2025/duckweed_coordinates.csv")
     all_coordinates$latitude = sapply(all_coordinates[,"GPS_S"], convert_dmm_to_dd)
     all_coordinates$longitude = sapply(all_coordinates[,"GPS_E"], convert_dmm_to_dd)
     
     ## assemble plotting window
     split.screen(rbind(c(0, 0.5, 0, 1), 
                        c(0.5, 1, 0, 1)))
     
     ## LANDOLTIA
     
       ## subset to landoltia
       landoltia_coor = all_coordinates[which(colnames(landoltia_hamdist) %in% all_coordinates[,"ID"]),]
       
       ## calculate geographic distance
       landoltia_geodist = round(distm(landoltia_coor[, c("longitude", "latitude")], fun = distVincentyEllipsoid),2)
       rownames(landoltia_geodist) = landoltia_coor[,"ID"]; colnames(landoltia_geodist) = landoltia_coor[,"ID"]
       
       ## reduce matrix to triangle
       landoltia_hamdist_vector = landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)]
       landoltia_geodist_vector = landoltia_geodist[lower.tri(landoltia_geodist, diag=FALSE)]
       
       ## convert to km
       landoltia_geodist_vector = landoltia_geodist_vector/1000
       
       ## assemble plot
       screen(1)
       par(mar=c(4,4,2,0.2))
       plot(landoltia_geodist_vector, landoltia_hamdist_vector,
            pch=21, bg=landoltia_col, cex=0.6, ylim=c(-0.01,max(c(lemna_hamdist,landoltia_hamdist))),
            xlab="Geographic distance km", 
            ylab="Genetic distance")
       lines(x = c(min(landoltia_geodist_vector), max(landoltia_geodist_vector)),
             y = c(mean(landoltia_hamdist_vector), mean(landoltia_hamdist_vector)), lwd=2)
       clip(min(landoltia_geodist_vector), max(landoltia_geodist_vector),0,0.8)
       abline(h=0.02, lty=3, lwd=2)
       abline(lm(landoltia_hamdist_vector ~ landoltia_geodist_vector), col="red", lwd=2)
       
       ## legend
       usr = par("usr"); clip(usr[1], usr[2], usr[3], usr[4])
       legend("topleft", inset=0.02, legend=c("mean", "model", "clone cutoff"),
              lty=c(1,1,3), col=c("black", "red", "black"), lwd=c(2,2,2))
       close.screen(1)
       
     ## LEMNA
       
       ## subset to lemna
       lemna_coor = all_coordinates[which(colnames(lemna_hamdist) %in% all_coordinates[,"ID"]),]
       
       ## calculate geographic distance
       lemna_geodist = round(distm(lemna_coor[, c("longitude", "latitude")], fun = distVincentyEllipsoid),2)
       rownames(lemna_geodist) = lemna_coor[,"ID"]; colnames(lemna_geodist) = lemna_coor[,"ID"]
       
       ## reduce matrix to triangle
       lemna_hamdist_vector = lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)]
       lemna_geodist_vector = lemna_geodist[lower.tri(lemna_geodist, diag=FALSE)]
       
       ## transform to km
       lemna_geodist_vector = lemna_geodist_vector/1000
       
       ## assemble plot
       screen(2)
       par(mar=c(4,0.2,2,4))
       plot(lemna_geodist_vector, lemna_hamdist_vector,
            pch=21, bg=lemna_col, cex=0.6, ylim=c(-0.01,max(c(lemna_hamdist,lemna_hamdist))),
            xlab="Geographic distance km", ylab="", yaxt = "n")
       lines(x = c(min(lemna_geodist_vector), max(lemna_geodist_vector)),
             y = c(mean(lemna_hamdist_vector), mean(lemna_hamdist_vector)), lwd=2)
       clip(min(lemna_geodist_vector), max(lemna_geodist_vector),0,0.8)
       abline(h=0.02, lty=3, lwd=2)
       abline(lm(lemna_hamdist_vector ~ lemna_geodist_vector), col="red", lwd=2)
       close.screen(2)
       
     ## Mantel tests
     
     ## LANDOLTIA
     vegan::mantel(landoltia_hamdist, landoltia_geodist, permutations = 1000)
     
     ## LEMNA
     vegan::mantel(lemna_hamdist, lemna_geodist, permutations = 1000)
     
     ## POPULATION: PCAs ####
     
     ## transform files
     landoltia_final_pca = vcfR2genlight(landoltia_final)
     landoltia_clone_no_missing_data = tab(landoltia_final_pca, NA.method = "mean")
     
     landoltia_final_no_clone_pca = vcfR2genlight(landoltia_final_no_clone)
     landoltia_no_clone_no_missing_data = tab(landoltia_final_no_clone_pca, NA.method = "mean")
     
     lemna_final_pca = vcfR2genlight(lemna_final)
     lemna_clone_no_missing_data = tab(lemna_final_pca, NA.method = "mean")
     
     lemna_final_no_clone_pca = vcfR2genlight(lemna_final_no_clone)
     lemna_no_clone_no_missing_data = tab(lemna_final_no_clone_pca, NA.method = "mean")
     
     
     par(mfrow=c(1,2))
     
     ## LANDOLTIA
       
         ## run PCA
         landoltia_dudipca_result = dudi.pca(landoltia_clone_no_missing_data, scannf = FALSE, nf = 3)
         landoltia_pca_plotter = landoltia_dudipca_result$li
         
         ## add clone identity
         landoltia_cloner = rep(NA, nrow(landoltia_pca_plotter))
         for (n in rownames(landoltia_pca_plotter)) {landoltia_cloner[which(rownames(landoltia_pca_plotter)==n)] = colnames(landoltia_clone_df)[which(landoltia_clone_df == n, arr.ind = TRUE)[2]]}
         landoltia_cloner = paste0("c",sub(".*?(\\d+)$", "\\1", landoltia_cloner))
         landoltia_pca_plotter$cloneID = landoltia_cloner
         
         ## assemble plpt
         plot(landoltia_pca_plotter[,"Axis1"], landoltia_pca_plotter[,"Axis2"], 
              xlab = paste0("PC1 (", round((dudipca_result$eig[1] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
              ylab = paste0("PC2 (", round((dudipca_result$eig[2] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
              main = "landoltia", pch = 21, 
              bg = ifelse(grepl("P10", rownames(dudipca_result$li)), P10_col,
                          ifelse(grepl("P14", rownames(dudipca_result$li)), P14_col, 
                                 ifelse(grepl("P19", rownames(dudipca_result$li)), P19_col,
                                        ifelse(grepl("P27", rownames(dudipca_result$li)), P27_col,
                                               ifelse(grepl("P36", rownames(dudipca_result$li)), P36_col, rest_col))))))
         text(landoltia_pca_plotter[,"Axis1"], landoltia_pca_plotter[,"Axis2"], 
              labels=landoltia_pca_plotter[,"cloneID"], cex=0.7, pos=sample(1:4,nrow(landoltia_pca_plotter), replace = TRUE),
              col="gray30")
         legend("bottomright", legend=c("P10","P14", "P19", "P27", "P36", "other"),
                pt.bg=c(P10_col,P14_col, P19_col, P27_col, P36_col, rest_col),
                pch=21, bty="n")
             
     ## LEMNA
     
         ## run PCA
         lemna_dudipca_result = dudi.pca(lemna_clone_no_missing_data, scannf = FALSE, nf = 3)
         lemna_pca_plotter = lemna_dudipca_result$li
         
         ## add clone identity
         lemna_cloner = rep(NA, nrow(lemna_pca_plotter))
         for (n in rownames(lemna_pca_plotter)) {lemna_cloner[which(rownames(lemna_pca_plotter)==n)] = colnames(lemna_clone_df)[which(lemna_clone_df == n, arr.ind = TRUE)[2]]}
         lemna_cloner = paste0("c",sub(".*?(\\d+)$", "\\1", lemna_cloner))
         lemna_pca_plotter$cloneID = lemna_cloner
         
         ## assemble plpt
         plot(lemna_pca_plotter[,"Axis1"], lemna_pca_plotter[,"Axis2"], 
              xlab = paste0("PC1 (", round((dudipca_result$eig[1] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
              ylab = paste0("PC2 (", round((dudipca_result$eig[2] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
              main = "lemna", pch = 21, 
              bg = ifelse(grepl("P10", rownames(dudipca_result$li)), P10_col,
                          ifelse(grepl("P14", rownames(dudipca_result$li)), P14_col, 
                                 ifelse(grepl("P19", rownames(dudipca_result$li)), P19_col,
                                        ifelse(grepl("P27", rownames(dudipca_result$li)), P27_col,
                                               ifelse(grepl("P36", rownames(dudipca_result$li)), P36_col, rest_col))))))
         text(lemna_pca_plotter[,"Axis1"], lemna_pca_plotter[,"Axis2"], 
              labels=lemna_pca_plotter[,"cloneID"], cex=0.7, pos=sample(1:4,nrow(lemna_pca_plotter), replace = TRUE),
              col="gray30")
         legend("topright", legend=c("P10","P14", "P19", "P27", "P36", "other"),
                pt.bg=c(P10_col,P14_col, P19_col, P27_col, P36_col, rest_col),
                pch=21, bty="n")
          
         
     ## mantel test and the like ####
     
     library(ecodist)
     MRM(as.dist(landoltia_hamdist) ~ as.dist(landoltia_geodist), nperm = 999)
     
     mantel(landoltia_hamdist,landoltia_geodist)
     
     
     
     
     
     
     ## STRUCTURE ANALYSIS ####
     
     ## set dir
     lemna_str_dir = "C:/Users/timte/Desktop/Brisbane/Chapter 1/STRUCTURE/Lemna/lemna_structure_results"
     
     ## get filenames
     file_names = list.files(path= lemna_str_dir, pattern="str_K.*_rep.*_f")
     
     ## extract K, rep and LnProb
     str_analysis_df = as.data.frame(matrix(ncol=3, nrow=0))
     for (n in file_names) {
     
       ## load file
       runner_text = readLines(paste0(lemna_str_dir,"/", n, sep=""), n=300)
       
       ## Get estimated Ln prob of Data
       lnprob_full = grep("Estimated Ln Prob of Data", runner_text, value=TRUE)
       lnprob = as.numeric(regmatches(lnprob_full, regexpr("[-]?[0-9]+\\.?[0-9]*", lnprob_full)))
       
       ## Get K and replicate
       k = as.numeric(sub("str_K(\\d+)_rep.*", "\\1", n))
       rep = as.numeric(sub(".*_rep(\\d+)_f", "\\1", n))
       
       ## save in df
       str_analysis_df = rbind(str_analysis_df, c(k, rep, lnprob))
     }
     colnames(str_analysis_df) = c("K", "rep", "LnProb")
     ## plot lnProb against K
     plot(str_analysis_df[,"K"], str_analysis_df[,"LnProb"],
          xlab="K", ylab="Ln probability", main="Lemna")
     
     ## stacked barplots
     
     ## extract table for barplots
     str_list = list()
     for (n in file_names) {
       
       ## load file
       runner_text = readLines(paste0(lemna_str_dir,"/", n, sep=""), n=300)
        
       ## extract table
       start = grep("Inferred ancestry of individuals:", runner_text)
       ancestry_lines = runner_text[(start+2):length(runner_text)]
       end = which(ancestry_lines == "")[1] - 1
       ancestry_lines = ancestry_lines[1:end]
       
       # Use regex to extract Label, Pop, and the cluster values
       matches <- str_match(ancestry_lines, "\\s*\\d+\\s+(\\S+)\\s+\\(\\d+\\)\\s+(\\d+)\\s+:\\s+(.+)")
       
       ## make dataframe
       df = data.frame(Label = matches[,2], Pop   = as.integer(matches[,3]), Clusters = matches[,4], stringsAsFactors = FALSE)
       cluster_matrix = str_split(df$Clusters, "\\s+", simplify = TRUE)
       cluster_matrix = apply(cluster_matrix, 2, as.numeric)
       df_final <- cbind(df[,1:2], cluster_matrix[,1:(ncol(cluster_matrix)-1)])
       colnames(df_final) <- c("Label", "Pop", paste0("Cluster", 1:(ncol(cluster_matrix)-1)))
       
       str_list[[which(file_names==n)]] = df_final
     }
     names(str_list) = file_names
     for (n in 1:length(str_list)) {
       
       runner_df = str_list[[n]]
       mat = as.matrix(runner_df[,2:ncol(runner_df)])
       png(filename = paste0("Lemna ",names(str_list[n]),".png"), res=100, width = 600, height=600)
       
       bp = barplot(t(mat[,2:ncol(mat)]), col=rainbow(ncol(mat)), border=NA, space=0,
                    xlab="Individuals", ylab="Ancestry proportion", main=names(str_list[n]),
                    xaxt = "n")
       axis(1, at = bp, labels = df_final[,"Label"], las = 2, cex.axis = 0.7)
       
       dev.off()
       
     }
     
     ## set dir
     landoltia_str_dir = "C:/Users/timte/Desktop/Brisbane/Chapter 1/STRUCTURE/Landoltia/landoltia_structure_results"
                      
     ## get filenames
     file_names = list.files(path= landoltia_str_dir, pattern="str_K.*_rep.*_f")
     
     ## extract K, rep and LnProb
     str_analysis_df = as.data.frame(matrix(ncol=3, nrow=0))
     for (n in file_names) {
       
       ## load file
       runner_text = readLines(paste0(landoltia_str_dir,"/", n, sep=""), n=300)
       
       ## Get estimated Ln prob of Data
       lnprob_full = grep("Estimated Ln Prob of Data", runner_text, value=TRUE)
       lnprob = as.numeric(regmatches(lnprob_full, regexpr("[-]?[0-9]+\\.?[0-9]*", lnprob_full)))
       
       ## Get K and replicate
       k = as.numeric(sub("str_K(\\d+)_rep.*", "\\1", n))
       rep = as.numeric(sub(".*_rep(\\d+)_f", "\\1", n))
       
       ## save in df
       str_analysis_df = rbind(str_analysis_df, c(k, rep, lnprob))
     }
     colnames(str_analysis_df) = c("K", "rep", "LnProb")
     ## plot lnProb against K
     plot(str_analysis_df[,"K"], str_analysis_df[,"LnProb"],
          xlab="K", ylab="Ln probability", main="Landoltia")
     
     ## stacked barplots
     
     ## extract table for barplots
     str_list = list()
     for (n in file_names) {
       
       ## load file
       runner_text = readLines(paste0(landoltia_str_dir,"/", n, sep=""), n=300)
       
       ## extract table
       start = grep("Inferred ancestry of individuals:", runner_text)
       ancestry_lines = runner_text[(start+2):length(runner_text)]
       end = which(ancestry_lines == "")[1] - 1
       ancestry_lines = ancestry_lines[1:end]
       
       # Use regex to extract Label, Pop, and the cluster values
       matches <- str_match(ancestry_lines, "\\s*\\d+\\s+(\\S+)\\s+\\(\\d+\\)\\s+(\\d+)\\s+:\\s+(.+)")
       
       ## make dataframe
       df = data.frame(Label = matches[,2], Pop   = as.integer(matches[,3]), Clusters = matches[,4], stringsAsFactors = FALSE)
       cluster_matrix = str_split(df$Clusters, "\\s+", simplify = TRUE)
       cluster_matrix = apply(cluster_matrix, 2, as.numeric)
       df_final <- cbind(df[,1:2], cluster_matrix[,1:(ncol(cluster_matrix)-1)])
       colnames(df_final) <- c("Label", "Pop", paste0("Cluster", 1:(ncol(cluster_matrix)-1)))
       
       str_list[[which(file_names==n)]] = df_final
     }
     names(str_list) = file_names
     for (n in 1:length(str_list)) {
       
       runner_df = str_list[[n]]
       mat = as.matrix(runner_df[,2:ncol(runner_df)])
       png(filename = paste0("Landoltia ",names(str_list[n]),".png"), res=100, width = 1000, height=1000)
       
       bp = barplot(t(mat[,2:ncol(mat)]), col=rainbow(ncol(mat)), border=NA, space=0,
                    xlab="Individuals", ylab="Ancestry proportion", main=names(str_list[n]),
                    xaxt = "n")
       axis(1, at = bp, labels = df_final[,"Label"], las = 2, cex.axis = 0.5)
       
       dev.off()
       
     }
     
     
     
     ## SUPPLEMENT PLOTS ####
     
     ## clone cutoff
     split.screen(rbind(c(0, 0.49, 0, 1),
                        c(0.51, 1, 0, 1)))
     
     screen(1)
     par(mar=c(4.5,4.5,2,0))
     hist(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)], breaks=100,
          main = "Landoltia", col=landoltia_col, xlab="", ylim=c(0,1100)); box()
     abline(v=0.02, col="red", lty=2)
     close.screen(1)
     
     
     screen(2)
     lemna_hist = hist(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)], breaks=100,
                       main = "Lemna", col=lemna_col, xlab="", ylim=c(0,1100), plot=FALSE)
     par(mar=c(4.5,0,2,4.5))
     plot(NULL, xlim=c(min(lemna_hist$breaks), max(lemna_hist$breaks)), ylim=c(0,1100), xlab="", yaxt="n", main="Lemna")
     plot(lemna_hist, add=TRUE, col=lemna_col)
     abline(v=0.02, col="red", lty=2)
     close.screen(2)
     close.screen(all.screens = TRUE)
     
     par(xpd = NA)  # allow drawing outside plot region
     text(x = 0.33, y = -0.2,labels = "Pairwise genetic distances")
     par(xpd = FALSE)
     
     ## rarefaction curves
     
     ## subset data, and compute rarecurves
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
     
     landoltia_boxplotter = sapply(apply(landoltia_clone_df, 2, function(x) unique(na.omit(substr(x,1,3)))),length)
     lemna_boxplotter = sapply(apply(lemna_clone_df, 2, function(x) unique(na.omit(substr(x,1,3)))),length)
     
     ## assemble plot
     par(mfrow=c(2,2))
     
     ## plot 1 - population rarefaction
     plot(1:length(landoltia_rareplot_all[[1]]),landoltia_rareplot_all[[1]], type="l", lwd=3, col=landoltia_col,
          xlab="sampled individuals", ylab="clones found", main="Meta-population")
     legend("topleft", inset=0.02, legend=c("Landoltia", "Lemna"), 
            col=c(landoltia_col, lemna_col), lty=1, lwd=3)
     lines(1:length(lemna_rareplot_all[[1]]),lemna_rareplot_all[[1]], type="l", lwd=3, col=lemna_col)
     
     ## plot 2 - clone distributions
     plot(1:length(landoltia_boxplotter), rev(sort(landoltia_boxplotter)), main="Clone distribution",
          xlab="clones", ylab="# of ponds each clone occurs in", col=landoltia_col, lwd=3, type="l",
          xlim=c(0,68), ylim=c(1,6))
     lines(1:length(lemna_boxplotter), rev(sort(lemna_boxplotter)),
           col=lemna_col, lwd=3)
     legend("topright", inset=0.02, legend=c("Landoltia", "Lemna"), 
            col=c(landoltia_col, lemna_col), lty=1, lwd=3)
     
     ## this needs work - lines are wrong, needs to be dots or something
     
     ## plot 3 - landoltia rarefaction
     plot(1:length(landoltia_rareplot_P10[[1]]),landoltia_rareplot_P10[[1]], col=P10_col, lwd=3, type="l",
          ylim=c(1,13), xlim=c(0,27), xlab="sampled individuals", ylab="clones found", main="Landoltia")
     legend("topleft", inset=0.02, legend=c("P10", "P14", "P19", "P27", "P36"), 
            col=c(P10_col, P14_col, P19_col, P27_col, P36_col), lty=1, lwd=3)
     lines(1:length(landoltia_rareplot_P10[[1]]),landoltia_rareplot_P10[[1]], col=P10_col, lwd=3)
     lines(1:length(landoltia_rareplot_P14[[1]]),landoltia_rareplot_P14[[1]], col=P14_col, lwd=3)
     lines(1:length(landoltia_rareplot_P19[[1]]),landoltia_rareplot_P19[[1]], col=P19_col, lwd=3)
     lines(1:length(landoltia_rareplot_P27[[1]]),landoltia_rareplot_P27[[1]], col=P27_col, lwd=3)
     lines(1:length(landoltia_rareplot_P36[[1]]),landoltia_rareplot_P36[[1]], col=P36_col, lwd=3)
     
     ## plot 4 - lemna rarefaction
     plot(1:length(lemna_rareplot_P10[[1]]),lemna_rareplot_P10[[1]], col=P10_col, lwd=3, type="l",
          ylim=c(1,4), xlim=c(0,27), xlab="sampled individuals", ylab="clones found", main="Lemna",
          yaxt="n")
     axis(2, at=c(1,2,3,4), labels=c(1,2,3,4))
     legend("topleft", inset=0.02, legend=c("P10", "P14", "P19", "P27", "P36"), 
            col=c(P10_col, P14_col, P19_col, P27_col, P36_col), lty=1, lwd=3)
     lines(1:length(lemna_rareplot_P10[[1]]),lemna_rareplot_P10[[1]], col=P10_col, lwd=3)
     lines(1:length(lemna_rareplot_P14[[1]]),lemna_rareplot_P14[[1]], col=P14_col, lwd=3)
     lines(1:length(lemna_rareplot_P19[[1]]),lemna_rareplot_P19[[1]], col=P19_col, lwd=3)
     lines(1:length(lemna_rareplot_P27[[1]]),lemna_rareplot_P27[[1]], col=P27_col, lwd=3)
     lines(1:length(lemna_rareplot_P36[[1]]),lemna_rareplot_P36[[1]], col=P36_col, lwd=3)
     
     
     
     
## scraps ####
     ## subsetting alive ####
     
     landoltia_alive = c("P13S1", 
                         "P14S7",
                         "P14S16", 
                         "P14S42", 
                         "P14S22", 
                         "P19S36",
                         "P19S16", 
                         "P19S52", 
                         "P27S13",
                         "P27S3",
                         "P28S3", 
                         "P36S1", 
                         "P36S41", 
                         "P36S34",
                         "P11S12", 
                         "P12S4", 
                         "P12S5", 
                         "P12S6", 
                         "P14S18",
                         "P19S11", 
                         "P19S23", 
                         "P19S34", 
                         "P19S54", 
                         "P23S4",
                         "P23S5", 
                         "P25S1", 
                         "P26S1", 
                         "P26S3")
     
     lemna_alive = c("P7S2",
                     "P10S28",
                     "P33S3",
                     "P11S8",
                     "P19S30",
                     "P27S4",
                     "P32S5",
                     "P14S43")
       
     landoltia_alive_hamdist = landoltia_hamdist[which(colnames(landoltia_hamdist) %in% landoltia_alive),which(colnames(landoltia_hamdist) %in% landoltia_alive)]
     lemna_alive_hamdist = lemna_hamdist[which(colnames(lemna_hamdist) %in% lemna_alive),which(colnames(lemna_hamdist) %in% lemna_alive)]
     
     write.csv(lemna_alive_hamdist, file="lemna_alive.csv")
     write.csv(landoltia_alive_hamdist, file="landoltia_alive.csv")
     
     ## plot landoltia
     col_pal = colorRampPalette(c(landoltia_col, "white"))
     colbreaks = 10
     clone_cutoff = 0.02
     image(1:ncol(landoltia_alive_hamdist), 1:nrow(landoltia_alive_hamdist),
           landoltia_alive_hamdist, axes = FALSE, frame=FALSE, xlab="", ylab="",
           col=c("black",col_pal(colbreaks-2)), main="Landoltia",
           breaks=c(0,seq(from=clone_cutoff, to=max(c(landoltia_hamdist,lemna_hamdist)), length.out = colbreaks-1)))
     mtext(colnames(landoltia_alive_hamdist), side=1, at=1:length(landoltia_alive), las=2)
     mtext(colnames(landoltia_alive_hamdist), side=2, at=1:length(landoltia_alive), las=2)
     
     ## lemna
     col_pal = colorRampPalette(c(lemna_col, "white"))
     image(1:ncol(lemna_alive_hamdist), 1:nrow(lemna_alive_hamdist),
           lemna_alive_hamdist, axes = FALSE, frame=FALSE, xlab="", ylab="",
           col=c("black",col_pal(colbreaks-2)), main="Lemna",
           breaks=c(0,seq(from=clone_cutoff, to=max(c(lemna_hamdist,lemna_hamdist)), length.out = colbreaks-1)))
     mtext(colnames(lemna_alive_hamdist), side=1, at=1:length(lemna_alive), las=2)
     mtext(colnames(lemna_alive_hamdist), side=2, at=1:length(lemna_alive), las=2)
     
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
     lemna_coordinates$longitude = sapply(lemna_coordinates[,"GPS_E"], convert_dmm_to_dd)
     rownames(lemna_coordinates) = lemna_coordinates[,"ID"]
     
     ## calculate geographic distance
     lemna_geodist_matrix = distm(lemna_coordinates[, c("longitude", "latitude")], fun = distHaversine)
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
     all_coordinates$longitude = sapply(all_coordinates[,"GPS_E"], convert_dmm_to_dd)
     
     ## assemble plot
     
     par(mfrow=c(1,2))
     
     ## full dataset
     
     ## subset to lemna
     lemna_coor = all_coordinates[which(colnames(lemna_hamdist) %in% all_coordinates[,"ID"]),]
     
     ## calculate geographic distance
     lemna_geodist = round(distm(lemna_coor[, c("longitude", "latitude")], fun = distVincentyEllipsoid),2)
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
     lemna_noclone_geodist = round(distm(lemna_noclone_coor[, c("longitude", "latitude")], fun = distVincentyEllipsoid),2)
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
     landoltia_noclone_geodist = round(distm(landoltia_noclone_coor[, c("longitude", "latitude")], fun = distVincentyEllipsoid),2)
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
     perm_clone = vector()
     perm_diversity = vector()
     for (n in 1:100000) {
       
       runner_sample = sample(colnames(landoltia_hamdist), 22)
       runner_matrix = landoltia_hamdist[runner_sample,runner_sample]
       perm_mean[n] = mean(runner_matrix[lower.tri(runner_matrix, diag=FALSE)])
       perm_clone[n] = ncol(hamdist_to_rarecurve(runner_matrix))
       perm_diversity[n] = diversity(hamdist_to_rarecurve(runner_matrix))
     }
     points(rep(6,100000), perm_mean)
     abline(h=quantile(perm_mean, probs = c(0.001, 0.999)))
     
     hist(perm_clone)
     hist(perm_diversity)
     
     
     
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
     
     
     
     ## FIGURE X: triplet sampling ####
     
     par(mfrow=c(1,2))
     
     ## LANDOLTIA triplets    
     landoltia_same_location = data.frame(site1 = c("P1S2",NA,NA),
                                          site2 = c("P2S2",NA,NA),
                                          site3 = c("P4S2",NA,NA),
                                          site4 = c("P10S1b", "P10S3b", NA),
                                          site5 = c("P10S3",NA,NA),
                                          site6 = c("P10S4","P10S6", NA),
                                          site7 = c("P10S4b",NA, NA),
                                          site8 = c("P10S7b", "P10S8b", "P10S9b"),
                                          site9 = c("P10S10b", "P10S11b", "P10S12b"),
                                          site10 = c("P10S13b","P10S14b", NA),
                                          site11 = c("P10S16", "P10S18", NA),
                                          site12 = c("P10S22", "P10S23", NA),
                                          site13 = c("P10S26", NA, NA),
                                          site14 = c("P10S31", NA, NA),
                                          site15 = c("P11S10", "P11S11", "P11S12"),
                                          site16 = c("P12S4", "P12S5", "P12S6"),
                                          site17 = c("P13S1", "P13S2", "P13S3"),
                                          site18 = c("P14S1", NA, NA),
                                          site19 = c("P14S7", "P14S8", "P14S9"),
                                          site20 = c("P14S12", "P14S13", NA),
                                          site21 = c("P14S16", "P14S17", "P14S18"),
                                          site22 = c("P14S22", "P14S23", "P14S24"),
                                          site23 = c("P14S28", "P14S29", "P14S30"),
                                          site24 = c("P14S37", "P14S38", "P14S39"),
                                          site25 = c("P14S40", "P14S42", NA),
                                          site26 = c("P14S46", "P14S47", "P14S48"),
                                          site27 = c("P15S1", "P15S2", "P15S3"),
                                          site28 = c("P16S4", "P16S5", "P16S6"),
                                          site29 = c("P17S4", "P17S6", NA),
                                          site30 = c("P18S4", "P18S6", NA),
                                          site31 = c("P19S4", "P19S5", NA),
                                          site32 = c("P19S10", "P19S11", NA),
                                          site33 = c("P19S16", "P19S17", "P19S18"),
                                          site34 = c("P19S22", "P19S23", "P19S24"),
                                          site35 = c("P19S25", "P19S26", "P19S27"),
                                          site36 = c("P19S34", "P19S35", "P19S36"),
                                          site37 = c("P19S40", "P19S41", "P19S42"),
                                          site38 = c("P19S43", "P19S44", "P19S45"),
                                          site39 = c("P19S52", "P19S53", "P19S54"),
                                          site40 = c("P20S1", "P20S2", "P20S3"),
                                          site41 = c("P22S1", "P22S2", "P22S3"),
                                          site42 = c("P23S4", NA, NA),
                                          site43 = c("P23S5", NA, NA),
                                          site44 = c("P24S1", "P24S2", "P24S3"),
                                          site45 = c("P25S1", "P25S2", "P25S3"),
                                          site46 = c("P26S1", "P26S2", "P26S3"),
                                          site47 = c("P27S1", "P27S2", "P27S3"),
                                          site48 = c("P27S7", "P27S8", "P27S9"),
                                          site49 = c("P27S13", "P27S14", "P27S15"),
                                          site50 = c("P27S22", "P27S23", "P27S24"),
                                          site51 = c("P27S26", "P27S28", NA),
                                          site52 = c("P27S30", "P27S31", "P27S32"),
                                          site53 = c("P27S33", "P27S34", "P27S35"),
                                          site54 = c("P27S36", "P27S37", "P27S38"),
                                          site55 = c("P27S49", "P27S51", NA),
                                          site56 = c("P28S1", "P28S2", "P28S3"),
                                          site57 = c("P32S1", "P32S2", "P32S3"),
                                          site58 = c("P34S1", "P34S2", "P34S3"),
                                          site59 = c("P36S1", "P36S2", "P36S3"),
                                          site60 = c("P36S7", "P36S8", "P36S9"),
                                          site61 = c("P36S16", "P36S17", "P36S18"),
                                          site62 = c("P36S19", "P36S20", "P36S21"),
                                          site63 = c("P36S25", "P36S26", "P36S27"),
                                          site64 = c("P36S34", "P36S35", "P36S36"),
                                          site65 = c("P36S40", "P36S41", NA))
     
     
     
     ## LEMNA triplets
     lemna_same_location = data.frame(site2 = c("P5S1", NA, NA),
                                      site3 = c("P6S1", "P6S3",NA),
                                      site4 = c("P7S1", "P7S2", "P7S3"),
                                      site5 = c("P10S1", "P10S2", NA),
                                      site6 = c("P10S2b", NA, NA),
                                      site7 = c("P10S5", NA, NA),
                                      site8 = c("P10S5b", "P10S6b", NA),
                                      site9 = c("P10S7", "P10S8", "P10S9"),
                                      site10 = c("P10S10", "P10S11", "P10S12"),
                                      site11 = c("P10S13", "P10S14", NA),
                                      site12 = c("P10S15b", NA, NA),
                                      site13 = c("P10S19", "P10S20", "P10S21"),
                                      site14 = c("P10S28", "P10S30",NA),
                                      site15 = c("P11S7", "P11S8", "P11S9"),
                                      site16 = c("P14S10", "P14S11", NA),
                                      site17 = c("P14S15", NA, NA),
                                      site18 = c("P14S19", "P14S20", "P14S21"),
                                      site19 = c("P14S25", "P14S26", "P14S27"),
                                      site20 = c("P14S31", "P14S32", NA),
                                      site21 = c("P14S34", "P14S35", "P14S36"),
                                      site22 = c("P14S43", "P14S44", "P14S45"),
                                      site23 = c("P14S49", "P14S50", NA),
                                      site24 = c("P14S55", "P14S56", "P14S57"),
                                      site25 = c("P16S2", "P16S3", NA),
                                      site26 = c("P17S1", "P17S2", "P17S3"),
                                      site27 = c("P18S2", NA, NA),
                                      site28 = c("P19S1", "P19S2", "P19S3"),
                                      site29 = c("P19S7", "P19S8", "P19S9"),
                                      site30 = c("P19S13", "P19A14", "P19S15"),
                                      site31 = c("P19S19", "P19S20", "P19S21"),
                                      site32 = c("P19S28", "P19S29", "P19S30"),
                                      site33 = c("P19S31", "P19S32", "P19S33"),
                                      site34 = c("P19S37", "P19S38", "P19S39"),
                                      site35 = c("P19S46", "P19S47", "P19S48"),
                                      site36 = c("P19S49", "P19S50", "P19S51"),
                                      site37 = c("P21S1","P21S2", NA),
                                      site38 = c("P22S4", "P22S6",NA),
                                      site39 = c("P23S3", NA, NA),
                                      site40 = c("P27S4", "P27S5", NA),
                                      site41 = c("P27S10", "P27S11", "P27S12"),
                                      site42 = c("P27S16", "P27S17", "P27S18"),
                                      site43 = c("P27S19", "P27S20", "P27S21"),
                                      site44 = c("P27S25", "P27S27", NA),
                                      site45 = c("P27S39", "P27S40", "P27S41"),
                                      site46 = c("P27S42", NA, NA),
                                      site47 = c("P27S43", "P27S44", "P27S45"),
                                      site48 = c("P27S46", "P27S47", "P27S48"),
                                      site49 = c("P27S52", NA, NA),
                                      site50 = c("P28S4", "P28S5", "P28S6"),
                                      site51 = c("P30S1", "P30S2", "P30S3"),
                                      site52 = c("P31S1", "P31S2", "P31S3"),
                                      site53 = c("P32S4", "P32S5", "P32S6"),
                                      site54 = c("P33S1", "P33S2", "P33S3"),
                                      site55 = c("P34S5", "P34S6", NA),
                                      site56 = c("P35S2", "P35S3", NA),
                                      site57 = c("P36S4", "P36S5", "P36S6"),
                                      site58 = c("P36S10", "P36S11", "P36S12"),
                                      site59 = c("P36S13", "P36S14", "P36S15"),
                                      site60 = c("P36S22", "P36S23", "P36S24"),
                                      site61 = c("P36S28", "P36S29", "P36S30"),
                                      site62 = c("P36S31", "P36S32", "P36S33"),
                                      
                                      site63 = c("P36S37", "P36S38", "P36S39")
     )
     
     
     
     ## identify sites with 2 and 3 plants 
     landoltia_sites_with_two_or_more = as.vector(which(apply(apply(landoltia_same_location, 2, function(x) is.na(x)), 2, sum) != 2))
     
     ## remove singles
     landoltia_same_location = landoltia_same_location[,c(landoltia_sites_with_two_or_more)]
     
     landoltia_saver = list()   
     for (n in 1:ncol(landoltia_same_location)) {
       
       runner_matrix = landoltia_hamdist[which(colnames(landoltia_hamdist) %in% landoltia_same_location[,n]),which(colnames(landoltia_hamdist) %in% landoltia_same_location[,n]), drop=FALSE]
       landoltia_saver[[n]] = runner_matrix[lower.tri(runner_matrix, diag=FALSE)]
       
     }    
     
     ## compute histograms with equal bins
     breaks = pretty(range(c(0,max(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)]))), n = 30)
     landoltia_local_samples = hist(unlist(landoltia_saver), breaks = breaks, plot = FALSE)
     landoltia_global_samples = hist(landoltia_hamdist[lower.tri(landoltia_hamdist, diag=FALSE)], breaks = breaks, plot = FALSE)
     
     # assemble triplet plot
     plot(landoltia_local_samples,
          col = scales::alpha("red", 0.5),
          xlab = "Genetic distance",
          ylab = "", main="",
          ylim = c(0, max(landoltia_local_samples$counts)),
          xlim = c(0, max(landoltia_hamdist+0.02)))
     axis(side=2, col.axis="red")
     
     ## add second histogram
     par(new = TRUE)
     plot(landoltia_global_samples,
          col = scales::alpha("blue", 0.5),
          axes = FALSE, xlab = "", ylab = "", main = "",
          ylim = c(0, max(landoltia_global_samples$counts)),
          xlim = c(0, max(landoltia_hamdist+0.02)))
     axis(side = 4, col.axis="blue")
     
     ## Lemna triplets
     lemna_same_location = data.frame(site2 = c("P5S1", NA, NA),
                                      site3 = c("P6S1", "P6S3",NA),
                                      site4 = c("P7S1", "P7S2", "P7S3"),
                                      site5 = c("P10S1", "P10S2", NA),
                                      site6 = c("P10S2b", NA, NA),
                                      site7 = c("P10S5", NA, NA),
                                      site8 = c("P10S5b", "P10S6b", NA),
                                      site9 = c("P10S7", "P10S8", "P10S9"),
                                      site10 = c("P10S10", "P10S11", "P10S12"),
                                      site11 = c("P10S13", "P10S14", NA),
                                      site12 = c("P10S15b", NA, NA),
                                      site13 = c("P10S19", "P10S20", "P10S21"),
                                      site14 = c("P10S28", "P10S30",NA),
                                      site15 = c("P11S7", "P11S8", "P11S9"),
                                      site16 = c("P14S10", "P14S11", NA),
                                      site17 = c("P14S15", NA, NA),
                                      site18 = c("P14S19", "P14S20", "P14S21"),
                                      site19 = c("P14S25", "P14S26", "P14S27"),
                                      site20 = c("P14S31", "P14S32", NA),
                                      site21 = c("P14S34", "P14S35", "P14S36"),
                                      site22 = c("P14S43", "P14S44", "P14S45"),
                                      site23 = c("P14S49", "P14S50", NA),
                                      site24 = c("P14S55", "P14S56", "P14S57"),
                                      site25 = c("P16S2", "P16S3", NA),
                                      site26 = c("P17S1", "P17S2", "P17S3"),
                                      site27 = c("P18S2", NA, NA),
                                      site28 = c("P19S1", "P19S2", "P19S3"),
                                      site29 = c("P19S7", "P19S8", "P19S9"),
                                      site30 = c("P19S13", "P19A14", "P19S15"),
                                      site31 = c("P19S19", "P19S20", "P19S21"),
                                      site32 = c("P19S28", "P19S29", "P19S30"),
                                      site33 = c("P19S31", "P19S32", "P19S33"),
                                      site34 = c("P19S37", "P19S38", "P19S39"),
                                      site35 = c("P19S46", "P19S47", "P19S48"),
                                      site36 = c("P19S49", "P19S50", "P19S51"),
                                      site37 = c("P21S1","P21S2", NA),
                                      site38 = c("P22S4", "P22S6",NA),
                                      site39 = c("P23S3", NA, NA),
                                      site40 = c("P27S4", "P27S5", NA),
                                      site41 = c("P27S10", "P27S11", "P27S12"),
                                      site42 = c("P27S16", "P27S17", "P27S18"),
                                      site43 = c("P27S19", "P27S20", "P27S21"),
                                      site44 = c("P27S25", "P27S27", NA),
                                      site45 = c("P27S39", "P27S40", "P27S41"),
                                      site46 = c("P27S42", NA, NA),
                                      site47 = c("P27S43", "P27S44", "P27S45"),
                                      site48 = c("P27S46", "P27S47", "P27S48"),
                                      site49 = c("P27S52", NA, NA),
                                      site50 = c("P28S4", "P28S5", "P28S6"),
                                      site51 = c("P30S1", "P30S2", "P30S3"),
                                      site52 = c("P31S1", "P31S2", "P31S3"),
                                      site53 = c("P32S4", "P32S5", "P32S6"),
                                      site54 = c("P33S1", "P33S2", "P33S3"),
                                      site55 = c("P34S5", "P34S6", NA),
                                      site56 = c("P35S2", "P35S3", NA),
                                      site57 = c("P36S4", "P36S5", "P36S6"),
                                      site58 = c("P36S10", "P36S11", "P36S12"),
                                      site59 = c("P36S13", "P36S14", "P36S15"),
                                      site60 = c("P36S22", "P36S23", "P36S24"),
                                      site61 = c("P36S28", "P36S29", "P36S30"),
                                      site62 = c("P36S31", "P36S32", "P36S33"),
                                      site63 = c("P36S37", "P36S38", "P36S39")
     )
     
     ## identify sites with 2 and 3 plants 
     lemna_sites_with_two_or_more = as.vector(which(apply(apply(lemna_same_location, 2, function(x) is.na(x)), 2, sum) != 2))
     
     ## remove singles
     lemna_same_location = lemna_same_location[,c(lemna_sites_with_two_or_more)]
     
     lemna_saver = list()   
     for (n in 1:ncol(lemna_same_location)) {
       
       runner_matrix = lemna_hamdist[which(colnames(lemna_hamdist) %in% lemna_same_location[,n]),which(colnames(lemna_hamdist) %in% lemna_same_location[,n]), drop=FALSE]
       lemna_saver[[n]] = runner_matrix[lower.tri(runner_matrix, diag=FALSE)]
       
     }    
     
     ## compute histograms with equal bins
     breaks = pretty(range(c(0,max(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)]))), n = 30)
     lemna_local_samples = hist(unlist(lemna_saver), breaks = breaks, plot = FALSE)
     lemna_global_samples = hist(lemna_hamdist[lower.tri(lemna_hamdist, diag=FALSE)], breaks = breaks, plot = FALSE)
     
     # assemble triplet plot
     plot(lemna_local_samples,
          col = scales::alpha("red", 0.5),
          xlab = "Genetic distance",
          ylab = "", main="",
          ylim = c(0, max(lemna_local_samples$counts)),
          xlim = c(0, max(lemna_hamdist+0.02)))
     axis(side=2, col.axis="red")
     
     ## add second histogram
     par(new = TRUE)
     plot(lemna_global_samples,
          col = scales::alpha("blue", 0.5),
          axes = FALSE, xlab = "", ylab = "", main = "",
          ylim = c(0, max(lemna_global_samples$counts)),
          xlim = c(0, max(lemna_hamdist+0.02)))
     axis(side = 4, col.axis="blue")
     
     
     
     ## Figure 4: Trees ####
     
     ## compute trees  
     landoltia_tree = bionj(landoltia_hamdist)
     landoltia_tree$edge.length[landoltia_tree$edge.length < 0] <- 0
     landoltia_chronotree = chronos(landoltia_tree,control = chronos.control(iter.max = 1000), model="relaxed")
     landoltia_chronotree$tip.label = substr(landoltia_chronotree$tip.label,1,3)
     
     lemna_tree = bionj(lemna_hamdist)
     lemna_tree$edge.length[lemna_tree$edge.length < 0] <- 0
     lemna_chronotree = chronos(lemna_tree,control = chronos.control(iter.max = 1000), model="relaxed")
     lemna_chronotree$tip.label = substr(lemna_chronotree$tip.label,1,3)
     
     ## plot trees
     
     par(mfrow=c(1,2))
     par(mar=c(0,0,0,0))
     plot.phylo(landoltia_chronotree, cex=0.6, type = "fan", tip.color=ifelse(grepl("P10", landoltia_chronotree$tip.label), "red",
                                                                              ifelse(grepl("P14", landoltia_chronotree$tip.label), "blue", 
                                                                                     ifelse(grepl("P19", landoltia_chronotree$tip.label), "darkgreen",
                                                                                            ifelse(grepl("P27", landoltia_chronotree$tip.label), "purple",
                                                                                                   ifelse(grepl("P36", landoltia_chronotree$tip.label), "orange", "black")))))
     )
     par(mar=c(0,0,0,0))
     plot.phylo(lemna_chronotree, cex=0.6, type = "fan", tip.color=ifelse(grepl("P10", lemna_chronotree$tip.label), "red",
                                                                          ifelse(grepl("P14", lemna_chronotree$tip.label), "blue", 
                                                                                 ifelse(grepl("P19", lemna_chronotree$tip.label), "darkgreen",
                                                                                        ifelse(grepl("P27", lemna_chronotree$tip.label), "purple",
                                                                                               ifelse(grepl("P36", lemna_chronotree$tip.label), "orange", "black")))))
                
     ) 
     
     
     
     
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
     
     
     ## LANDOLTIA compute euclidean distance ####
     
     ## euclidean distance
     landoltia_eucdist_raw = as.matrix(bitwise.dist(landoltia_final_genind, euclidean = TRUE))
     
     ## order euclidean distance matrix
     landoltia_eucdist_raw = landoltia_eucdist_raw[landoltia_ordered_names, landoltia_ordered_names]
     
     ## scale by maximum distance
     landoltia_eucdist = round(landoltia_eucdist_raw/(sqrt((nrow(landoltia_final@fix)*4))),5)
     
     
     ## LEMNA compute euclidean distance ####
     
     ## euclidean distance
     lemna_eucdist_raw = as.matrix(bitwise.dist(lemna_final_genind, euclidean = TRUE))
     
     ## order euclidean distance matrix
     lemna_eucdist_raw = lemna_eucdist_raw[lemna_ordered_names, lemna_ordered_names]
     
     ## scale by maximum distance
     lemna_eucdist = round(lemna_eucdist_raw/(sqrt((nrow(lemna_final@fix)*4))),5)
     
     
     ## 4x4 PCAs ####
     
     
     
     ## ade4::dudi.pca()
     par(mfrow=c(2,2))
     
     ## Landoltia with clones 
     dudipca_result = dudi.pca(landoltia_clone_no_missing_data, scannf = FALSE, nf = 3)
     plot(dudipca_result$li[,1], dudipca_result$li[,2], 
          xlim=c(-12,60), ylim=c(-35,15),
          xlab = paste0("PC1 (", round((dudipca_result$eig[1] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
          ylab = paste0("PC2 (", round((dudipca_result$eig[2] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
          main = "landoltia with clones (ade4)", pch = 21, 
          bg = ifelse(grepl("P10", rownames(dudipca_result$li)), P10_col,
                      ifelse(grepl("P14", rownames(dudipca_result$li)), P14_col, 
                             ifelse(grepl("P19", rownames(dudipca_result$li)), P19_col,
                                    ifelse(grepl("P27", rownames(dudipca_result$li)), P27_col,
                                           ifelse(grepl("P36", rownames(dudipca_result$li)), P36_col, rest_col))))))
     text(dudipca_result$li[,1], dudipca_result$li[,2],
          labels=rownames(dudipca_result$li), cex=0.5, pos=sample(1:4,nrow(dudipca_result$li), replace=TRUE),
          col="gray30")
     legend("bottomright", legend=c("P10","P14", "P19", "P27", "P36", "other"),
            pt.bg=c(P10_col,P14_col, P19_col, P27_col, P36_col, rest_col),
            pch=21, bty="n")
     
     ## landoltia without clones
     dudipca_result = dudi.pca(landoltia_no_clone_no_missing_data, scannf = FALSE, nf = 3)
     plot(dudipca_result$li[,1], dudipca_result$li[,2],
          xlab = paste0("PC1 (", round((dudipca_result$eig[1] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
          ylab = paste0("PC2 (", round((dudipca_result$eig[2] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
          main = "landoltia no clones (ade4)", pch = 21, 
          bg = ifelse(grepl("P10", rownames(dudipca_result$li)), P10_col,
                      ifelse(grepl("P14", rownames(dudipca_result$li)), P14_col, 
                             ifelse(grepl("P19", rownames(dudipca_result$li)), P19_col,
                                    ifelse(grepl("P27", rownames(dudipca_result$li)), P27_col,
                                           ifelse(grepl("P36", rownames(dudipca_result$li)), P36_col, rest_col))))))
     
     ## Lemna with clones
     dudipca_result = dudi.pca(lemna_clone_no_missing_data, scannf = FALSE, nf = 3)
     plot(dudipca_result$li[,1], dudipca_result$li[,2], 
          xlab = paste0("PC1 (", round((dudipca_result$eig[1] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
          ylab = paste0("PC2 (", round((dudipca_result$eig[2] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
          main = "lemna with clones (ade4)", pch = 21, 
          bg = ifelse(grepl("P10", rownames(dudipca_result$li)), P10_col,
                      ifelse(grepl("P14", rownames(dudipca_result$li)), P14_col, 
                             ifelse(grepl("P19", rownames(dudipca_result$li)), P19_col,
                                    ifelse(grepl("P27", rownames(dudipca_result$li)), P27_col,
                                           ifelse(grepl("P36", rownames(dudipca_result$li)), P36_col, rest_col))))))
     
     ## lemna without clones
     dudipca_result = dudi.pca(lemna_no_clone_no_missing_data, scannf = FALSE, nf = 3)
     plot(dudipca_result$li[,1], dudipca_result$li[,2], 
          xlab = paste0("PC1 (", round((dudipca_result$eig[1] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
          ylab = paste0("PC2 (", round((dudipca_result$eig[2] / sum(dudipca_result$eig)) * 100,2),"%)", sep=""), 
          main = "lemna no clones (ade4)", pch = 21, 
          bg = ifelse(grepl("P10", rownames(dudipca_result$li)), P10_col,
                      ifelse(grepl("P14", rownames(dudipca_result$li)), P14_col, 
                             ifelse(grepl("P19", rownames(dudipca_result$li)), P19_col,
                                    ifelse(grepl("P27", rownames(dudipca_result$li)), P27_col,
                                           ifelse(grepl("P36", rownames(dudipca_result$li)), P36_col, rest_col))))))
     
     
     
     ## adegenet::glPca()
     par(mfrow=c(2,2))
     
     ## landoltia with clones
     glpca_result = glPca(landoltia_final_pca, nf=3)
     plot(glpca_result$scores[,"PC1"], glpca_result$scores[,"PC2"],
          xlab = paste0("PC1 (", round((glpca_result$eig[1] / sum(glpca_result$eig)) * 100,2),"%)", sep=""), 
          ylab = paste0("PC1 (", round((glpca_result$eig[2] / sum(glpca_result$eig)) * 100,2),"%)", sep=""),
          main = "Landoltia with clones (adegenet)", pch = 21, 
          bg = ifelse(grepl("P10", rownames(glpca_result$scores)), P10_col,
                      ifelse(grepl("P14", rownames(glpca_result$scores)), P14_col, 
                             ifelse(grepl("P19", rownames(glpca_result$scores)), P19_col,
                                    ifelse(grepl("P27", rownames(glpca_result$scores)), P27_col,
                                           ifelse(grepl("P36", rownames(glpca_result$scores)), P36_col, rest_col))))))
     legend("topleft", legend=c("P10","P14", "P19", "P27", "P36", "other"),
            pt.bg=c(P10_col,P14_col, P19_col, P27_col, P36_col, rest_col),
            pch=21, bty="n")
     
     ## landoltia no clones
     glpca_result = glPca(landoltia_final_no_clone_pca, nf=3)
     plot(glpca_result$scores[,"PC1"], glpca_result$scores[,"PC2"],
          xlab = paste0("PC1 (", round((glpca_result$eig[1] / sum(glpca_result$eig)) * 100,2),"%)", sep=""), 
          ylab = paste0("PC1 (", round((glpca_result$eig[2] / sum(glpca_result$eig)) * 100,2),"%)", sep=""),
          main = "Landoltia no clones (adegenet)", pch = 21, 
          bg = ifelse(grepl("P10", rownames(glpca_result$scores)), P10_col,
                      ifelse(grepl("P14", rownames(glpca_result$scores)), P14_col, 
                             ifelse(grepl("P19", rownames(glpca_result$scores)), P19_col,
                                    ifelse(grepl("P27", rownames(glpca_result$scores)), P27_col,
                                           ifelse(grepl("P36", rownames(glpca_result$scores)), P36_col, rest_col))))))
     
     ## lemna with clones
     glpca_result = glPca(lemna_final_pca, nf=3)
     plot(glpca_result$scores[,"PC1"], glpca_result$scores[,"PC2"],
          xlab = paste0("PC1 (", round((glpca_result$eig[1] / sum(glpca_result$eig)) * 100,2),"%)", sep=""), 
          ylab = paste0("PC1 (", round((glpca_result$eig[2] / sum(glpca_result$eig)) * 100,2),"%)", sep=""),
          main = "lemna with clones (adegenet)", pch = 21, 
          bg = ifelse(grepl("P10", rownames(glpca_result$scores)), P10_col,
                      ifelse(grepl("P14", rownames(glpca_result$scores)), P14_col, 
                             ifelse(grepl("P19", rownames(glpca_result$scores)), P19_col,
                                    ifelse(grepl("P27", rownames(glpca_result$scores)), P27_col,
                                           ifelse(grepl("P36", rownames(glpca_result$scores)), P36_col, rest_col))))))
     
     ## lemna no clones
     glpca_result = glPca(lemna_final_no_clone_pca, nf=3)
     plot(glpca_result$scores[,"PC1"], glpca_result$scores[,"PC2"],
          xlab = paste0("PC1 (", round((glpca_result$eig[1] / sum(glpca_result$eig)) * 100,2),"%)", sep=""), 
          ylab = paste0("PC1 (", round((glpca_result$eig[2] / sum(glpca_result$eig)) * 100,2),"%)", sep=""),
          main = "lemna no clones (adegenet)", pch = 21, 
          bg = ifelse(grepl("P10", rownames(glpca_result$scores)), P10_col,
                      ifelse(grepl("P14", rownames(glpca_result$scores)), P14_col, 
                             ifelse(grepl("P19", rownames(glpca_result$scores)), P19_col,
                                    ifelse(grepl("P27", rownames(glpca_result$scores)), P27_col,
                                           ifelse(grepl("P36", rownames(glpca_result$scores)), P36_col, rest_col))))))
     
     
     
     
     draw_pie <- function(x0, y0, values, radius = 0.05, cols) {
       values <- values / sum(values)
       angles <- c(0, cumsum(values)) * 2 * pi
       
       for (i in seq_along(values)) {
         theta <- seq(angles[i], angles[i + 1], length.out = 100)
         x <- c(x0, x0 + radius * cos(theta), x0)
         y <- c(y0, y0 + radius * sin(theta), y0)
         
         polygon(x, y, col = cols[i], border = "black")
       }
     }
     
     plot(1:3, 1:3, type = "n", asp = 1)
     
     draw_pie(1, 1, c(2, 3, 4), cols = c("red", "blue", "green"))
     draw_pie(2, 2, c(5, 1, 2), cols = c("orange", "purple", "cyan"))
     draw_pie(3, 1.5, c(1, 1, 1), cols = c("black", "grey", "pink"))
     pie()
     
     
     
     
     theta <- seq(0, 2*pi, length.out = 200)
     x <- 0.5 + 0.1 * cos(theta)
     y <- 0.5 + 0.1 * sin(theta)
     polygon(x, y, col = c("red","blue"), border = "black")
     
     windows(type = "cairo")
     
     pdf("plot.pdf", width = 6, height = 6)
     plot(1, 1, type = "n")
     symbols(1, 1, circles = 0.1, inches = FALSE, add = TRUE)
     pie(c(2, 3, 4), radius = 0.1)
     dev.off()
     ## stacked barplot ####
     
     ## stacked barplot
     
     # # assembly plot and legend
     # barplot(stacked_barplotter, 
     #         col=c(lemna_col, landoltia_col, both_col), 
     #         space=0.1, ylab="Frequency", xlab="Plants per micro site", ylim=c(0,(max(colSums(stacked_barplotter))+1)))
     # box();
     # legend("topleft", inset=0.01, legend=c("lemna", "landoltia", "both"),
     #        fill=c(lemna_col, landoltia_col, both_col))
     