suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(dendextend))


sankoff <- function(hc, merge_values, leaf_values) {
    
    # merge values counting the numbers of mutations are initially all 0
    
    # iterate through all inner nodes of the dendrogram
    for(i in 1:nrow(hc$merge)) { 
        
        # initialize the values for the left and right child nodes of node i with Infinity
        lval <- c(Inf,Inf)
        rval <- c(Inf,Inf)
        
        # if the left child node is a leaf, set the left child value to the mutation numbers computed for this leaf
        if(hc$merge[i,1] < 0) {
            lval <- leaf_values[-hc$merge[i,1],]
        }
        # if the left child node is an inner node, set the left child value to the mutation numbers computed for this node
        else {
            lval <- merge_values[hc$merge[i,1],]
        }
        
        # if the right child node is a leaf, set the right child value to the mutation numbers computed for this leaf
        if(hc$merge[i,2] < 0) {
            rval <- leaf_values[-hc$merge[i,2],]
        } 
        # if the right child node is an inner node, set the right child value to the mutation numbers computed for this node
        else {
            rval <- merge_values[hc$merge[i,2],]
        }
        new_val <- c(0, 0)
        new_val[1] <- min(lval[1], lval[2] +1) + min(rval[1], rval[2] + 1)
        new_val[2] <- min(lval[1] + 1, lval[2]) + min(rval[1] + 1, rval[2])
        merge_values[i,] <- new_val
    }
    return(merge_values)
}

# algorithm counting the numbers of mutations that are necessary at minimum in each node of the dendrogram using the Sankoff algorithm
computeSankoffMerge <- function(hc, geno_matrix, first_position, second_position) {
    # because there could be more than just 0 and 1, but it's always bi-allelic
    # so I'm just taking the first value as 0 and the second as 1
    first_pos_genotype <- unique(geno_matrix[,as.character(first_position)])
    first_pos_genotype <- sort(first_pos_genotype)
    second_pos_genotype <- unique(geno_matrix[,as.character(second_position)])
    second_pos_genotype <- sort(second_pos_genotype)
    # extract the genotypes of the TF and TFBS from the genotype matrix
    tf_snvs <- geno_matrix[,as.character(first_position)]
    tfbs_snvs <- geno_matrix[,as.character(second_position)]
    
    # initialize the leaf mutation numbers to Infinity for both genotypes
    # the first genotype is 0 and the second genotype is 1
    leaf_values_tf <- matrix(Inf, nrow = length(tf_snvs), ncol = 2)
    leaf_values_tfbs <- matrix(Inf, nrow = length(tfbs_snvs), ncol = 2)
    
    # for the TF variant initialize the leaf mutation numbers for both genotypes
    for(c in 1:length(tf_snvs)) {
        # initialize both values with 0 if the leaf genotype is NA
        if(is.na(tf_snvs[c])) {
            print("WARNING: NA value in TF")
            leaf_values_tf[c,] = c(0,0)
        }
        else if(tf_snvs[c] == first_pos_genotype[1]) {
            leaf_values_tf[c,1] = 0
        } else if(tf_snvs[c] == first_pos_genotype[2]) {
            leaf_values_tf[c,2] = 0
        }
    }
    
    # for the TFBS variant initialize the leaf mutation numbers for both genotypes
    for(c in 1:length(tfbs_snvs)) {
        # initialize both values with 0 if the leaf genotype is NA
        if(is.na(tfbs_snvs[c])) {
            print("WARNING: NA value in TFBS")
            leaf_values_tfbs[c,] = c(0,0)
        }
        else if(tfbs_snvs[c] == second_pos_genotype[1]) {
            leaf_values_tfbs[c,1] = 0
        } else if(tfbs_snvs[c] == second_pos_genotype[2]) {
            leaf_values_tfbs[c,2] = 0
        }
    }
    
    # the first genotype is 0 and the second genotype is 1
    # construct matrices storing the mutation numbers computed with the Sankoff algorithm for the TF variant and the TFBS variant
    merge_values_tf <- matrix(0, nrow = nrow(hc$merge), ncol = 2)
    merge_values_tfbs <- matrix(0, nrow = nrow(hc$merge), ncol = 2)
    
    # apply the Sankoff algorithm to reconstruct the numbers of mutations necessary in each node
    merge_values_tf <- sankoff(hc, merge_values_tf, leaf_values_tf)
    merge_values_tfbs <- sankoff(hc, merge_values_tfbs, leaf_values_tfbs)
    
    # print(cbind(merge_values_tf, merge_values_tfbs))
    return(cbind(merge_values_tf, merge_values_tfbs))
    
}

# perform hierarchical clustering on genetic variants and label the leaves with the genotypes of the TF and the TFBS variant pair 
hierarchicalClustering <- function(geno_matrix, method, first_position, second_position, visualize) {
    
    first_pos_genotype <- unique(geno_matrix[,as.character(first_position)])
    first_pos_genotype <- sort(first_pos_genotype)
    second_pos_genotype <- unique(geno_matrix[,as.character(second_position)])
    second_pos_genotype <- sort(second_pos_genotype)
    # extract the genotypes of the TF and TFBS from the genotype matrix
    tf_labels <- geno_matrix[,as.character(first_position)]
    tfbs_labels <- geno_matrix[,as.character(second_position)]
    
    # perform hierarchical clustering based on the given genetic variants
    hc <- hclust(dist(geno_matrix), method=method)
    
    # convert genotypes: 0 to o and 1 to x
    tf_labels[tf_labels == first_pos_genotype[1]] <- 'o'
    tf_labels[tf_labels == first_pos_genotype[2]] <- 'x'
    tf_labels[is.na(tf_labels)] <- 'N'
    
    tfbs_labels[tfbs_labels == second_pos_genotype[1]] <- 'o'
    tfbs_labels[tfbs_labels == second_pos_genotype[2]] <- 'x'
    tfbs_labels[is.na(tfbs_labels)] <- 'N'
    
    # combine the labels for the TF variant and the TFBS variant
    hc$labels <- paste(tf_labels, tfbs_labels, sep=" ")
    
    # visualize dendrogram if chosen by the user
    if(visualize) {
        # assign colors to labels
        colors <- rep('black', 220)
        colors[hc$labels == 'x o'] <- 'blue'
        colors[hc$labels == 'o x'] <- 'yellow'
        colors[hc$labels == 'o o'] <- 'green'
        colors[hc$labels == 'x x'] <- 'red'
        
        # convert result of the hierarchical clustering to a dendrogram that can be visualized
        dend <- as.dendrogram(hc)
        # color the leafs of the dendrogram according to their genotypes
        labels_colors(dend) <- colors[order.dendrogram(dend)]
        
        # save visualizations of thIe constructed dendrograms to user-specified directory
        
        output_name = paste("pair_",first_position, "_",second_position, format(Sys.time(), "_%s"), ".pdf", sep='')
        pdf(output_name)
        plot(dend)
        dev.off()
    }
    
    # return the constructed dendrogram
    return(hc)
}


# algorithm counting the number of same direction mutations and opposite direction mutations in a dendrogram
# look how TF variant changes with TFBS seen as constant, can be applied for constant TF and changing TFBS by exchanging TF and TFBS
countMutations <- function(geno_matrix, TF_position, TFBS_position, hc, merge_mutations_tf, merge_mutations_tfbs) {
  
  # number of same direction mutations
  sym_tf <- 0
  
  # number of opposite direction mutations
  assym_tf <- 0
  
  # extract the genotypes of the TF and the TFBS variants
  tf_snvs <- geno_matrix[,as.character(TF_position)]
  tfbs_snvs <- geno_matrix[,as.character(TFBS_position)]
  
  # iterate through the dendrogram, starting at the root
  for(i in nrow(hc$merge):1) {
    
    # get the left and right child nodes of the current node i
    left <- hc$merge[i,1]
    right <- hc$merge[i,2]
    
    # the left child node is an inner node
    if(left > 0) {
      # the genotype of the TF variant is different in the left child node than in the parent node, meaning a mutation takes place
      if(merge_mutations_tf[i] != merge_mutations_tf[left]) {
        # if the genotypes of the TF variant and the TFBS variant are the same after the mutation, count a same direction mutation, else count an opposite direction mutation
        if(merge_mutations_tf[left] == merge_mutations_tfbs[left]) {
          sym_tf = sym_tf + 1
        } else {
          assym_tf = assym_tf + 1
        }
      }
    } 
    # the left child node is a leaf
    else {
      # cannot count mutation direction if the resulting genotypes are unknown
      if(is.na(tf_snvs[-left]) | is.na(tfbs_snvs[-left])) {
        next
      }
      # the genotype of the TF variant is different in the left child node than in the parent node, meaning a mutation takes place
      if(merge_mutations_tf[i] != tf_snvs[-left]) {
        # if the genotypes of the TF variant and the TFBS variant are the same after the mutation, count a same direction mutation, else count an opposite direction mutation
        if(tf_snvs[-left] == tfbs_snvs[-left]) {
          sym_tf = sym_tf + 1
        } else {
          assym_tf = assym_tf + 1
        }
      }
    }
    
    # the right child node is an inner node
    if(right > 0) {
      # the genotype of the TF variant is different in the right child node than in the parent node, meaning a mutation takes place
      if(merge_mutations_tf[i] != merge_mutations_tf[right]) {
        # if the genotypes of the TF variant and the TFBS variant are the same after the mutation, count a same direction mutation, else count an opposite direction mutation
        if(merge_mutations_tf[right] == merge_mutations_tfbs[right]) {
          sym_tf = sym_tf + 1
        } else {
          assym_tf = assym_tf + 1
        }
      }
    } 
    # the right child node is a leaf
    else {
      # cannot count mutation direction if the resulting genotypes are unknown
      if(is.na(tf_snvs[-right]) | is.na(tfbs_snvs[-right])) {
        next
      }
      # the genotype of the TF variant is different in the right child node than in the parent node, meaning a mutation takes place
      if(merge_mutations_tf[i] != tf_snvs[-right]) {
        # if the genotypes of the TF variant and the TFBS variant are the same after the mutation, count a same direction mutation, else count an opposite direction mutation
        if(tf_snvs[-right] == tfbs_snvs[-right]) {
          sym_tf = sym_tf + 1
        } else {
          assym_tf = assym_tf + 1
        }
      }
    }
  } 
  
  # return the counted numbers of same and opposite direction mutations
  return(c(sym_tf, assym_tf))
}


dendrogramMutation <- function(geno_matrix, TF_position, TFBS_position, hc, merge_values_tf, merge_values_tfbs, binomial_tail) {
    
    first_pos_genotype <- unique(geno_matrix[,as.character(TF_position)])
    first_pos_genotype <- sort(first_pos_genotype)
    second_pos_genotype <- unique(geno_matrix[,as.character(TFBS_position)])
    second_pos_genotype <- sort(second_pos_genotype)
    
    merge_mutations_tf <- c()
    # reconstruct the most likely genotypes (the genotypes that require the least mutations) 
    # for all inner nodes of the TF variant dendrogram based on the numbers of mutations (merge values) computed by the Sankoff algorithm
    for(i in nrow(merge_values_tf):1) {
        if(merge_values_tf[i,1] > merge_values_tf[i,2]) {
            merge_mutations_tf[i] <- first_pos_genotype[2]
        } else if(merge_values_tf[i,1] < merge_values_tf[i,2]) {
            merge_mutations_tf[i] <- first_pos_genotype[1]
        }
        # if both genotypes require the same number of mutations, 
        # take the first genotype if the considered node is the root of the dendrogram and take the genotype of the parent node otherwise
        else {
            if(i == nrow(merge_values_tf)) {
                merge_mutations_tf[i] <- first_pos_genotype[1]
            } else {
                # same as parent
                parent <- which(hc$merge == i, TRUE)[1]
                merge_mutations_tf[i] <- merge_mutations_tf[parent]
            }
        }
    } 
    # print(merge_mutations_tf)
    merge_mutations_tfbs <- c()
    
    # reconstruct the most likely genotypes (the genotypes that require the least mutations) 
    # for all iner nodes of the TFBS variant dendrogram based on the numbers of mutations (merge values) computed by the Sankoff algorithm
    for(i in nrow(merge_values_tfbs):1) {
        if(merge_values_tfbs[i,1] > merge_values_tfbs[i,2]) {
            merge_mutations_tfbs[i] <- second_pos_genotype[2]
        } else if(merge_values_tfbs[i,1] < merge_values_tfbs[i,2]) {
            merge_mutations_tfbs[i] <- second_pos_genotype[1]
        } 
        # if both genotypes require the same number of mutations, take the first genotype if the considered node is the root of the dendrogram and take the genotype of the parent node otherwise
        else {
            if(i == nrow(merge_values_tfbs)) {
                merge_mutations_tfbs[i] <- second_pos_genotype[1]
            } else {
                #same as parent
                parent <- which(hc$merge == i, TRUE)[1]
                merge_mutations_tfbs[i] <- merge_mutations_tfbs[parent]
            }
        }
    } 
    
    # print(merge_mutations_tfbs)
    # look how tf changes with tfbs seen as constant
    mutation_count_tf <- countMutations(geno_matrix, TF_position, TFBS_position, hc, merge_mutations_tf, merge_mutations_tfbs)
    print("The mutation count of tf")
    print(mutation_count_tf)
    # look how tfbs changes with tf seen as constant
    mutation_count_tfbs <- countMutations(geno_matrix, TFBS_position, TF_position, hc, merge_mutations_tfbs, merge_mutations_tf)
    print("The mutation count of tfbs")
    print(mutation_count_tfbs)
    # apply right-tailed, two-tailed, or left-tailed binomial tests to the counted numbers of same and opposite direction mutations for the TF and the TFBS
    binom_tf <- NULL
    binom_tfbs <- NULL
    if(binomial_tail == "right-tailed") {
        # print("calculating binom.test")
        binom_tf <- binom.test(mutation_count_tf, 0.5, alternative='greater')
        binom_tfbs <- binom.test(mutation_count_tfbs, 0.5, alternative='greater')
    } else if(binomial_tail == "two-tailed"){
        binom_tf <- binom.test(mutation_count_tf, 0.5)
        binom_tfbs <- binom.test(mutation_count_tfbs, 0.5)
    } else if(binomial_tail == "left-tailed") {
        binom_tf <- binom.test(mutation_count_tf, 0.5, alternative='less')
        binom_tfbs <- binom.test(mutation_count_tfbs, 0.5, alternative='less')
    }
    
    # extract the p-values for the TF and the TFBS from the results of the binomial tests
    pval_tf <- binom_tf$p.value
    pval_tfbs <- binom_tfbs$p.value
    
    
    return(c(pval_tf,pval_tfbs))
}


# find combinations between two positions in the aligned proteins
# hardcoded as it was faster thanto write 

combinations <- function(first_position_aa, second_position_aa){
    if ((length(first_position_aa) == 2) & (length(second_position_aa) == 2)){
        return(list(list(first_position_aa, second_position_aa)))
    }
    
    else if((length(first_position_aa) == 2) & (length(second_position_aa) == 3)){
        return(list(list(first_position_aa, c(second_position_aa[1], second_position_aa[2])),
                    list(first_position_aa, c(second_position_aa[1], second_position_aa[3]))))
    }
    
    else if((length(first_position_aa) == 2) & (length(second_position_aa) >= 4)){
        return(list(list(first_position_aa, c(second_position_aa[1], second_position_aa[2])),
                    list(first_position_aa, c(second_position_aa[1], second_position_aa[3])),
                    list(first_position_aa, c(second_position_aa[1], second_position_aa[4]))))
    }
    
    else if((length(first_position_aa) == 3) & (length(second_position_aa) == 2)){
        return(list(list(c(first_position_aa[1], first_position_aa[2]), second_position_aa),
                    list(c(first_position_aa[1], first_position_aa[3]), second_position_aa)))
    }
    
    else if((length(first_position_aa) == 3) & (length(second_position_aa) == 3)){
        return(list(list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[3])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[3]))))
    }
    
    
    else if((length(first_position_aa) == 3) & (length(second_position_aa) >= 4)){
        return(list(list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[3])),
                    list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[4])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[3])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[4]))))
    }
    
    else if ((length(first_position_aa) >= 4) & (length(second_position_aa) == 2)){
        return(list(list(c(first_position_aa[1], first_position_aa[2]), second_position_aa),
                    list(c(first_position_aa[1], first_position_aa[3]), second_position_aa),
                    list(c(first_position_aa[1], first_position_aa[4]), second_position_aa)))
    }
    
    
    else if((length(first_position_aa) >= 4) & (length(second_position_aa) == 3)){
        return(list(list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[3])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[3])),
                    list(c(first_position_aa[1], first_position_aa[4]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[4]), c(second_position_aa[1], second_position_aa[3]))))
    }
    
    else if((length(first_position_aa) >= 4) & (length(second_position_aa) >= 4)){
        return(list(list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[3])),
                    list(c(first_position_aa[1], first_position_aa[2]), c(second_position_aa[1], second_position_aa[4])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[3])),
                    list(c(first_position_aa[1], first_position_aa[3]), c(second_position_aa[1], second_position_aa[4])),
                    list(c(first_position_aa[1], first_position_aa[4]), c(second_position_aa[1], second_position_aa[2])),
                    list(c(first_position_aa[1], first_position_aa[4]), c(second_position_aa[1], second_position_aa[3])),
                    list(c(first_position_aa[1], first_position_aa[4]), c(second_position_aa[1], second_position_aa[4]))))
    }
    
}

############################################# 
# program options
option_list = list(
    make_option(c("-v", "--vcf"), type="character", default=NULL, 
                help="VCF file containing all genetic variants", metavar="character"),
    # make_option(c("-f", "--functions"), type="character", default=NULL,
    #               help="Path to the functions script"),
    make_option(c("-p", "--pairs"), type="character", default=NULL, 
                help="tab-delimited file containing the positions of a variant pair (one tab-delimited variant pair per line)"),
    make_option(c("-t", "--sequences_table"), type="character", default=NULL,
                help="The protein sequences TSV generated with the python code"),
    make_option(c("-o", "--output"), type="character", default="output.tsv",
                help="You can specify the output table name"),
    make_option(c("-b", "--binomial"), type="character", default="right-tailed", 
              help="should the binomial test be two-tailed, left-tailed, or right-tailed (default: right-tailed)", metavar="character"),
    make_option(c("-s", "--start"), type="integer", default=1,
                help="Where to start calculating in the pairs table (in case there are too many pairs and you need run in parallel)"),
    make_option(c("-l", "--slice"), type="integer", default=0,
                help="How many pairs to take after the start position")
); 
# parse program parameters from command line
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

method <- 'ward.D2'

if(is.null(opt$vcf)){
    print_help(opt_parser)
    stop("A VCF file needs to be supplied.", call.=FALSE)
}

# check if a tab-delimited file containing the variant pair positions has been supplied by the user
if(is.null(opt$pairs)){
    print_help(opt_parser)
    stop("A file containing variant pairs needs to be supplied.", call.=FALSE)
}

# check if an output file path has been supplied
if(is.null(opt$sequences_table)){
    print_help(opt_parser)
    stop("The protein sequences table needs to be provided", call.=FALSE)
}

# my VCF file
protein_vcf <- read.table(opt$vcf, header = T, sep = '\t')

# reading the protein table, it usually takes T as true and F as false, so I need the last two arguments
proteins <- read.table(opt$sequences_table, header = T, sep = '\t', stringsAsFactors=FALSE, colClasses = c("character"),
                       check.names = F)


possible_pairs <- read.table(opt$pairs, header=T, sep='\t')

protein_vcf$gt_GT_allele <- as.character(protein_vcf$gt_GT_allele)

# I'm using all the VCF to generate one geno matrix
protein_geno_matrix <- acast(protein_vcf, Indiv ~ POS,
                             value.var = 'gt_GT', fill = NaN,
                             fun.aggregate=function(x) as.integer(median(x,na.rm=TRUE)))


# testing two positions
# protein_geno_matrix <- acast (protein_vcf[protein_vcf$POS == 2 | protein_vcf$POS == 10,],
#                               Indiv ~ POS, value.var = 'gt_GT', fill = NaN, fun.aggregate = function(x) as.integer(median(x, na.rm=TRUE)))
stopifnot(nrow(proteins) == nrow(protein_geno_matrix))
protein_geno_matrix <- protein_geno_matrix[order(rownames(protein_geno_matrix)),]

rownames(proteins) <- proteins[,"samples_name"]
# I need to remove this later once I remove the sample name column from the python code
proteins <- proteins[,-1]

proteins <- proteins[order(rownames(proteins)),]
stopifnot(rownames(proteins) == rownames(protein_geno_matrix))


mutation_pvalues <- matrix(0,nrow=nrow(possible_pairs),ncol=2)


# just a random first row to initialize the table that I will remove later
output_df <- data.frame("pos1"=0,"pos2"=0, "a11"="A", "a12"="A", "a21"="A",
                        "a22"="A", "num_samples"=10, "num_samples11"=10,
                        "num_samples12"=10, "num_samples21"=10, "num_samples22"=10,
                        "pval1"="0.4", "pval2"="0.4", stringsAsFactors=FALSE)

# I'm splitting the pairs to run every 1000 pairs in parallel to finish faster
start_of_loop = as.integer(opt$start)
slice = as.integer(opt$slice)

if (start_of_loop == 1 & slice == 0){
    end_of_loop = nrow(possible_pairs)

} else if (start_of_loop > (nrow(possible_pairs) - slice)){
    end_of_loop = nrow(possible_pairs)

} else {
    end_of_loop = start_of_loop + slice
}

start_print("starting the loop")
for(i in c(start_of_loop:end_of_loop)) {
# for (i in c(1, 10, 530, 200)){
    
    first_position <- possible_pairs[i,1]
    second_position <- possible_pairs[i,2]

    first_position_aa <- unique(protein_vcf[protein_vcf$POS == first_position,][,6])
    first_position_aa <- strsplit(first_position_aa, ',')[[1]]
    second_position_aa <- unique(protein_vcf[protein_vcf$POS == second_position,][,6])
    second_position_aa <- strsplit(second_position_aa, ',')[[1]]

    combinations_list <- combinations(first_position_aa, second_position_aa)
    # print(paste("at ", first_position, " and ", second_position))
    # second loop for possible combinations between the same positions in case they had more than 2 possible aa for that position
    # the combinations are needed because we want two letters to have a 0,1 value when counting mutations
    # so if first or second positions had more than 2 proteins, I take combinations
    for(j in c(1:length(combinations_list))){
        # This will return all possible combinations between the proteins in two positions
        first_position_aa <- combinations_list[[j]][[1]]
        second_position_aa <- combinations_list[[j]][[2]]
        # print(paste(first_position_aa, second_position_aa))
        
        # some filtering here based on the letters
        # (one) is a (true,false) vectore for which samples has the two AAs I'm looking for or not
        # (two) is the same but for the second position
        # then I'm taking the intersection (only samples that has both combination) for the current matrix
        # The column names in proteins is the same as python indexing (0 indexing) because all positions
        # are produced by the python preprocessing code, so when it says position 2 it's position 3 if we're counting from 0
        one <- proteins[,as.character(first_position)] == first_position_aa[1] | proteins[,as.character(first_position)] == first_position_aa[2]
        two <- proteins[,as.character(second_position)] == second_position_aa[1] | proteins[,as.character(second_position)] == second_position_aa[2]
        
        current_matrix <- protein_geno_matrix[one & two,]
        if (nrow(current_matrix) == 0){
            next
        }
        
        hc <- hierarchicalClustering(current_matrix, method, first_position, second_position, F)
        # count the numbers of mutations that are necessary at minimum in each node of the dendrogram using the Sankoff algorithm and reconstruct the most likely genotypes of the inner nodes
        merge_values <- computeSankoffMerge(hc, current_matrix, first_position, second_position)
        
        # count the numbers of same and opposite direction mutations for both the TF and the TFBS variant within the dendrogram and compute p-values
        pval <- dendrogramMutation(current_matrix, first_position, second_position, hc,
                                   merge_values[,1:2], merge_values[,3:4], opt$binomial)
        
        # store the computed p-values in the result matrix
        # mutation_pvalues[i,] <- pval
        sample_num11 <- proteins[as.character(first_position)] == first_position_aa[1] & proteins[as.character(second_position)] == second_position_aa[1]
        # nsample_num11 <- as.double(nrow(protein_geno_matrix[sample_num11,]))
        nsample_num11 <- length(which(sample_num11))
        
        sample_num12 <- proteins[as.character(first_position)] == first_position_aa[1] & proteins[as.character(second_position)] == second_position_aa[2]
        # nsample_num12 <- as.double(nrow(protein_geno_matrix[sample_num12,]))
        nsample_num12 <- length(which(sample_num12))

        sample_num21 <- proteins[as.character(first_position)] == first_position_aa[2] & proteins[as.character(second_position)] == second_position_aa[1]
        # nsample_num21 <- as.double(nrow(protein_geno_matrix[sample_num21,]))
        nsample_num21 <- length(which(sample_num21))

        sample_num22 <- proteins[as.character(first_position)] == first_position_aa[2] & proteins[as.character(second_position)] == second_position_aa[2]
        # nsample_num22 <- as.double(nrow(protein_geno_matrix[sample_num22,]))
        nsample_num22 <- length(which(sample_num22))
        # print(sample_num22)
        # print(which(sample_num22))
        # print(nsample_num22)
        # print(dim(sample_num22))
        # print(first_position_aa)
        # print(first_position)
        # print(second_position_aa)
        # print(second_position)
        # if (nsample_num22 == 0){
        #     print(sample_num22)
        # }
        
        # filling the output data frame
        output_df[nrow(output_df) + 1, ] <- list(first_position, second_position,
                                                 first_position_aa[1], first_position_aa[2],
                                                 second_position_aa[1], second_position_aa[2],
                                                 as.double(dim(current_matrix)[1]),
                                                 nsample_num11, nsample_num12, nsample_num21, nsample_num22,
                                                 format(pval[1], scientific = FALSE),
                                                 format(pval[2], scientific = FALSE))
    }

}
# format(pval[1], scientific=FALSE)
output_df <- output_df[-1,]
write.table(output_df, opt$output ,row.names=FALSE, sep="\t", quote = FALSE)
