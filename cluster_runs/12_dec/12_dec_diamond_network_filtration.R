args = commandArgs(trailingOnly=TRUE)

# 1 - alpha
# 2 - iterations

if (!require("pacman")) install.packages("pacman", repos = getOption("repos"))
pacman::p_load(tidyverse, igraph, parallel)

seeds_10_weeks <- readRDS("seeds10weeks.rds")
STRING_mouse_700_removed_duplicates_extended <- readRDS("STRING_mouse_700_removed_duplicates_extended.rds")


# filtering genes for a given cell type (arg 2)

#===============================================================================
#                calculate connectivity significance p values 
#===============================================================================
# ks: number of connections to all seeds (initial & diamond), ks0 (number of connections to inital seeds),
# N: number of all the genes in the network
# s: numer of all seed genes (inital seeds & diamond genes that become seeds), s0: number of initial seed genes
# k: number of all connections (degree)
# alpha.: seed weight

get_CS <- function(ks, N, s, s0, k, alpha){ 
  m <- s + (alpha-1)*s0
  x <- ks
  n <- N - s
  K <- k 
  
  p_val <- phyper(q = x-1, m = m, n = n, k = K, lower.tail = FALSE)
  p_val
}

#================================================================================
#                       EXTENDED DIAMOND FUNCTION
#================================================================================

# function which chooses diamond genes not solely based on p value but also based on ks/k ratio (if the min p value is the same for multiple genes)
# seeds: vector of seed gene IDs, PPI: igraph object, iterations: number of diamond genes to add to the module
diamond_extended <- function(seeds, PPI, iterations, alpha) {
  
  
  initial_seeds <- seeds # a vector of inital seed genes that is always the same
  
  s0 <- length(initial_seeds) # number of inital seeds
  N <- gorder(PPI) # number of all the vertices in the network
  
  all_seeds <- initial_seeds  # after every iteration one diamond gene should be added to all seeds
  
  diamond_genes <- c()  #after every iteration one diamond gene should be added to all seeds.
  pvalues <- c() # after every iteration a pvalue of added diamond gene is added to this vector
  number_of_links_to_seeds <- c() # after every iteration a ks of added diamond gene is added to this vector
  number_of_all_links <- c() # after every iteration a k(degree) of added diamond gene is added to this vector
  
 
  while(length(diamond_genes) < iterations){
    
    s <- length(all_seeds) # number of all seeds
    
    # find neighbors of all seeds
    neighbors_all <- igraph::adjacent_vertices(graph = PPI, v = all_seeds)
    neighbors_all <- unlist(sapply(X = neighbors_all, FUN = as_ids, simplify = TRUE, USE.NAMES = TRUE))
    neighbors_all <- unique(neighbors_all)
    
    # exclude from neighbors_all all the genes that are already in the module:
    neighbors  <- neighbors_all[!neighbors_all %in% all_seeds]
    
    # find neighbors of neighbors
    
    neighbors_of_neighbors<- igraph::adjacent_vertices(graph = PPI, v = neighbors)
    neighbors_of_neighbors<- sapply(X = neighbors_of_neighbors, FUN = as_ids, simplify = TRUE, USE.NAMES = TRUE)
    
    # count how many neighbors of each neighbor are seeds (inital seeds and all seeds) --> make ks and ks0 vecs
    
    ks_of_neighbors <- sapply(neighbors_of_neighbors, function(x){
      kb <- length(dplyr::intersect(x, all_seeds))
      return(kb)
    }
    , simplify = TRUE, USE.NAMES =TRUE)
    
    ks0_of_neighbors <- sapply(neighbors_of_neighbors, function(x){
      kb <- length(dplyr::intersect(x, initial_seeds))
      return(kb)
    }
    , simplify = TRUE, USE.NAMES =TRUE)
    
  
    # get degree of all neighbors:
    
    neighbors_degrees <- degree(graph = PPI, v = neighbors)
    
    # add weights to number of connections to seeds & total number of connections
    
    ks_final <- ks_of_neighbors + (alpha-1)*ks0_of_neighbors
    k_final <- neighbors_degrees +(alpha-1)*ks0_of_neighbors
    
    #reduce the number of genes for which the p value calculation will be performed
    ks_k_df <- tibble(gene = names(ks_final), ks = ks_final, k = k_final) %>%
               # classify the nodes based on their ks and rank the node with lowest k highest within that class
               group_by(ks) %>% 
               filter(k == min(k)) %>% 
               # classify top ranks of each class by their degree k and choose the ones with highest ks
               ungroup() %>% 
               group_by(k) %>%
               filter(ks == max(ks))
    
    # calculate connectivity significance p value for genes that remain in ks_k_df
    
    remaining_genes_ks <- ks_k_df$ks
    names(remaining_genes_ks) <- ks_k_df$gene
    
    remaining_genes_k <- ks_k_df$k
    names(remaining_genes_k) <- ks_k_df$gene
    
 
    CS_pvals <- mapply(ks = remaining_genes_ks, k = remaining_genes_k, N = N,
                                                s = s, s0 = s0, alpha = alpha, FUN = get_CS, SIMPLIFY = TRUE, USE.NAMES = TRUE)
    
    if(length(which(CS_pvals == min(CS_pvals))) == 1){ # if there is only one gene with the min pvalue, include this gene into the module
      
      new_gene <- names(CS_pvals)[which(CS_pvals == min(CS_pvals))]
      new_gene_ks <- remaining_genes_ks[new_gene]
      new_gene_k <- remaining_genes_k[new_gene]
      p_value <- CS_pvals[new_gene]
    

    }else{# if there are more genes with min pvalue, get ks/k ratio for those genes and include the gene with highest ks/k ratio in the module
      
      equal_genes <- names(CS_pvals)[which(CS_pvals == min(CS_pvals))]
      ks_equal_genes <- remaining_genes_ks[equal_genes]
      k_equal_genes <- remaining_genes_k[equal_genes]
      
      ks_k_ratio <- ks_equal_genes/k_equal_genes
      
      new_gene <- names(ks_k_ratio)[which(ks_k_ratio == max(ks_k_ratio))]
      new_gene_ks <- remaining_genes_ks[new_gene]
      new_gene_k <- remaining_genes_k[new_gene] 
      p_value <- CS_pvals[new_gene]
      
    }

    # add new gene to the list of diamond genes ( and corresponding CS p-value, ks and k)
    
    diamond_genes <- c(diamond_genes, new_gene)
    pvalues <- c(pvalues, p_value)
    number_of_links_to_seeds <- c(number_of_links_to_seeds, new_gene_ks) 
    number_of_all_links <- c(number_of_all_links, new_gene_k) 
    
    
    # add new gene to the list of seed genes
    
    all_seeds <- c(all_seeds, new_gene)
    
    
  } # end of while loop
  
  added_genes_df <- data.frame(gene = diamond_genes,
                              degree = number_of_all_links,
                              connectivity = number_of_links_to_seeds,
                              pvalue = pvalues)
  
  module_genes_vec <- all_seeds
  
  results_list <- list(module_genes = module_genes_vec,
                       added_genes = added_genes_df)
  
  
  return(results_list)
  
  
}

STRING_mouse_700_vertices_extended <- STRING_mouse_700_removed_duplicates_extended %>%
        select(protein1, protein2) %>%
        flatten_chr() %>%
        unique()

string_igraph <- graph_from_data_frame(STRING_mouse_700_removed_duplicates_extended, directed = FALSE, vertices = STRING_mouse_700_vertices_extended)

# filtering of genes with p-values
RML_data_26_May_ENTREZ_NAs_filtered <- readRDS("RML_data_26_May_ENTREZ_NAs_filtered.rds")

# get ENTREZ vecctors

cell_type_genes <- sapply(c("vGluT2", "Gad2", "Cx43"), function(x){
  RML_data_26_May_ENTREZ_NAs_filtered$ENTREZ[which(RML_data_26_May_ENTREZ_NAs_filtered$Cell_Type == x & RML_data_26_May_ENTREZ_NAs_filtered$Timepoint == "18 Weeks")]
}
, USE.NAMES = TRUE, simplify = TRUE
)

string_igraphs_filtered <- lapply(cell_type_genes, function(x){
  igraph::induced_subgraph(graph = string_igraph, vids = intersect(x, STRING_mouse_700_vertices_extended), impl = "auto")
} 
)


run_diamond <- function(iteration_no){
                result <- mcmapply(
                                 FUN = diamond_extended,
                                 seeds = seeds_10_weeks[1:3],
                                 PPI = string_igraphs_filtered,
                                 MoreArgs = list(
                                         iterations = iteration_no,
                                         alpha = as.numeric(args[1])),
                                 SIMPLIFY = FALSE, USE.NAMES = TRUE,
                                 mc.preschedule = TRUE, mc.set.seed = TRUE,
                                 mc.silent = FALSE, mc.cores = getOption("mc.cores", 3L),
                                 mc.cleanup = TRUE, affinity.list = NULL)
                
                return(result)
}

b <- Sys.time()
d <- run_diamond(as.numeric(args[2]))
a <- Sys.time()
a-b
saveRDS(d, paste0("all_cells_diamond_",args[2],"x_alpha_",args[1],".rds"))

# just test
# b <- Sys.time()
# d <- run_diamond(as.numeric(3))
# a <- Sys.time()
# a-b
# 
# run_diamond <- function(iteration_no){
#   result <- mcmapply(
#     FUN = diamond_extended,
#     seeds = seeds_10_weeks[1:3],
#     PPI = string_igraphs_filtered,
#     MoreArgs = list(
#       iterations = iteration_no,
#       alpha = as.numeric(10)),
#     SIMPLIFY = FALSE, USE.NAMES = TRUE,
#     mc.preschedule = TRUE, mc.set.seed = TRUE,
#     mc.silent = FALSE, mc.cores = getOption("mc.cores", 3L),
#     mc.cleanup = TRUE, affinity.list = NULL)
#   
#   return(result)
# }