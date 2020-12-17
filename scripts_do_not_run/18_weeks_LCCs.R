if (!require("pacman")) install.packages("pacman", repos = getOption("repos"))
pacman::p_load(tidyverse, igraph, matrixStats, parallel)

# file.create("log_randomization_11_october.txt")
# sink(file = "log_randomization_11_october.txt")

STRING_mouse_700_removed_duplicates <- readRDS("STRING_mouse_700_removed_duplicates.rds")
RML_data_26_May_ENTREZ_NAs_filtered <- readRDS("RML_data_26_May_ENTREZ_NAs_filtered.rds")

#extract all nodes from STRING PPI network:

STRING_mouse_700_vertices <- STRING_mouse_700_removed_duplicates %>%
        flatten_chr() %>%
        unique()

RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE <- gsub(" ", "", RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE)

RML_data_filtered_network_genes <- RML_data_26_May_ENTREZ_NAs_filtered %>%
        dplyr::filter(ENTREZ %in% STRING_mouse_700_vertices, !grepl("mt-", Gene_Symbol)) %>%
        dplyr::group_by(Comparison_DE)


#function that calculates LCC size of 1 random graph:
fun <- function(x){
        
        n_seed = x
        vertices <- sample(STRING_mouse_700_vertices, size = n_seed, replace = FALSE) #select n_seed random vertices from STRING network
        
        edges <- STRING_mouse_700_removed_duplicates %>%
                filter(protein1 %in% vertices & protein2 %in% vertices)
        
        igraph <- graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
        
        CCs <- components(igraph) # get connected components of igraph
        LCC_size <- max(CCs$csize) # get size of the LCC
        return(LCC_size)
}

random_STRING_graph <- function(x, n_permutations){   # x = number of top genes
        #function that calculates LCC size of 1 random graph:
        
        #do n_permutations iterations of fun() to get a vector of LCC sizes:
        
        LCC_vec <- replicate(n = n_permutations, fun(x), simplify = "array")
        return(LCC_vec)
}

#get vectors with seed genes:

get_seeds <- function(treshold, treshold_value, cell_time){
        if(treshold == "FDR"){
                df <- RML_data_filtered_network_genes %>% filter(Comparison_DE == cell_time, FDR < treshold_value) %>% as.data.frame()
                vec <- c(as.vector(as.character(df$ENTREZ)))
                return(vec)
        } else{
                df <- RML_data_filtered_network_genes %>% filter(Comparison_DE == cell_time, P_value < treshold_value) %>% as.data.frame()
                vec <- c(as.vector(as.character(df$ENTREZ)))
                return(vec)
        }
}

get_top_seeds_based_on_cutoff <- function(cutoff_type, cutoff_values, cell_time){
        
        my_list <- sapply(cutoff_values, function(x, y = cutoff_type, z = cell_time) {
                
                df <- RML_data_filtered_network_genes %>% filter(Comparison_DE == z) %>% 
                        dplyr::slice_min(UQ(rlang::sym(y)), n = x) %>% as.data.frame()
                
                vec <- c(as.vector(as.character(df$ENTREZ)))
                
                return(vec)
                
        }, simplify = FALSE, USE.NAMES = FALSE)
        
        # df <- RML_data_filtered_network_genes %>% filter(Comparison_DE == cell_time) %>% dplyr::slice_min(UQ(rlang::sym(cutoff_type)), n = cutoff_values) %>% as.data.frame()
        # vec <- c(as.vector(as.character(df$ENTREZ)))
        return(unlist(my_list))
        
}


# top_100_seeds <- RML_data_filtered_network_genes %>%  dplyr::filter(!grepl("18Week|S1", Comparison_DE)) %>% dplyr::slice_min(P_value, n = 100)

# top_100_seeds_list <- vector(mode="list", length = 5L)

# # names(top_100_seeds_list) <- unique(top_100_seeds$Comparison_DE)

# for (i in unique(top_100_seeds$Comparison_DE)){
#         top_100_seeds_list[[i]] <- top_100_seeds$ENTREZ[top_100_seeds$Comparison_DE == i]
#         names(top_100_seeds_list[[i]]) <- i
# }

Comparison_DE_vec18 <- c("18Weeks__vGluT2", "18Weeks__Gad2", "18Weeks__Cx43" ,"18Weeks__PV", "18Weeks__SST" )
# treshold_value_vec <- c(rep(0.001, 3), rep(0.001, 2))
treshold_vec <- rep("FDR", 5)
# names_vec <- Comparison_DE_vec18

vector_for_extraction_18_week <- RML_data_filtered_network_genes %>% filter(FDR <= 0.1, grepl("18W", Comparison_DE)) %>% summarize(total = n()) %>% select(total) %>% as.data.frame()
vector_for_extraction_18_week <- as.numeric(vector_for_extraction_18_week$total)
names(vector_for_extraction_18_week) <- RML_data_filtered_network_genes %>% filter(FDR <= 0.1, grepl("18W", Comparison_DE)) %>% summarize(total = n()) %>% select(Comparison_DE) %>% as.data.frame %>% .$Comparison_DE
vector_for_extraction_18_week <- vector_for_extraction_18_week[-4]

vector_for_extraction_18_week[3] <- 1000


seeds_allCells_18_top_FDR_except_PV <- mapply(cell_time = names(vector_for_extraction_18_week), cutoff_type = treshold_vec, cutoff_values = vector_for_extraction_18_week,
                                              FUN = get_top_seeds_based_on_cutoff, USE.NAMES = TRUE, SIMPLIFY = FALSE)

names(seeds_allCells_18_top_FDR_except_PV) <- paste(Comparison_DE_vec18, as.numeric(sapply(seq(100,1000,100), rep, 5), sep = "_size_"))

#funtion that calculates LCC size of genes in ENTREZ_vector: 

LCC_size_from_seed_vec <- function(ENTREZ_vector){
        
        seed_vertices <- ENTREZ_vector
        seed_edges <- STRING_mouse_700_removed_duplicates %>%
                filter(protein1 %in% seed_vertices & protein2 %in% seed_vertices)
        
        seeds_igraph <- graph_from_data_frame(seed_edges, directed = FALSE, vertices = seed_vertices)
        
        seeds_CCs <- components(seeds_igraph)
        seeds_LCC_size <- max(seeds_CCs$csize)
        return(seeds_LCC_size)
        
}

get_components <- function(ENTREZ_vector){
        
        seed_vertices <- unique(ENTREZ_vector)
        seed_edges <- STRING_mouse_700_removed_duplicates %>%
                filter(protein1 %in% seed_vertices & protein2 %in% seed_vertices)
        
        seeds_igraph <- igraph::graph_from_data_frame(seed_edges, directed = FALSE, vertices = seed_vertices)
        
        seeds_CCs <- igraph::components(seeds_igraph)
        #seeds_LCC_size <- max(seeds_CCs$csize)
        return(seeds_CCs)
        
}

# LCCs_all_18 <- sapply(components_18_ig_0.001, function(x){
#         
#         return(names(x$membership)[which(x$membership == which(
#                 components_18_ig_0.001$`18Weeks__vGluT2`$csize == max(components_18_ig_0.001$`18Weeks__vGluT2`$csize)
#                 ))])
#         
# }, USE.NAMES = T, simplify = T)

# saveRDS(LCCs_all_18, "LCCs_all_18.rds" )

allCells_LCC_sizes <- sapply(X = seeds_allCells_18_top_, FUN = LCC_size_from_seed_vec, simplify = FALSE, USE.NAMES = TRUE)

RML_data_26_May_ENTREZ_NAs_filtered

gene_set_size_vec <- unlist(lapply(seeds_allCells_18_top_, length))

gene_set_size_vec <- seq(100,1000,100)

random_graphs_LCCsizes <- function(n_permutations){
        
        #get list of vectors of LCC sizes:
        
        random_LCC_vectors <- mcmapply(x = gene_set_size_vec, n_permutations = n_permutations,
                                       FUN = random_STRING_graph, SIMPLIFY = FALSE, USE.NAMES = TRUE,
                                       mc.preschedule = TRUE, mc.set.seed = TRUE,
                                       mc.silent = FALSE, mc.cores = getOption("mc.cores", 10L),
                                       mc.cleanup = TRUE, affinity.list = NULL)
        
        # random_LCC_vectors <- mapply(x = gene_set_size_vec, n_permutations = n_permutations,
        #                                FUN = random_STRING_graph, SIMPLIFY = FALSE, USE.NAMES = TRUE)
        
        return(random_LCC_vectors)
        
}

start_time <- Sys.time()                                                                                                                                                                                    
aaarandom_LCCs_for_specific_gene_set_sizes_18w <- random_graphs_LCCsizes(n_permutations = 100)                                                                                                                
end_time <- Sys.time()                                                                                                                                                                                     
start_time - end_time 

saveRDS(random_LCCs_for_specific_gene_set_sizes_18w, "random_LCCs_for_specific_gene_set_sizes_18w.rds")

# get the right full randomization vectors for each celltype

# Full_randomization_vecs_allCells <- list()
# 
# for(i in 1:length(seeds_allCells_18_top_)){
#         
#         n <- length(seeds_allCells_18_top_[[i]])
#         vec <- random_LCCs_for_specific_gene_set_sizes_100k[i]
#         Full_randomization_vecs_allCells[i] <- vec
#         names(Full_randomization_vecs_allCells)[i] <- names(seeds_allCells)[i]
#         
# }


######### DEGREE PRESERVING RANDOMIZATION ########################################################################


# make an igraph object with all the nodes from the STRING PPI network:

entire_PPI_network_igraph <- graph_from_data_frame(STRING_mouse_700_removed_duplicates, directed = FALSE, vertices = STRING_mouse_700_vertices)

#get degrees of all vertices in the STRING network

all_vertices_degrees <- degree(entire_PPI_network_igraph, v = V(entire_PPI_network_igraph), loops = FALSE)

# arrange vertex degees into intervals

degree_breaks <- seq(min(all_vertices_degrees), 1240, 10)
degree_interval <- cut(all_vertices_degrees, breaks = degree_breaks, include.lowest = TRUE, right = FALSE)

#df with names of all vertices, their degrees and the coresponding degree interval

degree_df <- data.frame(vertex = names(all_vertices_degrees), degree = all_vertices_degrees, interval = degree_interval)

# a function that samples n_cat random vertices with degrees in the interval_chr interval:

sample_random <- function(interval_chr, n_cat){
        
        filtered_df <- degree_df
        filtered_df$interval <- as.character(filtered_df$interval)
        filtered_df_2 <- filtered_df %>% 
                dplyr::filter(interval == interval_chr) %>%
                as.data.frame()
        
        filtered_df_vertices <- as.vector(as.character(filtered_df_2$vertex))
        random_vertices <- sample(filtered_df_vertices, size = n_cat, replace = FALSE)
        return(random_vertices)
}


# function that samples vertices that have similar degrees as vertices in ENTREZ_vector and measures their LCC size

sample_similar_vertices <- function(ENTREZ_vector){
        x <- degree_df %>% filter(vertex %in% ENTREZ_vector)
        x$interval <- factor(x$interval, levels = unique(x$interval))
        y <- as.vector(table(x$interval))
        names(y) <- as.character(levels(x$interval))
        
        interval_categories_df <- data.frame(interval = names(y), n_genes = y)
        intervals_vec <- as.vector(interval_categories_df$interval)
        intervals_vec_genes <- as.vector(interval_categories_df$n_genes)
        
        
        random_vertices_vec <- mapply(interval_chr = intervals_vec, n_cat = intervals_vec_genes, FUN = sample_random, SIMPLIFY = FALSE, USE.NAMES = TRUE)
        random_vertices_vec <- unlist(random_vertices_vec)
        
        
        random_edges <- STRING_mouse_700_removed_duplicates %>% filter(protein1 %in% random_vertices_vec & protein2 %in% random_vertices_vec)
        
        random_igraph <- graph_from_data_frame(random_edges, directed = FALSE, vertices = random_vertices_vec)
        CCs <- components(random_igraph) 
        LCC_size <- max(CCs$csize) 
        return(LCC_size)
        
}


# function that does n_permutations of the sample_similar_vertices() and returns a vector of random LCC sizes

LCC_sizes_similarDegree_graphs <- function(ENTREZ_vector, n_permutations){
        LCC_vec <- replicate(n = n_permutations, sample_similar_vertices(ENTREZ_vector), simplify = "array")
        return(LCC_vec)
}

start_time <- Sys.time()

aaaallCells_similarDegree_graphs_18w <- mcmapply(ENTREZ_vector = seeds_allCells_18_top_, n_permutation = 100, 
                                                 FUN = LCC_sizes_similarDegree_graphs,  SIMPLIFY = FALSE, USE.NAMES = TRUE,
                                                 mc.preschedule = TRUE, mc.set.seed = TRUE,
                                                 mc.silent = FALSE, mc.cores = getOption("mc.cores", 10L),
                                                 mc.cleanup = TRUE, affinity.list = NULL)
end_time <- Sys.time()

start_time - end_time

saveRDS(allCells_similarDegree_graphs_18w, "allCells_similarDegree_graphs_18w.rds")
