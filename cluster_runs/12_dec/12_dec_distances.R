args = commandArgs(trailingOnly=TRUE)

if (!require("pacman")) install.packages("pacman", repos = getOption("repos"))
pacman::p_load(tidyverse, igraph, parallel, matrixStats)

seeds_10_weeks <- readRDS("seeds10weeks.rds")
STRING_mouse_700_removed_duplicates_extended <- readRDS("STRING_mouse_700_removed_duplicates_extended.rds")
RML_data_26_May_ENTREZ_NAs_filtered <- readRDS("RML_data_26_May_ENTREZ_NAs_filtered.rds")


#extract all nodes from STRING PPI network:

STRING_mouse_700_vertices <- STRING_mouse_700_removed_duplicates_extended %>% 
        dplyr::select(protein1, protein2) %>%
        flatten_chr() %>%
        unique()

# filter RML datatable with p values so that all remaining genes will be in the STRING PPI network:

RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE <- gsub(" ", "", RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE)

RML_data_filtered_network_genes <- RML_data_26_May_ENTREZ_NAs_filtered %>%
  dplyr::filter(ENTREZ %in% STRING_mouse_700_vertices) %>%
  dplyr::group_by(Comparison_DE)

# entire STRING PPI network igraph object:

entire_PPI_network_igraph <- graph_from_data_frame(STRING_mouse_700_removed_duplicates_extended, directed = FALSE, vertices = STRING_mouse_700_vertices)

#get degrees of all vertices in the STRING network

all_vertices_degrees <- degree(entire_PPI_network_igraph,
                               v = V(entire_PPI_network_igraph),
                               loops = TRUE)

# arrange vertex degrees into intervals

degree_breaks <- seq(min(all_vertices_degrees), 1240, 10)

degree_interval <- cut(all_vertices_degrees, 
                       breaks = degree_breaks, 
                       include.lowest = TRUE, 
                       right = FALSE)

# df with names of all vertices, their degrees and the corresponding degree interval

degree_df <- data.frame(vertex = names(all_vertices_degrees),
                        degree = all_vertices_degrees,
                        interval = degree_interval)

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

# function that samples vertices that have similar degrees as vertices in ENTREZ_vector and gets their ds values

similar_vertices_ds <- function(x, ENTREZ_vector){
        
  x <- degree_df %>% filter(vertex %in% ENTREZ_vector)
  x$interval <- factor(x$interval, levels = unique(x$interval))
  y <- as.vector(table(x$interval))
  names(y) <- as.character(levels(x$interval))
  
  random_vertices_vec <- mapply(interval_chr = names(y), n_cat = y, FUN = sample_random, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  random_vertices_vec <- unlist(random_vertices_vec)
  
  all_distances <- distances(entire_PPI_network_igraph,
                             v = random_vertices_vec,
                             to = random_vertices_vec,
                             weights = NA) # calculate the shortest paths between all pairs of vertices
 
  all_distances_removed_zeros <- replace(all_distances, all_distances == 0, NA) # replace 0 distances with NAs (0 = shortest path length from vertex to itself)
  
  #get all the shortest distances to the nearest seed gene
  
  ds_vec <- rowMins(all_distances_removed_zeros, na.rm = TRUE)
  ds_vec <- ds_vec[which(!is.infinite(ds_vec))]
  
  return(ds_vec)

}

# repeat similar_vertices_ds() n_randomizations times to get a list of vectors with shortest distances

multiple_random_ds_vecs <- function(n_randomizations){   # x = vector with ENTREZ gene IDs

  ds_list <- mclapply(1:n_randomizations,
                      FUN = similar_vertices_ds,
                      ENTREZ_vector = seeds_10_weeks[[args[1]]],
                      mc.preschedule = TRUE, mc.set.seed = TRUE,
                      mc.silent = FALSE, mc.cores = getOption("mc.cores", 10L),
                      mc.cleanup = TRUE, mc.allow.recursive = TRUE,
                      affinity.list = NULL)
  return(ds_list)
  
}
b <- Sys.time()
a <- multiple_random_ds_vecs(n_randomizations = as.numeric(args[2]))
c <- Sys.time()
c-b
saveRDS(a, paste0("random_distances_", args[2], "_", args[1], ".rds"))

