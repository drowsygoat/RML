if (!require("pacman")) install.packages("pacman", repos = getOption("repos"))
pacman::p_load(tidyverse, igraph, matrixStats, parallel)

# file.create("log_randomization_11_october.txt")
# sink(file = "log_randomization_11_october.txt")

RML_data_26_May_ENTREZ_NAs_filtered <- readRDS("../RML_data_26_May_ENTREZ_NAs_filtered.rds")
RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE <- gsub(" ", "", RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE)
STRING_mouse_700_removed_duplicates_extended <- readRDS("../STRING_mouse_700_removed_duplicates_extended.rds")

#extract all nodes from STRING PPI network:

STRING_mouse_700_vertices <- STRING_mouse_700_removed_duplicates_extended %>% dplyr::select(protein1, protein2) %>%
        flatten_chr() %>%
        unique()

# filter RML datatable with p values so that all remaining genes will be in the STRING PPI network:
RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE <- gsub(" ", "", RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE)

RML_data_filtered_network_genes <- RML_data_26_May_ENTREZ_NAs_filtered %>%
        dplyr::filter(ENTREZ %in% STRING_mouse_700_vertices) %>%
        dplyr::group_by(Comparison_DE)



Comparison_DE_vec_10W <- c("10Weeks__vGluT2", "10Weeks__Gad2", "10Weeks__Cx43" ,"10Weeks__PV", "10Weeks__SST" )
treshold_value_vec <- c(rep(0.1, 3), rep(0.01, 2))
treshold_vec <- c(rep("FDR", 3), rep("P_value", 2))
names_vec <- gsub("10Weeks__", "", Comparison_DE_vec_10W)
# networks for all cells 
allCells_unweighted_networks <- vector("list", 3)

for(i in 1:length(Comparison_DE_vec_10W)){
        all_ENTREZ <-RML_data_filtered_network_genes %>%
                ungroup() %>%
                filter(Comparison_DE == Comparison_DE_vec_10W[i]) %>% 
                dplyr::select(ENTREZ) %>% 
                flatten_chr() %>%
                unique() 
        
        vertices <- all_ENTREZ[all_ENTREZ %in% STRING_mouse_700_vertices]
        edges <- STRING_mouse_700_removed_duplicates_extended%>% filter(protein1 %in% vertices, protein2 %in% vertices)
        network <- graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
        
        allCells_unweighted_networks[[i]] <- network
        names(allCells_unweighted_networks)[i] <- gsub("10Weeks__", "", Comparison_DE_vec_10W[i])
        
        
}



all_cells_degree_dfs <- vector("list", 3)

for(i in 1:length(Comparison_DE_vec_10W)){
        
        degrees <- degree(allCells_unweighted_networks[[i]], v = V(allCells_unweighted_networks[[i]]), loops= TRUE)
        
        # arrange vertex degrees into intervals
        breaks <- c()
        j <- max(degrees)
        
        while (j >= min(degrees)){
                x <- sort(degrees, decreasing = TRUE)
                k <- j
                if(length(x[x %in% k:j]) < 100){
                        while (length(x[x %in% k:j]) < 100){
                                k <- k-1
                        }
                        breaks[length(breaks)+1] <- k              
                        j <- k-1
                }else{
                        breaks[length(breaks)+1] <- j
                        j <- j -1
                }
        }
        
        breaks <- c(sort(breaks), max(degrees))
        
        degree_interval <- cut(degrees, breaks = breaks, include.lowest = TRUE, right = FALSE)
        degree_df <- data.frame(vertex = names(degrees), degree = degrees, interval = degree_interval)
        
        
        all_cells_degree_dfs[[i]] <- degree_df
        names(all_cells_degree_dfs)[i] <- gsub("10Weeks__", "", Comparison_DE_vec_10W[i])
        
}



DP_LCC_size_randomization <- function(degree_df, n_randomizations, ENTREZ_vector, PPI){
        
        
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
        
        
        similar_vertices_LCCs <- function(ENTREZ_vector){
                x <- degree_df %>% filter(vertex %in% ENTREZ_vector)
                x$interval <- factor(x$interval, levels = unique(x$interval))
                y <- as.vector(table(x$interval))
                names(y) <- as.character(levels(x$interval))
                
                
                random_vertices_vec <- mapply(interval_chr = names(y), n_cat = y, FUN = sample_random, SIMPLIFY = FALSE, USE.NAMES = TRUE)
                random_vertices_vec <- unlist(random_vertices_vec)
                
                random_subgraph <- induced_subgraph(PPI, vids = random_vertices_vec, impl = "auto")
                
                C_random_subgraph <- components(random_subgraph)
                LCC_size <- max(C_random_subgraph$csize)
                
                return(LCC_size)
                
                
        }
  
        
        LCC_vec <- replicate(n = n_randomizations, similar_vertices_LCCs(ENTREZ_vector), simplify = TRUE)
        
        
        return(LCC_vec)
        
}



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

seeds_allCells<- mapply(cell_time = Comparison_DE_vec_10W, treshold = treshold_vec, treshold_value = treshold_value_vec,
                        FUN = get_seeds, USE.NAMES = TRUE, SIMPLIFY = FALSE)

names(seeds_allCells) <- names_vec


start_time <- Sys.time()

random_LCCs_all_Cells_1000x <- mapply(degree_df = all_cells_degree_dfs, 
				      ENTREZ_vector = seeds_allCells, 
                                      PPI = allCells_unweighted_networks, n_randomizations = 1000,
                                      FUN = DP_LCC_size_randomization, 
                                      SIMPLIFY = FALSE, USE.NAMES = TRUE)

end_time <- Sys.time()
start_time - end_time

saveRDS(random_LCCs_all_Cells_1000x, paste0(round(runif(1)*10000, 0),"_random_LCCs_all_Cells_1000x", ".rds"))


