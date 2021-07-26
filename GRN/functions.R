##----------------------------ggplot style---------------------------------------------
theme_min <- function(base_size = 11, base_family = "",
                      border_colour = "black",
                      border_size = 1) {
  
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = border_colour, size = border_size),
      axis.ticks = element_line(colour = border_colour),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black", size = rel(1.2)),
      strip.text.y = element_text(colour = "black", size = rel(1.2)),
      title = element_text(size = rel(0.9)),
      axis.text = element_text(colour = "black", size = rel(1.2)),
      axis.title = element_text(colour = "black", size = rel(1.5)),
      legend.title = element_text(colour = "black", size = rel(1.2)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}




#----------------------------TF information table---------------------------------------------
#function to add a column to data-table containing the HTML code necessary to display the motif logo
addMotifPic <- function(subset_data){ #need to comment this and test to see if it works when I have wifi
  subset_data <- mutate(subset_data, motif_logo = bestMotif)
 
   for (i in 1:nrow(subset_data)){
    motifID <- subset_data$bestMotif[i]
    motifHTML <- glue('<img src=\"http://motifcollections.aertslab.org/v9/logos/{motifID}.png\" height=\"100\" ></img>')
    subset_data$motif_logo[i] = motifHTML
  }

  return(subset_data)
}

#------------------------------Bubble plot--------------------------------------
#' Prepare input for bubble_plot
#' 
#' Load gene expression data from feather, tidy & optionally scale expression,
#' and return the dataframe required as input for the bubble_plot() function
#' 
#' @param gene Character vector, one or more genes of interest to plot
#' @param scale Logical, whether or not to linearly scale gene expression across
#' clusters to [0,1] to improve visualization. Default: TRUE
#' @param show_mean Logical, whether or not to display the mean expression of
#' given genes in a new bubble plot line. Default: FALSE
#' TODO: Use this to provide an option to download the underlying data.
#' 
#' @examples 
#' bubble_prep("Dlx1")
#' bubble_prep <- function(tf,
#'                         scale = TRUE,
#'                         show_mean = FALSE,
#'                         region) {
#'   
#'   # Load the mean expression of genes across clusters, given gene of interest
#'   file_path <- glue("data/joint_{region}/joint_cortex.regulon_activity_per_cluster.feather")
#'   exp <- read_feather(path = file_path,
#'                       columns = c("Cluster", tf)) #%>% 
#'     #dendogram_order from hydra file, commented out for now because no dendogram
#'     #filter(Cluster %in% dendrogram_order)
#'   
#'   # Scale expression of each gene linearly across clusters to [0, 1]
#'   if (scale) {
#'     
#'     exp <- exp %>%
#'       as.data.frame() %>%
#'       tibble::column_to_rownames(var = "Cluster") %>%
#'       apply(2, scales::rescale, to = c(0, 1)) %>%
#'       as.data.frame %>%
#'       tibble::rownames_to_column(var = "Cluster")
#'     
#'   }
#'   
#'   # Convert to long / tidy format with columns: Cluster, Gene, Expression
#'   exp <- exp %>%
#'     gather(., "TF", "Activity", 2:ncol(.))
#'   
#'   # Probably don't need this in the TF plot
#'   # # Load the proportion of cells in each cluster in which each gene was detected,
#'   # # and convert to long / tidy format with columns: Cluster, Gene, Pct1
#'   # pct1 <- read_feather("data/joint_mouse/pct1_per_ID_20190715_cluster.feather",
#'   #                      columns = c("Cluster", gene)) %>%
#'   #   gather(., "Gene", "Pct1", 2:ncol(.))
#'   
#'   # # Join with cluster metadata
#'   # df <- left_join(exp, pct1, by = c("Cluster", "Gene"))  %>% 
#'   #   left_join(metadata, by = c("Cluster" = "Cluster_nounderscore"))
#'   
#'   # Tidy data for plotting
#'   df <- df %>%
#'     
#'     # Order genes to match order input by user
#'     mutate(TF = factor(TF, levels = rev(TF))) %>% 
#'     arrange(TF) %>% 
#'     
#'     # Pad gene names so that the plot takes up a more standardized
#'     # width; to roughly the the # of characters in the gene w/ longest name
#'     # However, letters take up more pixels than spaces, so do less padding
#'     # for genes with longer names
#'     # TODO: Test the (commented) third line inside mutate() and adjust padding as required
#'     mutate(TF_padded = case_when(
#'       str_length(TF) <= 5 ~ str_pad(TF, 15, side = 'right', pad = " "),
#'       between(str_length(TF), 5, 8) ~ str_pad(TF, 12, side = 'right', pad = " ")
#'       #, str_length(Gene) > 8 ~ str_pad(Gene, 9, side = 'right', pad = " ")
#'     )
#'     ) %>% 
#'     mutate(TF_padded = factor(TF_padded, levels = unique(.$TF_padded))) %>% 
#'     
#'     # Order the clusters on the x-axis to match the dendrogram image
#'     mutate(Cluster = factor(Cluster, levels = dendrogram_order)) %>%
#'     
#'     filter(!is.na(Cluster)) %>% 
#'     
#'     # Convert NAs (undetected genes) to 0s -- this ensures all
#'     # clusters have a value for all genes, so that all clusters are plot,
#'     # even if the gene was undetected
#'     replace_na(list(Expression = 0, Pct1 = 0)) %>% 
#'     
#'     # Keep columns
#'     select(Gene, Cluster, Sample, Cell_type, Cell_class, N_cells, Expression, Pct1, Sample, Colour, Gene_padded)
#'   
#'   # Create & append set of rows containing mean expression over all selected genes
#'   if(show_mean) {
#'     
#'     # Create mean expression rows, preserving information for tooltip
#'     mean_exp <- df %>% 
#'       group_by(Cluster, Sample, Cell_type, N_cells, Cell_class, Colour) %>%
#'       summarize(#Gene = "MEAN", 
#'         #Cluster = Cluster,
#'         #Sample = Sample,
#'         #Cell_type = Cell_type,
#'         #Cell_class = Cell_class,
#'         #N_cells = N_cells,
#'         Expression = mean(Expression) 
#'         #Pct1 = mean(Pct1),
#'         #Colour = Colour,
#'         #Gene_padded = "MEAN"
#'       ) %>% 
#'       # Remove the Pct1 value from the mean expression
#'       # and label the mean expression
#'       mutate(Pct1 = 1, Gene = "MEAN", Gene_padded = "MEAN") 
#'     
#'     # Add the rows containing mean expression to the original dataframe,
#'     # removing duplicate rows and ordering them once more by user input,
#'     # except the mean which is placed at the bottom
#'     gene_order_padded <- levels(df$Gene_padded)
#'     df <- bind_rows(df, mean_exp) %>% 
#'       distinct(.) %>%
#'       mutate(Gene_padded = factor(Gene_padded, levels = c("MEAN", gene_order_padded)))
#'     
#'   }
#'   
#'   return(df)
#'   
#' }
#' 
#' 
#' #' Bubbleplot of gene expression
#' #' 
#' #' Generate a bubble plot for genes of interest across clusters in the mouse
#' #' dendrogram, where bubble colour encodes the mean expression in each cluster
#' #' and bubble size encodes the proportion of cells where each gene is detected
#' #'
#' #' @param df Data frame as returned by bubble_prep(), with require columns Cluster,
#' #' Gene_padded, Pct1, and Expression
#' #'
#' #' @return ggplot2 object
#' #'
#' #' @examples
#' #' bubble_prep("Dlx1") %>% bubbleplot()
#' #' 
#' #' @export
#' bubble_plot <- function(df, max_point_size) {
#'   
#'   # Generate plot
#'   p1 <- df %>% 
#'     ggplot(aes(x = Cluster, y = Gene_padded)) +
#'     geom_point(aes(size = Pct1, colour = Expression), alpha = 0.8) +
#'     scale_size_area(max_size = max_point_size) +
#'     scale_color_gradientn(colours = tail(rdbu, 70)) +
#'     theme_min() +
#'     ylab(NULL) +
#'     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
#'                                      colour = joint_mouse_palette, size = rel(0.7)),
#'           panel.grid.major.x = element_line(colour = "grey90"),
#'           panel.border = element_blank(),
#'           axis.ticks.x = element_blank(),
#'           axis.ticks.y = element_blank(),
#'           # Do not show the legend because it is included in the static
#'           # dendrogram image displayed above the bubbleplot
#'           legend.position = "bottom") +
#'     # Put gene labels on the right hand side to improve alignment
#'     scale_y_discrete(position = "right")
#'   
#'   return(p1)
#'   
#' }
#----------------------------ggNet visualisation---------------------------------------------
#function to create an igraph object
#' @param tf the user input vector of transcription factors
#' @param tf_target_gene_info a data frame containing information for each association 
#' between a TF and a target gene; created in data_prep.R, loaded in global.R and used depending 
#' on the region input in input_new()
make_network <- function(tf, tf_target_gene_info, gene_list){
  #add a step to select only the transcription factors that are in the list 
  #create edgelist
  edges <- tf_target_gene_info %>% select(TF, gene, nMotifs, Genie3Weight.weight) %>%
    #and filter it to only the transcription factors that are the input
    filter(TF %in% tf)
  #print(edges)
  #create node list
  #the nodes in this case are all the TFs from user input and all the genes that are regulated by the 
  #transcription factors 
  unique_gene_targets <- unique(edges$gene) 
  unique_TF <- unique(edges$TF)
  nodes <- as_tibble(append(unique_TF, unique_gene_targets, after = length(unique_TF)))
  nodes <- unique(nodes)
  
  #check if the input gene is in the network
  gene_list_in_network <- gene_list[gene_list %in% unique_gene_targets]
  
  #making a basic igraph object with an attribute for each gene indicating if it is a gene 
  #target or transcription factor
  net <- graph_from_data_frame(d=edges, vertices = nodes) %>% 
    set_vertex_attr("Gene_Type", index = unique_gene_targets, "Target Genes") %>%
    set_vertex_attr("Gene_Type", index = gene_list_in_network, "Input Target Genes") %>%
    set_vertex_attr("Gene_Type", index = unique_TF, "TF") 
  
}  
plot_network <- function(net, labelNodes, tf){ 
  #print(net)
  if(labelNodes){
    set.seed(0.926)
    ggnet2(net, color = "Gene_Type", alpha = "Gene_Type", size = "Gene_Type", shape = "Gene_Type", 
           label = labelNodes, label.size = 3, edge.alpha = 0.3,
           mode = "fruchtermanreingold", layout.par = list(niter = 1000),
           size.palette = c("Target Genes" = 2, "TF" = 6, "Input Target Genes" = 6),
           alpha.palette = c("Target Genes" = "1", "TF" = "1", "Input Target Genes" = "1"),
           shape.palette = c("Target Genes" = "19", "TF" = "19", "Input Target Genes" = "19"),
           palette = c("Target Genes" = "grey", "TF" = "lightblue", "Input Target Genes" = "orange"),
           group = 1, text = names(V(net))) +
      guides(size = FALSE, alpha = FALSE, shape = FALSE, color = guide_legend(title = "Gene Type"))
  }
  else{
    set.seed(0.926)
    ggnet2(net, color = "Gene_Type", alpha = "Gene_Type", size = "Gene_Type", shape = "Gene_Type", 
           label = tf, label.size = 3, edge.alpha = 0.3,
           mode = "fruchtermanreingold", layout.par = list(niter = 1000),
           size.palette = c("Target Genes" = 2, "TF" = 6, "Input Target Genes" = 6),
           alpha.palette = c("Target Genes" = "1", "TF" = "1", "Input Target Genes" = "1"),
           shape.palette = c("Target Genes" = "19", "TF" = "19", "Input Target Genes" = "19"),
           palette = c("Target Genes" = "grey", "TF" = "lightblue", "Input Target Genes" = "orange"),
           group = 1, text = names(V(net))) +
      guides(size = FALSE, alpha = FALSE, shape = FALSE, color = guide_legend(title = "Gene Type")) 
  }
  
}

# --------------------------OBSOLETE-cytoscape network visualization----------------------------------------------
# function to create network


#' Create rcytoscape network data
#' 
#' Takes a vector input that contains user selected TFs and output a list of nodeData and edgeData
#' which will be used for createCytoscapeJsNetwork
#' A good to visualize correlations among multiple TFs
#' 
#' @param tf one single tf name character
#' @param TF_target_gene TF_target_gene data, specific for cortex/pon
#' @param unique_TF unique_TF data, specific for cortex/pon
#'
#' @return a list of nodeData and edgeData that are required for generating a rcytoscapejs network object
#' 
#' @examples 
#' TF <- c("Arx","Lef1")
#' # Note that TF_target_gene and unique_TF will be saved in data_cortex list, by data_prep.R
#' nodeData <- create_network(TF, TF_target_gene_pon, unique_TF)$nodes
#' edgeData <- create_network(TF, TF_target_gene_pon, unique_TF)$edges
#' network <- createCytoscapeJsNetwork(nodeData, edgeData)
#' rcytoscapejs2(network$nodes, network$edges)
#' 
# create_network <- function(tf, TF_target_gene, unique_TF, pathway_genes = c(),
#                            shrink_gray = FALSE){ 
#   TF_interest <- filter(TF_target_gene, TF %in% tf)[["TF"]]
#   gene_target <- filter(TF_target_gene, TF %in% tf)[["gene"]]
#   
#   source <- TF_interest
#   target <- gene_target
#   
#   id <- c(TF_interest, gene_target)
#   name <- id
#   nodeData <- data.frame(id,name, stringsAsFactors = FALSE)
#   edgeData <- data.frame(source, target, stringsAsFactors = FALSE)
#   
#   #unique_TF <- unique(TF_target_gene[["TF"]])
#   
#   mutual_target <- edgeData %>% 
#     # a character vector that indicates the nodes that are target of multiple selected TFs
#     count(target) %>%
#     filter(n > 1 & !target %in% tf ) %>%
#     .[[1]]
#   
#   nodeData <- nodeData %>%
#     # you can customize the color using the case_when structure easily,
#     # check the tfs in id column that exist in your vector, then you can control its size,
#     # shape and color easily
#     mutate(color = case_when(id %in% tf ~ "#9d4097", # orange
#                              # orange nodes are tfs that are active in this region
#                              id %in% pathway_genes ~ "green",
#                              id %in% unique_TF ~ "#D6604D", 
#                              id %in% mutual_target ~ "#4fafc6",
#                              TRUE ~ "lightgrey")) %>%
#     mutate(height = case_when(id %in% tf ~ "100",
#                            TRUE ~ "70")) %>%
#     mutate(width = case_when(id %in% tf ~ "100",
#                             TRUE ~ "70"))
#   
#   if(shrink_gray){
#     nodeData <- nodeData %>%
#       mutate(height = case_when(color %in% "lightgrey" ~ "40",
#                                 TRUE ~ "70")) %>%
#       mutate(width = case_when(color %in% "lightgrey" ~ "40",
#                                TRUE ~ "70"))
#     
#   }
#   
#   return(list(nodes = nodeData,
#               edges = edgeData
#   ))
# }
# --------------------------------Helper functions----------------------------------------------------
#' Identify transcription factor data type
#' 
#' Generate a tibble that has two columns indicating whether the tf has ext type,
#' the ext column labels whether that data is ext type
#'
#' @param TF_name_activity_tibble a dataframe/tibble that saves all names of transcription 
#' factor with suffix ext and weights
#' E2f1_extended (133g) and E2f1 (122g) are two examples of tf data detected
#' tf with no '_extended' attached on are those with high confidence annotation
#' tf with extended are those with relatively low confidence annotation
#' The point is to use the high annotation data if we have it, if not, we use
#' the ext type tf data
#' 
#' @return a tibble/dataframe with columns 'type' and 'ext'
#' type shows the tf name without any suffix -- ex. "E2f1"
#' ext shows whether that data is ext/regular
#' regular tf has its own name while extended data has 'ext' in that column
#' This feature is important in following functions, has_regular, tf_regular ...
#' to help retrieve the correct TF name to be used in other dataframes(tf activity in cell..)
#' 
#' @examples 
#' TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))
#' TF_and_ext <- identify_tf(TF_active)
#' 
identify_tf <- function(TF_name_activity_tibble){
  TF_name_activity_tibble %>% 
    rename(name = value) %>%
    mutate(no_space = str_replace_all(name, " \\(.+\\)$", ""))%>% #remove space and the weight
    mutate(TF_type = str_replace(no_space, "_extended","_ext")) %>%
    mutate(ext = str_replace(TF_type, ".+ext", "ext")) %>%
    mutate(type = str_replace(TF_type, "_ext", "")) %>%
    select(name, type, ext)
}

#' has_regular, has_ext
#' 
#' Boolean functions that determine the identity of TF data in the activity datasets
#'
#' @param TF character vector, containing one or more TF names
#' @param TF_and_ext The dataframe/tibble generated using identify_tf, that has three cols:
#' name, type, ext
#'
#' @return A boolean that checks if the datasets has the regular/ext data of that tf input
#'
#' @examples
#' 
#' has_regular("Arx", TF_and_ext) # False
#' has_ext("Arx", TF_and_ext) # True
#' 
#' 
has_regular <- function(TF, TF_and_ext){
  has_regular <- filter(TF_and_ext, type==TF & ext==TF)
  nrow(has_regular)!=0 #boolean
}
has_ext <- function(TF, TF_and_ext){
  is_ext <- filter(TF_and_ext, type==TF & ext=="ext")
  nrow(is_ext)!=0
}

#' tf_exist
#' 
#' This function also uses TF_and_ext, loaded in data_prep.R, can also generate using identify_tf
#' 
#' @param TF character vector, containing one or more TF names
#' @param TF_and_ext TF_and_ext data, specific for cortex/pons
#'
#' @return if the dataset doesn't have data for a tf input, it returns the tf input name; if it has them
#' all, this returns TRUE
#'
#' @examples 
#' TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))
#' TF_and_ext <- identify_tf(TF_active)
#' TF <- "Arx"
#' tf_exist(TF)
#' 
#' 
tf_exist <- function(TF, TF_and_ext){
  for(tf in TF){
    if(has_regular(tf, TF_and_ext) || has_ext(tf, TF_and_ext)){}
    else{return (tf)} # or regurn FALSE
  }
  return (TRUE)
}

# read the corresponding data by tf's identity

#' tf_regular, tf_ext
#'
#' @param TF character vector, containing one or more TF names
#' @param TF_and_ext The dataframe/tibble generated using identify_tf, that has three cols:
#' name, type, ext
#'
#' @return The best represented tf with suffix, still using the same logic: 
#' Use the high annotation data if we have it, if not, we use the ext type tf data
#' 
#'
#' @examples
#' TF_active <- as_tibble(read_rds("data/joint_cortex/joint_cortex.active_regulons.Rds"))
#' TF_and_ext <- identify_tf(TF_active)
#' TF <- c("Arx","Lef1")
#' tf_regular(TF, TF_and_ext)
#' 
tf_regular <- function(TF, TF_and_ext){
  filter(TF_and_ext, type==TF & ext==TF)[[1,1]]
}
tf_ext <- function(TF, TF_and_ext){
  filter(TF_and_ext, type==TF & ext=="ext")[[1,1]]
}

check_tf_input <- function(TF, TF_ref){
  TF_in_data <- TF[TF %in% TF_ref]
  TF_not_data <- TF[!(TF %in% TF_ref)]
  l <- list(
    "TF_in_data" = TF_in_data,
    "TF_not_data" = TF_not_data
  )
}

transform_tf_input <- function(tf, tf_and_ext){
  tf_to_read <- character(0)
  for(TF in tf){ # tf is input tf list, could contain many tfs
    
    if(has_regular(TF, tf_and_ext)){
      tf_to_read[TF] <- tf_regular(TF, tf_and_ext)
    } # a helper to read the corresponding data
    
    else{
      tf_to_read[TF] <- tf_ext(TF, tf_and_ext)
    }
  }
  names(tf_to_read) <- NULL
  return(tf_to_read)
}

# --------------------------------Create data for plots--------------------------------
# NOTE: TF_and_ext is a dataframe (loaded already) that created in order to identify 
# whether the TF data is a regular TF (with high confidence annotation) 
# or ext type(with lower confidence)
# Dogma: if we have regular TF type, we use that data; if we only have ext data, use ext
# if we have no data related to this tf, we will either give an error message or do nothing


#' create Cell/Cluster activity data
#' 
#' Make activity data used in tab2 by either Cell or Cluster, the method would be provided by
#' user's input in Shiny app
#' This function uses feather file that will be read by a certain col to maximize speed,
#' so we switch the paths of the feather file for different brain region
#' 
#' @param tf character vector, containing one or more TF names
#' @param method either by Cell --> use cell data, or by cluster --> use cluster data, 
#' this should be a string indicating the column name
#' @param TF_and_ext TF_and_ext data, specific for cortex/pons
#' @return a dataframe that has a column containing all the cell names and columns of the input tfs
#' the corresponding activity
#' data value (NES) 
#'
#' @examples
#' create_activity_data("Arx", "Cell", "cortex", TF_and_ext)
#' create_activity_data("Pax6", "Cluster", "pons", TF_and_ext_pon)


create_activity_data <- function(tf, method, region, TF_and_ext,
                                 timepoint = NULL, per_sample = FALSE,
                                 bad_cells = ""){ 
  # use the feature of feather data to read certain col to optimize speed
  #if(tf_exist(tf, TF_and_ext) != TRUE){return("TF does not exist")}
  tp <- timepoint
  
  #building path of the feather object to read if the per sample toggle is on
  #method should be joint if per_sample is true
  #print(per_sample)
  #print("ahhh")
  if(per_sample == TRUE){
  
    reg <- switch(region, "cortex" = "ct", "pons" = "po")
  
    #time <- substring(timepoint, 3) commented out for now
    meth <- str_to_lower(method) #per_sample should only be true for method: "joint" or "Cell"
                                 #if its "Cell", turn to lower case and use for DR plots
                                 #if its "joint", then its for per-sample heatmap and load feather
                                   # with per_cluster
    
    if(identical(method, "joint")){
      meth <- "cluster"
    }
    
    path <- glue("data/{reg}_{timepoint}/{reg}_{timepoint}.regulon_activity_per_{meth}.feather")
  }
  else{
    if(!region %in% c("cortex", "pons")) return("Wrong usage: region should be either cortex/pons")
    if(method == "joint"){
      method2 <- "joint_cluster"
    }
    else{
      method2 <- str_to_lower(method)
    }
 
    # set up the path of the feather file to read dependingo on region and cluster or cell
    path <- glue('data/joint_{region}/joint_{region}.regulon_activity_per_{method2}.feather')
  }
  # case-insensitive checking and reading in the first column which corresponds to cell or cluster label
  #why even do this if method2 object has all lowercase already?
  if(str_detect(method,"(?i)Cell")){
    cell_col <- read_feather(path, "Cell")
  }
  else if(str_detect(method,"(?i)Cluster")){
    cell_col <- read_feather(path,"Cluster")
  }
  else if(str_detect(method,"(?i)joint")){
    if(per_sample){
      ID <- "ID_20201028"
      cell_col <- read_feather(path, ID)
      colnames(cell_col) <- "Cluster"
      method <- "Cluster"
    }
    else{
      cell_col <- read_feather(path, "Joint_cluster")
      method <- "Joint_cluster"
    }
  }
  else{return("Wrong usage, method should be Cell/Cluster")} #this should never return
  
  
  
  # add certain tf activity data to the Cell column
  #loops through each factor in tf input and checks to see if the entry has regular or ext forms and extracts
  #data prioritizing the regular factor data and not extended
  #puts each regular or ext TF name in a list 
  activity <- cell_col
  #print(tf)
  for(TF in tf){ # tf is input tf list, could contain many tfs
    tf_to_read <- TF #wouldnt this line combined with the if else if statement double up the data for TF 
                    #with both regular and ext forms 
    if(has_regular(TF, TF_and_ext)){
      tf_to_read <- tf_regular(TF, TF_and_ext) # a helper to read the corresponding data
    }
    else if(has_ext(TF, TF_and_ext)){
      tf_to_read <- tf_ext(TF, TF_and_ext)
    }
    else{
      next # means we don't have that data, we jump over it and do nothing
    }
    #both outcomes are the same here? why is there an if else statement
    #reads the data from file specificed in path using the tf_to_read list 
    #should check what the tf_to_read list actually says 
    
    col <- read_feather(path,tf_to_read)
    
    activity <- add_column(activity,col)
  }
  activity %>%
    select(method, everything()) # move method col to start
  
  if(!(is.null(tp)) & per_sample == FALSE & identical(method, "Cluster")){ #when timepoint has an input and the input is not All
    #split column, select rows for the corresponding timepoints using filter 
    if(identical(tp, "F-All") || identical(tp, "P-All")){
    }
    else{
      activity <- separate(activity, Cluster, into = c("Timepoint", "Cluster"), sep = " ")
      activity <- filter(activity, Timepoint == tp)
      activity <- unite(activity, Cluster, Timepoint, Cluster, sep = " ")
    }
  }
  else if(identical(method, "Cell") & per_sample){ #gets rid of blacklisted cells in the per sample data
    #print(method)
     activity <- activity %>%
       filter(!(Cell %in% bad_cells))
    #print(str(activity))
  }
  activity
}


#--------------------------------Heatmap----------------------------------
# This function takes a colour palette as input,
# and creates the data formats needed to annotate
# a pheatmap with some colours

makePheatmapAnno <- function(palette, column) {
  
  palette <- palette[unique(names(palette))]
  # Make dataframe, retrieve the data frame
  anno_row <- data.frame(cluster = names(palette))
  names(anno_row) <- column
  rownames(anno_row) <- anno_row[[1]]
  
  # Make list containing colours
  side_colors <- list(cluster = palette)
  names(side_colors) <- column
  
  return(list(anno_row = anno_row,
              side_colors = side_colors))
  
}


#' Plot heatmap by cluster/cells
#'
#' @param tf 
#' @param method 
#' @param region 
#' @param TF_and_ext 
#' @param brain_data either forebrain_data or pon_data, eventually will be saved by data_prep.R
#' and loaded at the beginning of app.R as an element in a list
#'
#' @return
#' @export
#'
#' @examples
#' plot_heatmap(c("Arx","Lef1"), "Cluster","cortex", TF_and_ext,forebrain_data)
#' plot_heatmap(c("Arx","Lef1"), "Cell","cortex", TF_and_ext,forebrain_data)
#' plot_heatmap(c("Pax6","Lef1"), "Cluster","pons", TF_and_ext_pon, pon_data)
#' plot_heatmap(c("Pax6","Lef1"), "Cell","pons", TF_and_ext_pon,pon_data)
#' 
plot_heatmap <- function(tf, method, region, TF_and_ext, #brain_data, cell_plot_num = 300, 
                         timepoint = NULL, per_sample = FALSE){
  # sanity checking
  if(!region %in% c("cortex", "pons")) return("Wrong usage: region should be either cortex/pons")
  if(!method %in% c("Cell","Cluster", "joint")) return("Wrong usage, method should be Cell/Cluster")
  
  # if(method == "Cell"){
  #   # 1. create the activity data for plotting 
  #   act_cell <- create_activity_data(tf, "Cell",region, TF_and_ext) %>%
  #     mutate(Cluster = gsub("_"," ",brain_data[["Sample_cluster"]])) %>%
  #     filter(!grepl("BLACKLIST", Cluster)) %>% # filter out bad samples
  #     sample_n(cell_plot_num) %>%  # randomly sample it
  #     tibble::column_to_rownames(var = "Cell") # make that column name as row name ...
  #   
  #   anno_row_cell <- select(act_cell, Cluster)
  #   # change the anno_row, since we change the color palettes
  #   new_anno_row_cell <- anno_row_cell %>%
  #   mutate(Cluster = gsub(pattern = ".* ", replacement = "", Cluster))
  #   rownames(new_anno_row_cell) <- rownames(anno_row_cell) # re-assign the rownames
  #   
  #   act <- select(act_cell, -Cluster) # must remove Cluster data before plotting
  #   
  #   # customized for plotting by cell
  #   anno_col <- new_anno_row_cell # assign to the same variable for plotting
  #   cell_width_plot <- 2
  #   show_colname_plot <- FALSE
  # }
  if(method == "Cluster"){
    act <- create_activity_data(tf, "Cluster",region, TF_and_ext, timepoint)
      #sample_n(cluster_plot_num) %>% # randomly sample it
    #str(act)
    act <- column_to_rownames(act, var = "Cluster") # make that column name as row name ...
    
    # change the anno_row, since we change the color palettes
    new_anno_row <- hm_anno$anno_row %>%
      mutate(Cluster = gsub(pattern = ".* ", replacement = "", Cluster))
    rownames(new_anno_row) <- rownames(hm_anno$anno_row) # re-assign the rownames
    # note that the rownames correspond to the col names of the matrix t(act_cluster)
    # customized for plotting by cluster
    anno_col <- new_anno_row # this is loaded by data_prep.R
    #print(anno_col)
    cell_width_plot <- 20
    cell_height_plot <- 20
    if (identical(timepoint, "F-All") || identical(timepoint, "P-All")){
      cell_width_plot <- 7
      cell_height_plot <- 10
    }
    show_colname_plot <- TRUE
    title <- glue('Transcription Factor Regulon Activity at Developmental Time: {timepoint}')
  }
  else if(method == "joint"){ #plot heat map by joint cluster 
    
    act <- create_activity_data(tf, "joint", region, TF_and_ext, per_sample = per_sample,
                                timepoint = timepoint)
    #print(act)
    if(per_sample == TRUE){
      col_to_row <- "Cluster"
      new_anno_row <- hm_anno$anno_row %>%
        mutate(Cluster = gsub(pattern = ".* ", replacement = "", Cluster))
      rownames(new_anno_row) <- rownames(hm_anno$anno_row)
      
    }
    else{
      col_to_row <- "Joint_cluster"
      new_anno_row <- act %>% mutate(Cluster = Joint_cluster) %>%
        column_to_rownames("Joint_cluster") %>% select(Cluster)
      #print(new_anno_row)
      
      #%>%
       # column_to_rownames("Cluster")
      #row.names(new_anno_row) <- new_anno_row$Cluster
    }
  
    act <- act %>%
      column_to_rownames(var = col_to_row) 
  
    # change the anno_row, since we change the color palettes
    # new_anno_row <- hm_anno$anno_row %>%
    #   mutate(Cluster = gsub(pattern = ".* ", replacement = "", Cluster))
    # rownames(new_anno_row) <- rownames(hm_anno$anno_row) # re-assign the rownames
    # note that the rownames correspond to the col names of the matrix t(act_cluster)
    # customized for plotting by cluster
    anno_col <- new_anno_row # this is loaded by data_prep.R
    #print(anno_col)
    cell_width_plot <- 20
    cell_height_plot <- 20
    show_colname_plot <- TRUE
    title <- "Transcription Factor Regulon Activity per Cluster"
  }
  cluster_row <- FALSE
 #do not do row clustering if there is only one TF selected
  if(length(tf) > 1){
   cluster_row <- TRUE
  }
  pheatmap::pheatmap(t(act),
                     show_colnames = show_colname_plot,
                     scale = "none",
                     border_color = NA,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     main = title,
                     cluster_rows = cluster_row,
                     annotation_col = anno_col,
                     # change the default color annotation
                     annotation_colors = master_palette, # loaded by data_prep.R
                     annotation_legend = FALSE,
                     cellwidth = cell_width_plot,
                     cellheight = cell_height_plot)
} 

#----------------------------Dimension reduction---------------------------------------------
#' Make UMAP clustering scatterplot
#'
#' @param tf_number Either 1 or 2. In the tf input vector we get from user in Shiny app, there could be
#' multiple tfs, but we only support plotting two tfs since the scatterplot is big
#' @param overall_brain_data metadata (forebrain_data or pon_data), saved in data_cortex
#' and data_pons
#' @param cell_activity_data made by create_activity_data() function given the tf input
#' @param sample_number we eliminate half of the cell samples to relieve the burden of 
#' the RAM to speed up plotting, since we have over 37000 cells(samples), we randomly sample 
#' 13000 to optimize speed, but one can also specify this value to see fewer or more sample points
#'
#' @return a UMAP scatter plot that shows in which cluster(region) the tf expresses the most
#'
#' @examples
#' tf <- c("Arx","Lef1")
#' activity_test_tf1 <- create_activity_data(tf, "Cell","cortex", data_cortex$TF_and_ext)
#' plot_UMAP(tf_number = 1,data_cortex$overall, activity_test_tf1)
#' 
plot_UMAP <- function(tf_number, cell_metadata, cell_activity_data, dim_red_type){ #cell_metadata is the tsv with the 
  #embedding coordinates for each cell in the dataset
  # if(tf_number == 1) tf_plot <- 2 # number of col, the first col is Cell, so start from 2
  # else if(tf_number == 2) tf_plot <- 3 
  # else{return(
  #   "Wrong usage, now we only support plotting two tfs since the scatterplot is big"
  # )}
  
  tf_plot <- tf_number + 1 #replaces the above control flow

  
  activity_tf <- cell_activity_data[,tf_plot][[1]] #extracts the TF activity from the cell_activity_data
  #and appends it to the cell_metadata to make cell meta with activity to plot
  
  cell_meta_with_activity <- mutate(cell_metadata, activity_tf = activity_tf)
  
  x_axis <- switch(dim_red_type, "umap" = "UMAP1", "tsne" = "tSNE_1", "pca" = "PC1")
  y_axis <- switch(dim_red_type, "umap" = "UMAP2", "tsne" = "tSNE_2", "pca" = "PC2")
  
  ggplot(data = cell_meta_with_activity, mapping = aes_string(x = x_axis, y = y_axis))+
    geom_point(aes(color = activity_tf), alpha = 0.2)+
    scale_color_gradient(low = "grey", high = "red")+
    theme_min() + labs(color = 'TF Activity')
  
}
#if I want to include labeled clusters, then I need to map cells to the clusters 
#place a label at the mean of the umap coordinates for the cells that belong in that cluster

color_by_cluster <- function(cell_metadata, cluster_palette, dim_red_type, cluster_label, 
                             per_sample = FALSE){
  
  x_axis <- switch(dim_red_type, "umap" = "UMAP1", "tsne" = "tSNE_1", "pca" = "PC1")
  y_axis <- switch(dim_red_type, "umap" = "UMAP2", "tsne" = "tSNE_2", "pca" = "PC2")
  
  # Store the center points (medians) of each cluster and add labels courtesy of Bhavyaa 
  
  var_group <- "Joint_cluster"
  if(per_sample){
    var_group <- "ID_20190715_with_blacklist"
  }
  # centers <- cell_metadata %>%
  #   group_by(get(var_group)) %>%
  #   summarise(center_x = median(get(x_axis)),
  #             center_y = median(get(y_axis)))
  
  #print("step1")
  #print(centers)
  
  gg <- ggplot(data = cell_metadata, mapping = aes_string(x = x_axis,y = y_axis))+
    geom_point(aes(color = get(var_group)), alpha = 0.2) + theme_min() + theme(legend.position="bottom") + 
    guides(fill=guide_legend(nrow=5, byrow=TRUE)) + scale_color_manual(values = cluster_palette) +
    labs(color = 'Cluster Label')
  
 # print("step2")
  if(cluster_label){
    
    centers <- cell_metadata %>%
      group_by(get(var_group)) %>%
      summarise(center_x = median(get(x_axis)),
                center_y = median(get(y_axis)))

    gg <- gg + ggrepel::geom_label_repel(data = centers,
                                         aes(x = center_x, y = center_y),
                                         label = centers$'get(var_group)',
                                         size = 4,
                                         segment.color = 'grey50',
                                         fontface = 'bold',
                                         alpha = 0.8,
                                         segment.alpha = 0.8,
                                         label.size = NA,
                                         force = 2,
                                         segment.size = 0.5,
                                         arrow = arrow(length = unit(0.01, 'npc')))
  }
 #print("step3")
  return(gg)
}

#----------------------------Time Course Ribbon Plot---------------------------------------------
#need to maybe change the colors, select out the numbers and rename legend, just cosmetic things

#' make cell metadata of certain region, cortex/pon
#'
#' @param cell_metadata a dataframe, forebrain_data or pons_data, saved in data_cortex / data_pons
#' @param part  a string, either "cortex" or "pons"
#' @return a dataframe/tibble of metadata with columns Age, Cell, Prefix and Cluster, used for plotting
#' timeseries
#'
#' @examples
#' cell_metadata_cortex <- read_tsv("data/joint_cortex/joint_cortex.metadata.tsv")
#' cell_metadata_cortex <- create_cell_metadata(cell_metadata_cortex)
#' 
create_metadata_timeseries <- function(cell_metadata, part){
  if(part == "cortex") level <- c("Forebrain E12.5",
                                  "Forebrain E15.5",
                                  "Forebrain P0",
                                  "Forebrain P3",
                                  "Forebrain P6")
  else if (part == "pons") level <- c("Hindbrain E12.5",
                                      "Pons E15.5",
                                      "Pons P0",
                                      "Pons P3",
                                      "Pons P6")
  else{(return("Wrong usage, input either cortex or pons"))}
  
  cell_metadata %>% 
    select(Age = Sample, Cell, Cluster = Sample_cluster) %>% 
    # In this case, we remove the "prefix" of the Cluster column, so that we are
    # simply left with the abbreviation representing the cell type, so that 
    # we can link the cells of the same cell type across ages
    separate(Cluster, into = c("Prefix", "Cluster"), sep = "_") %>% 
    mutate(Age = factor(Age, levels = level)) %>% 
    arrange(Cell)
  
}

# create_cell_metadata_pon <- function(metadata_part){
#   metadata_part %>% 
#     select(Age = orig.ident, Cell, Cluster = ID_20190715_with_blacklist_and_refined) %>% 
#     # In this case, we remove the "prefix" of the Cluster column, so that we are
#     # simply left with the abbreviation representing the cell type, so that 
#     # we can link the cells of the same cell type across ages
#     separate(Cluster, into = c("Prefix", "Cluster"), sep = "_") %>% 
#     mutate(Age = factor(Age, levels = c("Hindbrain E12.5",
#                                         "Pons E15.5",
#                                         "Pons P0",
#                                         "Pons P3",
#                                         "Pons P6"))) %>% 
#     arrange(Cell)
#   
# }


#' Translate transcription factor name version
#' 
#' Since we have two datasets with different tf name format, (with ext and weights or not)
#' and those tf names are essential for retrieving data in various datasets,
#' this function use a clean tf and return a tf with ext and weight type
#'
#' @param tf a character vector of a tf without any suffix. EX: "Arx" 
#' @param tf_dataframe a one column dataframe that contains all the TF names with suffix 
#' , get from the rownames of the cell binary activity data for timeseries tab3
#' 
#' @return If the dataframe contains the tf input, it return a best represented tf name
#'  with ext and weight suffix. If not, it returns FALSE
#'
#' @examples
#' tf_df <- data_cortex$ # as_tibble(rownames(activity))
#' translate_tf("Arx",tf_df)  # Arx_extended (21g)
#' translate_tf("Brahl",tf_df) # NULL
#' translate_tf(c("Lef1","Arx"),tf_df) # "Lef1 (22g)"         "Arx_extended (21g)"
translate_tf <- function(tf_list, tf_dataframe){
  tf_info <- identify_tf(tf_dataframe)
  l <- c()
  for(TF in tf_list){
    if(has_regular(TF, tf_info)){
      l <- c(l, tf_regular(TF,tf_info)) # a helper to read the corresponding data
    }
    else if(has_ext(TF, tf_info)){
      l <- c(l, tf_ext(TF,tf_info))
    }
    else{
      next
    }
  }
  if(is.null(l)) return (FALSE) # means we don't have that data at all
  else{return (l)} # return the list
}



#' Plot timeseries
#' @author Selin Jessa and Anthony Ma, most credit to Selin and Anthony puts codes into the function
#' @param TF a character vector that contains one or multiple TF names, that may need to be 
#' transformed by translate_tf function to change its string form
#' @param cell_metadata cell_metadata_cortex, loaded in data_prep.R, a dataframe
#' with Age, Cell, Prefix, Cluster columns, this data is specific to forebrain cells/pon cells
#' @param activity binary_activity data, loaded in data_prep.R
#'
#' @return a plot that displays the percentage level of that tf along with several the time points
#'
#' @examples
#' tf_df <- as_tibble(rownames(activity))
#' TF <- translate_tf("Lef1",tf_df)
#' binary_activity <- data_cortex$binary_activity
#' cell_metadata_cortex <- create_metadata_timeseries(data_cortex$cell_metadata, "cortex")
#' plot_timeseries(TF,cell_metadata_cortex, binary_activity)
#' 
plot_timeseries <- function(TF,cell_metadata, activity, make_plotly = FALSE, show_legend = TRUE){
  
  activity <- activity[TF, ] %>% 
    {data.frame("TF" = .)} %>% 
    tibble::rownames_to_column(var = "Cell") %>% # the original activity vector has names
    arrange(Cell)
  
  if(!all(cell_metadata$Cell == activity$Cell)) return (-1)
  # Add the TF activity to the new dataframe
  ribbon_df <- cell_metadata
  ribbon_df$TF <- activity$TF
  
  ribbon_df <- ribbon_df %>% 
    filter(!grepl("BLACKLIST", Cluster))
  ribbon_df_celltype_frac <- ribbon_df %>% 
    group_by(Age) %>% 
    # Total cells at each age
    mutate(total = n()) %>% 
    group_by(Age, Cluster) %>%
    # Proportion of TF+ cells per cluster, per age
    mutate(frac = sum(TF > 0) / total) %>% 
    distinct(Age, Cluster, frac) %>% 
    ungroup()
  
  ribbon_df_cum_frac <- ribbon_df %>% 
    group_by(Age) %>% 
    summarize(cumfrac = sum(TF > 0) / n()) %>% 
    ungroup()
  
  timepoints2 <- ribbon_df$Age
  clusters <- ribbon_df$Cluster
  
  df = data.frame(cluster = rep(unique(clusters), length(unique(timepoints2))),
                  stage = do.call(c, lapply(as.character(unique(timepoints2)), rep, times = length(unique(clusters)))))
  
  df$ranking = match(df$cluster, names(colours))
  df = df[order(df$stage, df$ranking),]
  
  df <- left_join(df, select(ribbon_df_celltype_frac, cluster = Cluster, stage = Age, frac)) %>% 
    mutate(frac = replace_na(frac, 0)) %>% 
    left_join(select(ribbon_df_cum_frac, stage = Age, cumfrac))
  
  df$xpos = match(df$stage, unique(timepoints2))
  
  #view(df)
  
  plot <- df %>%
    ggplot(aes(x = xpos, y = frac, fill = cluster)) +
    geom_area(stat = "identity", show.legend = show_legend) +
    scale_fill_manual(values = colour_palette, drop = FALSE, name = "") +
    scale_x_continuous(breaks = seq_along(unique(df$stage)),
                       labels = c("E12.5", "E15.5", "P0", "P3", "P6"),
                       limits = c(1, length(unique(df$stage)))) +
    labs(x = "Developmental Age", y = "Proportion", title = TF) +
    guides(fill = guide_legend(ncol = 5)) +
    theme_min() + 
    theme(legend.position = "bottom") 
  
  if(make_plotly) {
    return (ggplotly(plot, tooltip = "cluster") %>% style(hoveron = "points + fills"))
  }
  else{return(plot)}
}
#--------------------------Active_specific------------------------
active_specific_prep <- function(sample, cluster){
  path <- glue("data/{sample}/{sample}.active_specific_tf.Rds")
  data <- readRDS(path) 
  #print()
  tf_table <- data$tf_table[[cluster]] %>% 
    mutate(is_ext = grepl("extended", data$tf_table[[cluster]]$TF))
  FC_df <- data$FC_df
  AUC_df <- data$AUC_df 
  
  send_back <- list("tf_table" = tf_table,
                    "FC_df" = FC_df,
                    "AUC_df" = AUC_df)
  #names(data) <- cluster
}
plot_scatter <- function(data, fc, cluster){

  to_label <- data %>% filter(AUC_FC > fc) %>% mutate(why = gsub("_extended", "+", TF)) %>%
    mutate(TF = gsub("\\(.+\\)", "", why))
  #print(to_label)
  
  ggplot(data = data, mapping = aes(AUC_out, AUC_in)) +
    geom_point(aes(color = AUC_FC, shape = is_ext), size = 4, alpha = 0.6) + 
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(8, "RdBu"))) + 
    scale_shape_manual(values = c(19, 17)) +
    geom_abline(slope = fc, intercept = 0, show.legend = TRUE) +
    labs(color = 'AUC Fold\nChange', shape = 'Extended\nRegulon', 
         y = glue('Average AUC in {cluster}'), x = 'Average AUC in All Other Clusters',
         title = cluster) +
    theme_min() +
    theme(plot.title = element_text(size = 14, face="bold"),
          plot.margin = unit(c(1.75,0,0,0), "cm")) +
    guides(size = FALSE) + 
    ggrepel::geom_label_repel(data = to_label,
                              aes(x = AUC_out, y = AUC_in),
                              label = to_label$TF,
                              size = 4,
                              segment.color = 'grey50',
                              fontface = 'bold',
                              alpha = 0.8,
                              segment.alpha = 0,
                              label.size = NA,
                              force = 0.5,
                              segment.size = 0.5,
                              arrow = NULL)
}
plot_bar_list <- function(data, tf){
  
  data <- data %>% mutate(Cluster = gsub("_", " ", Cluster))
  #print(data)
  palette <- hm_anno$side_colors$Cluster[names(hm_anno$side_colors$Cluster) %in% data$Cluster]
  #print(palette)
  #print(length(data$Cluster))
  #print(length(palette))
  #(data)
  #purrr::map(tf, ~print(.x) )
  # 1. loop over genes, then each index can be referred to with .x
   purrr::map(.x = tf, .f = ~data %>% select(.x, Cluster) %>% 
                       arrange(desc(.x)) %>% head(30) %>%
               ggplot(aes(x = get(.x), y = reorder(Cluster, get(.x)))) +
               geom_bar(aes(colour = Cluster, fill = Cluster), stat = "identity") +
               scale_color_manual(values = palette) +
               scale_fill_manual(values = palette) +
               ggtitle(.x) + theme_min() + theme(legend.position = "none", plot.title = element_text(size = 14, face="bold")) +
               labs(y = "Cluster", x = glue('{.x} AUC'))
              ) %>% 
    # 2. pipe the output (a list) into plot_grid, using {} to indicate it's not the first argument
    {plot_grid(plotlist = .)}
  
}


