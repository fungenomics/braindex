##----------------------------ggplot style---------------------------------------------
#theme for plots, adopted from Selin in the Clusters App
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
#function to add a column to data-table containing the HTML necessary to display the binding motif logo
#for use in the Table tab
#' @param subset_data dataframe containing user selected TFs and its predicted regulatory relationships
addMotifPic <- function(subset_data){
  subset_data <- mutate(subset_data, motif_logo = bestMotif)
 
   for (i in 1:nrow(subset_data)){
    motifID <- subset_data$bestMotif[i]
    motifHTML <- glue('<img src=\"http://motifcollections.aertslab.org/v9/logos/{motifID}.png\" height=\"100\" ></img>')
    subset_data$motif_logo[i] = motifHTML
  }

  return(subset_data)
}


#----------------------------ggNet visualisation---------------------------------------------
#function to create an igraph object
#' @param tf the user input vector of transcription factors
#' @param tf_target_gene_info a data frame containing information for each association 
#' between a TF and a target gene; created in data_prep.R, loaded in global.R and used depending 
#' on the region input in input_new()
make_network <- function(gene_input, tf_target_gene_info, gene_list, network_by_target_gene = FALSE){

  #create edgelist
  edges <- tf_target_gene_info %>% select(TF, gene, nMotifs, starts_with("Genie3Weight")) #%>%
  #  filter(TF %in% tf)
  
  #and filter it to only the transcription factors or genes that are the input
  
  if(network_by_target_gene){edges <- edges %>% filter(gene %in% gene_input)} #at this point, this is a matrix relating genes nad 
  else{edges <- edges %>% filter(TF %in% gene_input)} #TFs 
 
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

#plots igraph network using ggNet
#' @param net igraph object made with make_network function
#' @param labelNodes User input that allows all nodes to be labelled; if false, only TF nodes are labelled
#' @param tf list of TFs input by user; used to label just the TF nodes if labelNodes is false

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


#'@param TF a user selected TF list
#'@param TF_ref a list of TFs that are in the dataset
#'
#'subsets TF based on if each element appears in TF_ref
check_tf_input <- function(TF, TF_ref){
  TF_in_data <- TF[TF %in% TF_ref]
  TF_not_data <- TF[!(TF %in% TF_ref)]
  l <- list(
    "TF_in_data" = TF_in_data,
    "TF_not_data" = TF_not_data
  )
}

#'@param tf a user selected TF list, gene symbols
#'@param TF_and_ext a list of TFs in the dataset and whether the regulon is a regular or extended regulon
#'
#'Converts gene symbols in tf to the regulon name that is used in the TF activity feather files  
transform_tf_input <- function(tf, tf_and_ext){
  tf_to_read <- character(0)
  for(TF in tf){ # tf is input tf list, could contain many tfs
    
    if(has_regular(TF, tf_and_ext)){
      tf_to_read[TF] <- tf_regular(TF, tf_and_ext) #if the TF has a non-extended/regular regulon, then the regular regulong name is used
    } 
    
    else{
      tf_to_read[TF] <- tf_ext(TF, tf_and_ext) #gets the extended regulon name for the TF
    }
  }
  names(tf_to_read) <- NULL
  return(tf_to_read)
}


#' Determine if a background colour is dark enough to warrant white text
#' 
#' @param hex_color String, colour in hex colour format e.g. #000000
#' 
#' @return TRUE if the colour is dark enough (arbitrary)
dark <- function(hex_color) {
  
  red <- substr(hex_color, 2, 2)
  green <- substr(hex_color, 4, 4)
  blue <- substr(hex_color, 6, 6)
  dark_nums <- c(0:8)
  
  if ((red %in% dark_nums && blue %in% dark_nums) || 
      (red %in% dark_nums && green %in% dark_nums) ||
      (green %in% dark_nums && blue %in% dark_nums)) {
    
    return(TRUE)
    
  } else {
    
    return(FALSE)
    
  }
}

# --------------------------------Create data for plots--------------------------------
# NOTE: TF_and_ext is a dataframe (loaded already) that created in order to identify 
# whether the TF data is a regular TF (with high confidence annotation) 
# or ext type(with lower confidence)
# Dogma: if we have regular TF type, we use that data; if we only have ext data, use ext
# if we have no data related to this tf, we will either give an error message or do nothing



#' 
#' @param tf character vector, containing one or more TF names
#' @param method either by Cell --> use cell data, or by cluster --> use cluster data, 
#' this should be a string indicating the column name
#' @param region cortex or pons
#' @param TF_and_ext TF_and_ext data, specific for cortex/pons
#' @param timepoint user input developmental timepoint, only used if looking at per sample data
#' @param per_sample boolean
#' @param bad_cells list of cells that belong to excluded/blacklisted clusters that needs to be filtered out 
#' @return a dataframe that has a column containing all the cell names and columns of the input tfs
#' the corresponding activity
#' data value (NES) 
#'
#' 
#' This function subsets a TF activity feather (basically a dataframe with cell/cluster identifiers as rows
#' and regulons as columns) based on the user TF input, region of interest, timepoint (if looking at 
#' per-sample data), and the cell/cluster identification method (methd can be joint for joint space cluster, 
#' cell for TF activity per cell, or cluster for TF activity per sample clustering)



create_activity_data <- function(tf, method, region, TF_and_ext,
                                 timepoint = NULL, per_sample = FALSE,
                                 bad_cells = ""){ 
  tp <- timepoint
  
  #building path of the feather file to read if the per sample toggle is on

  # set up the path of the feather file to read depending on region and cluster or cell
  if(per_sample == TRUE){
  
    reg <- switch(region, "cortex" = "ct", "pons" = "po")
  

    meth <- str_to_lower(method) #per_sample should only be true for method: "joint" or "Cell"
                                 #if its "Cell", turn to lower case and use for DR plots
                                 #if its "joint", then its for per-sample heatmap and load feather
                                   # with per_cluster
    
    if(identical(method, "joint")){
      meth <- "cluster"
    }
    
    path <- glue("data/{reg}_{timepoint}/{reg}_{timepoint}.regulon_activity_per_{meth}.feather")
  }
  # set up the path of the feather file to read depending on region and cluster or cell
  #per_sample is false here so read rether files for the joint_extended data
  else{
    
    if(!region %in% c("cortex", "pons")) return("Wrong usage: region should be either cortex/pons") #sanity check
    
    if(method == "joint"){
      path <- glue('data/joint_{region}_extended/joint_{region}_extended.regulon_activity_per_cluster.joint_extended.feather')
    }
    else if (method == "Cluster"){
      path <- glue('data/joint_{region}_extended/joint_{region}_extended.regulon_activity_per_cluster.per_sample.feather')
    }
    #method is Cell 
    else{
      path <- glue('data/joint_{region}_extended/joint_{region}_extended.regulon_activity_per_cell.feather')
    }
 

  }
  # case-insensitive checking and reading in the first column which corresponds to cell or cluster label
  # the first step just reads in the cell or cluster labels column -> collectively called cell_col
  #the TF columns are read after 
  if(str_detect(method,"(?i)Cell")){
    cell_col <- read_feather(path, "Cell")
  }
  else if(str_detect(method,"(?i)Cluster")){
    cell_col <- read_feather(path,"ID_20210710")
    colnames(cell_col) <- "Cluster"
  }
  else if(str_detect(method,"(?i)joint")){
    #per sample data has a different name for the cluster column compared to the joint_extended data
    if(per_sample){
      ID <- "ID_20201028"
      cell_col <- read_feather(path, ID)
      colnames(cell_col) <- "Cluster"
      method <- "Cluster"
    }
    else{
      cell_col <- read_feather(path, "ID_20210710_joint_clustering")
      colnames(cell_col) <- "Joint_cluster"
      method <- "Joint_cluster"
    }
  }
  else{return("Wrong usage, method should be Cell/Cluster")} #this should never return when called from the app
  
  
  
  #Read TF activity columns and add to the cell/cluster identifier column read above
  #loops through each factor in tf input and checks to see if the entry has regular or ext forms and extracts
  #data prioritizing the regular regulon data and not extended
  #puts each regular or ext TF name in a list 
  
  #activity is the matrix that will contain a subset of the entire feather depending on what TFs the user selected
  activity <- cell_col

  
  #buids tf_to_read vector: a list of regulon names that should be read from the feather file
  
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

    #reads the data from file specified in path using the tf_to_read list  
    
    col <- read_feather(path,tf_to_read)
  
    activity <- add_column(activity,col)
  }
  activity %>%
    select(method, everything()) # move the cell/cluster identifier column to start
  
  if(!(is.null(tp)) & per_sample == FALSE & identical(method, "Cluster")){ #when timepoint has an input and the input is not All
    #split column, select rows for the corresponding timepoints using filter 
    if(identical(tp, "F-All Time-Points") || identical(tp, "P-All Time-Points")){
    }
    else{
      activity <- separate(activity, Cluster, into = c("Timepoint", "Cluster"), sep = "_")
      activity <- filter(activity, Timepoint == tp)
      activity <- unite(activity, Cluster, Timepoint, Cluster, sep = "_")
    }
  }
  else if(identical(method, "Cell") & per_sample){ #gets rid of blacklisted cells in the per sample data

     activity <- activity %>%
       filter(!(Cell %in% bad_cells))
 
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
#' @param timepoint
#' @param per_sample 
#' 
#' params same as create_activity_data inputs because the function is called in this one
#'

#
#' 
plot_heatmap <- function(tf, method, region, TF_and_ext,
                         timepoint = NULL, per_sample = FALSE){
  # sanity checking
  if(!region %in% c("cortex", "pons")) return("Wrong usage: region should be either cortex/pons")
  if(!method %in% c("Cell","Cluster", "joint")) return("Wrong usage, method should be Cell/Cluster")
  
  
  #generates heatmap for the per sample clusters in the joint_extended datasets
  if(method == "Cluster"){
    act <- create_activity_data(tf, "Cluster",region, TF_and_ext, timepoint) 

    act <- column_to_rownames(act, var = "Cluster") # make that column name as row name ...
    
    # hm_anno$anno_row contains the cluster labels of all clusters in the joint_extended datasets
    # used to colour the cluster by a palette
    new_anno_row <- hm_anno$anno_row 
    
    rownames(new_anno_row) <- gsub(" ", "_", rownames(hm_anno$anno_row)) # re-assign the rownames
    # note that the rownames correspond to the col names of the matrix t(act_cluster)
    # customized for plotting by cluster
    
 
    anno_col <- new_anno_row %>%
      mutate('Broad Cluster' = recode(rownames(new_anno_row), !!!lvl2_cluster_labels)) %>%
      mutate('Broader Cluster' = recode(rownames(new_anno_row), !!!lvl1_cluster_labels))
    
    rownames(anno_col)

    show_colname_plot <- TRUE
    title <- glue('Transcription Factor Regulon Activity at Developmental Time: {timepoint}')
  }
  else if(method == "joint"){ #plot heat map by joint cluster for the joint_extended dataset, or for the per-sample analyses
    
    act <- create_activity_data(tf, "joint", region, TF_and_ext, per_sample = per_sample,
                                timepoint = timepoint) 
    
    if(per_sample == TRUE){
      #filter out the exclude clusters
      act <- act %>% filter(!grepl("EXCLUDE", act$Cluster))

      col_to_row <- "Cluster"
      
      
      #generate the labels for the clusters in the data-set so that the palette can be properly displayed
      new_anno_row <- act %>% mutate(rownames = Cluster) %>%
        column_to_rownames("rownames") %>% select(Cluster) %>%
        mutate('Broad Cluster' = recode(act$Cluster, !!!lvl2_cluster_labels)) %>%
        mutate('Broader Cluster' = recode(act$Cluster, !!!lvl1_cluster_labels))
      
      
      
    }
    else{
      col_to_row <- "Joint_cluster"
      new_anno_row <- act %>% mutate(Cluster = Joint_cluster) %>%
        column_to_rownames("Joint_cluster") %>% select(Cluster)
    }
  
    act <- act %>%
      column_to_rownames(var = col_to_row) 

    # note that the rownames correspond to the col names of the matrix t(act_cluster)
    # customized for plotting by cluster
    anno_col <- new_anno_row # this is loaded by data_prep.R

    show_colname_plot <- TRUE
    title <- "Transcription Factor Regulon Activity per Cluster"
  }
  cluster_row <- FALSE
 
  #do not do row heirarchal clustering if there is only one TF selected
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
                     cellwidth = 20,
                     cellheight = 20)
} 

#----------------------------Dimension reduction---------------------------------------------
#' Make DR clustering scatterplot
#'
#' @param tf_number Either 1 or 2. In the tf input vector we get from user in Shiny app, there could be
#' multiple tfs, but we only support plotting two tfs since the scatterplot is big
#' @param cell_metadata metadata containing DR coordinates , saved in data_cortex
#' and data_pons
#' @param cell_activity_data made by create_activity_data() function given the tf input
#' @param dim_red_type allow user to select between UMAP, tSNE, PCA DR methods
#'
#' @return a DR scatter plot coloured by TF activity
#'

plot_dr <- function(tf_number, cell_metadata, cell_activity_data, dim_red_type){ #cell_metadata is the tsv with the 
  #embedding coordinates for each cell in the dataset

  
  tf_plot <- tf_number + 1 #replaces the above control flow

  
  activity_tf <- cell_activity_data[,tf_plot][[1]] #extracts the TF activity from the cell_activity_data
  #and appends it to the cell_metadata to make cell meta with activity to plot
  
  TF <- colnames(cell_activity_data)[tf_plot]
  
  cell_meta_with_activity <- mutate(cell_metadata, activity_tf = activity_tf)
  
  x_axis <- switch(dim_red_type, "umap" = "UMAP1", "tsne" = "tSNE_1", "pca" = "PC1")
  y_axis <- switch(dim_red_type, "umap" = "UMAP2", "tsne" = "tSNE_2", "pca" = "PC2")
  
  ggplot(data = cell_meta_with_activity, mapping = aes_string(x = x_axis, y = y_axis))+
    geom_point(aes(color = activity_tf), alpha = 0.2)+
    scale_color_gradient(low = "grey90", high = "red")+
    theme_min() + labs(color = 'TF Activity') + ggtitle(TF)
  
}

#' @param cluster_label boolean: include ggrepel label of clusters or not
#' @param cell_metadata metadata containing DR coordinates
#' @param per_sample indicates if the function call is for a per_sample dataset 
#' @param dim_red_type allow user to select between UMAP, tSNE, PCA DR methods
#'
#' @return a DR scatter plot coloured by cluster label

color_by_cluster <- function(cell_metadata, dim_red_type, cluster_label, 
                             per_sample = FALSE){
  
  x_axis <- switch(dim_red_type, "umap" = "UMAP1", "tsne" = "tSNE_1", "pca" = "PC1")
  y_axis <- switch(dim_red_type, "umap" = "UMAP2", "tsne" = "tSNE_2", "pca" = "PC2")
  
  # Store the center points (medians) of each cluster and add labels courtesy of Bhavyaa 
  
  var_group <- "Joint_cluster"
  if(per_sample){
    var_group <- "ID_20201028_with_exclude"
  }
  
  gg <- ggplot(data = cell_metadata, mapping = aes_string(x = x_axis,y = y_axis))+
    geom_point(aes(color = get(var_group)), alpha = 0.2) + theme_min() + theme(legend.position="none") +
    scale_color_manual(values = master_palette$Cluster)
  
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
  return(gg)
}

#----------------------------Time Course Ribbon Plot---------------------------------------------


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
create_metadata_timeseries <- function(cell_metadata, part, general_cluster_labels){
  if(part == "cortex") level <- c("Forebrain E10.5",
                                  "Forebrain E12.5",
                                  "Forebrain E13.5",
                                  "Forebrain E15.5",
                                  "Forebrain E16.5",
                                  "Forebrain E18.5",
                                  "Forebrain P0",
                                  "Forebrain P3",
                                  "Forebrain P6")
  else if (part == "pons") level <- c("Hindbrain E10.5",
                                      "Hindbrain E12.5",
                                      "Hindbrain E13.5",
                                      "Pons E15.5",
                                      "Pons E16.5",
                                      "Pons E18.5",
                                      "Pons P0",
                                      "Pons P3",
                                      "Pons P6")
  else{(return("Wrong usage, input either cortex or pons"))} #sanity check
  
  #print(level)
  cell_metadata %>% 
    select(Age = Sample, Cell, Cluster = Sample_cluster) %>% 
    # In this case, we remove the "prefix" of the Cluster column, so that we are
    # simply left with the abbreviation representing the cell type, so that 
    # we can link the cells of the same cell type across ages
    #separate(Cluster, into = c("Prefix", "Cluster"), sep = "_") %>% 
    #filter(!grepl("EXCLUDE", Cluster)) %>%
    mutate(broad_cluster = recode(Cluster, !!!general_cluster_labels)) %>%
    mutate(Age = factor(Age, levels = level)) %>% 
    arrange(Cell)
  
}


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
plot_timeseries <- function(TF, cell_metadata, activity, make_plotly = FALSE, show_legend = TRUE){
  
  cell_names <- colnames(activity)
  
  activity <- activity[TF, ] %>% 
    {data.frame("TF" = .)} 
  
  activity$Cell <- cell_names
  
  activity <- activity %>% arrange(Cell)
  
  #there is one cell in the pons metadata file that is not in the binary activity matrix so im gonna remove it
  #not sure why its here
  
  
  cell_diff <- setdiff(cell_metadata$Cell, activity$Cell)
  if(length(cell_diff) > 0) {cell_metadata <- cell_metadata %>% filter(cell_metadata$Cell != cell_diff)}
  
  
  
  if(!all(cell_metadata$Cell == activity$Cell)) return (-1)
  # Add the binarized TF activity to the new dataframe
  ribbon_df <- cell_metadata
  ribbon_df$TF <- activity$TF
  
  ribbon_df <- ribbon_df %>% 
    filter(!grepl("EXCLUDE", Cluster)) %>%
    select(-Cluster) 
  
  ribbon_df_celltype_frac <- ribbon_df %>% 
    group_by(Age) %>% 
    # Total cells at each age
    mutate(total = n()) %>% 
    group_by(Age, broad_cluster) %>%
    # Proportion of TF+ cells per cluster, per age
    mutate(frac = sum(TF > 0) / total) %>% 
    distinct(Age, broad_cluster, frac) %>% 
    ungroup() %>% group_by(Age)
  
  ribbon_df_clusters_complete <- ribbon_df_celltype_frac %>%
    mutate(broad_cluster = factor(broad_cluster, levels = unique(.$broad_cluster))) %>%
    complete(broad_cluster, nesting(Age), fill = list(frac = 0))
  
  ribbon_df_clusters_complete$xpos = group_indices(ribbon_df_clusters_complete)
  
  plot <- ribbon_df_clusters_complete %>%
    ggplot(aes(x = xpos, y = frac, fill = broad_cluster)) +
    geom_area(stat = "identity", show.legend = show_legend) +
    scale_fill_manual(values = palette_broad_clusters, drop = FALSE, name = "") +
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9),
                       labels = c("E10.5", "E12.5", "E13.5", "E15.5", "E16.5", "E18.5", "P0", "P3", "P6"),
                       limits = c(1, 9)) +
    labs(x = "Developmental Age", y = "Proportion", title = TF) +
    guides(fill = guide_legend(ncol = 5)) +
    theme_min() + 
    theme(legend.position = "bottom") 
  
  if(make_plotly) {
    return (ggplotly(plot, tooltip = "broad_cluster") %>% style(hoveron = "points + fills"))
  }
  else{return(plot)}
}
#--------------------------Active_specific------------------------

#' @param sample name of the data-set that the user wants to look at, used as a part of the path of file to read
#' @param cluster cluster of interest 

#'
#' @return a list of 3 matrices: tf_table corresponds to the table for a specific cluster and used in the By_cluster subtab
#' FC_df and AUC_df are matriices containing AUC and FC values for each gene in each cluster   

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
}

#' @param data data generated by active_specific_prep. Specifically, the tf_table part
#' @param fc Fold change cutoff selected by the user, TFs above cutoff are labelled and displayed in the table, genes below
#' are not
#' @param cluster cluster of interest 

#'
#' @return a scatter plot of TFs

plot_scatter <- function(data, fc, cluster){

  #data used for ggrepel
  #the mutate(gsub) statements are to replace the "_extended (21g)" portion with just a "+" to shorten the label
  to_label <- data %>% filter(AUC_FC > fc) %>% mutate(why = gsub("_extended", "+", TF)) %>%
    mutate(TF = gsub("\\(.+\\)", "", why))

  
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

#' @param data data generated by active_specific prep. Specifically, the AUC_df data frame
#' @param tf cluster of interest 

#'
#' @return a list of bar plots that displays the TF ranked by AUC value in each cluster 
plot_bar_list <- function(data, tf){
  
  data <- data %>% mutate(Cluster = gsub("_", " ", Cluster))
    
  palette <- master_palette$Cluster[names(master_palette$Cluster) %in% data$Cluster]

  
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

#-------------Bubble-plot--------------------
#code adapeted from the clusters app

bubble_prep <- function(sample, tf, dend_order, scale){
  path <- glue("data/{sample}/{sample}.active_specific_tf.Rds")
  data <- readRDS(path)
  
  AUC <- data$AUC_df %>% select(Cluster, tf) %>%
    filter(Cluster %in% dend_order)
  
  FC <- data$FC_df %>% select(Cluster, tf) %>%
    filter(Cluster %in% dend_order)
  
  # Scale activity of each tf linearly across clusters to [0, 1]
  if (scale) {
    
    AUC <- AUC %>%
      as.data.frame() %>%
      select(-Cluster) %>%
      apply(2, scales::rescale, to = c(0, 1)) %>% 
      as.data.frame %>%
      mutate(Cluster = rownames(AUC)) %>% 
      select(Cluster, everything())
    
    # FC <- FC %>% 
    #   as.data.frame() %>%
    #   select(-Cluster) %>%
    #   apply(c(1,2), log2) %>%
    #   as.data.frame %>%
    #   mutate(Cluster = rownames(AUC)) %>% 
    #   select(Cluster, everything())
    #   
    
  }
  
  
  # Convert to long / tidy format with columns: Cluster, TF, AUC
  AUC <- AUC %>% 
    gather(., "TF", "AUC", 2:ncol(.))
  
  #print(AUC)
  
  FC <- FC %>%
    gather(., "TF", "FC", 2:ncol(.))
  
  #print(AUC)
  
  # print(FC)
  
  df <- left_join(AUC, FC, by = c("Cluster", "TF"))
  #print(df)
  
  # Tidy data for plotting
  df <- df %>%
    
    # Order genes to match order input by user
    mutate(TF = factor(TF, levels = rev(tf))) %>% 
    arrange(TF) %>% 
    
    # Pad gene names so that the plot takes up a more standardized
    # width; to roughly the the # of characters in the gene w/ longest name
    # However, letters take up more pixels than spaces, so do less padding
    # for genes with longer names
    # TODO: Fix alignment of bubble plot w/ dendrogram for long gene names (issue #7)
    mutate(TF_padded = case_when(
      str_length(TF) <= 5 ~ str_pad(TF, 15, side = 'right', pad = " "),
      str_length(TF) > 5 ~ str_pad(TF, 12, side = 'right', pad = " ")
    )
    ) %>% 
    mutate(TF_padded = factor(TF_padded, levels = unique(.$TF_padded))) %>% 
    
    # Order the clusters on the x-axis to match the dendrogram image
    #slice(match(dend_order_joint_cortex_extended, Cluster)) %>%
    
    # Keep columns
    select(TF, Cluster, AUC, FC, TF_padded)
  
 # df$Cluster <- factor(df$Cluster, levels = dend_order)
  
  label_palette <- hm_anno$side_colors$Cluster
  names(label_palette) <- gsub(" ", "_", names(hm_anno$side_colors$Cluster))
  label_palette <- label_palette[dend_order]

  
  
  return(list("data" = df, "label_palette" = label_palette))
}

plot_bubble <- function(data, label_palette, dend_order){
  
  data$Cluster <- factor(data$Cluster, levels=dend_order)
  
  p1 <- data %>% 
    ggplot(aes(x = Cluster, y = TF_padded)) +
    geom_point(aes(size = FC, colour = AUC), alpha = 0.8) +
    scale_size_area(max_size = 4) +
    scale_color_gradientn(colours = tail(rdbu, 70)) +
    theme_min() +
    ylab(NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                                     colour = label_palette,
                                     size = rel(0.7)),
          panel.grid.major.x = element_line(colour = "grey90"),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key.width = unit(1, "cm"),
          legend.position = "bottom") +
    # Put gene labels on the right hand side to improve alignment
    scale_y_discrete(position = "right")
  
  gene_labels <- cowplot::get_y_axis(plot = p1, position = "right")
  
  p1 <- p1 + scale_y_discrete(labels = NULL)
  
  return(list(plot = p1, labels = gene_labels))
}




