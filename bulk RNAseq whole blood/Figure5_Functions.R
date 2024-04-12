### Generate norm_anno 

generate_norm_anno <- function(dds_object = dds){
  
  norm_anno <- as.data.frame(counts(dds_object, normalized=T))
  norm_anno$GENEID <- row.names(norm_anno)
  
  # add gene annotation extracted from the gtf file
  gene_annotation <- tx_annotation[!duplicated(tx_annotation$GENEID),c("GENEID", "SYMBOL", "GENETYPE")]
  gene_annotation <- gene_annotation[match(rownames(norm_anno), gene_annotation$GENEID), ]
  
  # # check if row names of the normalized table and the gene annotation match perfectly
  # all(rownames(norm_anno) == gene_annotation$GENEID)
  
  # add additional gene annotation downloaded from biomart
  mart <- biomaRt::useMart(host = "https://may2021.archive.ensembl.org",
                           biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl")
  biomart<- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description",
                                 "chromosome_name", "transcript_length", "gene_biotype"), mart = mart)
  colnames(biomart) <- c("GENEID", "Gene.stable.ID", "Gene.description", "Chromosome.scaffold.name", "length", "genetype")
  
  # biomart <- read.delim(file.path(dir, "Data", "biomart_180914.txt"), stringsAsFactors = FALSE)
  idx <- match(unlist(lapply(strsplit(gene_annotation$GENEID, split = "[.]"), `[[`, 1)), biomart$Gene.stable.ID)
  gene_annotation$DESCRIPTION <- biomart$Gene.description[idx]
  gene_annotation$CHR <- biomart$Chromosome.scaffold.name[idx]
  
  # merge expression table and annotation
  norm_anno <- merge(norm_anno,
                     gene_annotation,
                     by = "GENEID")
  rownames(norm_anno) <- norm_anno$GENEID
  
  tmp <- list("gene_annotation" = gene_annotation, 
              "norm_anno" = norm_anno)
  return(tmp)
}

### PCA function

plotPCA <- function(pca_input = dds_vst,
                    pca_sample_table = sample_table,
                    ntop=500,
                    xPC=1,
                    yPC=2,
                    color,
                    anno_colour,
                    shape="NULL",
                    point_size=3,
                    title="PCA",
                    label = NULL,
                    label_subset = NULL){
  
  if(is.character(pca_input)){
    vst_matrix <- as.matrix(removedbatch_dds_vst)
  }else if(!is.data.frame(pca_input)){
    vst_matrix <- as.matrix(assay(pca_input))
  }else{
    vst_matrix <- pca_input
  }
  
  if(ntop=="all"){
    pca <- prcomp(t(vst_matrix))
  }else{
    # select the ntop genes by variance
    select <- order(rowVars(vst_matrix), decreasing=TRUE)[c(1:ntop)]
    pca <- prcomp(t(vst_matrix[select,]))
  }
  
  #calculate explained variance per PC
  explVar <- pca$sdev^2/sum(pca$sdev^2)
  # transform variance to percent
  percentVar <- round(100 * explVar[c(xPC,yPC)], digits=1)
  
  # Define data for plotting
  pcaData <- data.frame(xPC=pca$x[,xPC],
                        yPC=pca$x[,yPC],
                        color = pca_sample_table[[color]],
                        name= as.character(pca_sample_table$ID),
                        stringsAsFactors = F)
  
  #plot PCA
  if(is.factor(pcaData$color) || is.character(pcaData$color)|| is.integer(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        stat_ellipse(aes(fill = color), geom = "polygon", alpha = 0.05) +#
        geom_point(size =point_size) 
    }else{
      pcaData$shape = pca_sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
        stat_ellipse(geom = "polygon", aes(fill = color), alpha = 0.05) +
        geom_point(size =point_size) +
        scale_shape_discrete(name=shape)
      
    }
    
    if(anno_colour[1] == "NULL"){
      pca_plot <- pca_plot + scale_color_discrete(name=color)
    }else{
      pca_plot <- pca_plot + scale_color_manual(values=anno_colour, name=color, aesthetics = c("colour", "fill"))
    }
    
  }else if(is.numeric(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color) +
        scale_fill_gradientn(colours = bluered(100),name=color)
    }else{
      pcaData$shape = pca_sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color) +
        scale_fill_gradientn(colours = bluered(100),name=color) +
        scale_shape_discrete(name=shape)
    }
  }
  
  # adds a label to the plot. To label only specific points, put them in the arument label_subset
  if (!is.null(label) == TRUE){
    pcaData$label <- pca_sample_table[[label]]
    if(!is.null(label_subset) == TRUE){
      pcaData_labeled <- pcaData[pcaData$label %in% label_subset,]
    } else {
      pcaData_labeled <- pcaData
    }
    pca_plot <- pca_plot +
      geom_text_repel(data = pcaData_labeled, aes(label = label), nudge_x = 2, nudge_y = 2, colour = "black")
  }
  
  pca_plot <- pca_plot+
    xlab(paste0("PC ",xPC, ": ", percentVar[1], "% variance")) +
    ylab(paste0("PC ",yPC,": ", percentVar[2], "% variance")) +
    coord_fixed()+
    theme_bw()+
    theme(aspect.ratio = 1,
          legend.position = "bottom")+
    ggtitle(title)
  
  pca_plot
  ggMarginal(pca_plot, groupFill = T, type = "density")
}

#isolated GSVA stat function
GSVA_stats <- function(input = norm_anno,
                       sample_table,
                       voi,
                       signatures,
                       method = "gsva",
                       kcdf = "Gaussian",
                       mx.diff = T,
                       abs.ranking = F,
                       statistics = TRUE, 
                       comparisons = my_comparisons,
                       plot = F){
  
  ## Select reference gene set & create signature list
  signature_list <- list()
  for (i in 1:length(signatures)) {
    print(i)
    tmp<-input[input$SYMBOL%in%signatures[[i]],]$GENEID
    signature_list[[i]]<-tmp
  }
  # Reorder sample table to match the order of the count matrix
  count.matrix <- input[, as.character(colnames(input)) %in% as.character(sample_table[["ID"]])]
  idx <- match(as.character(colnames(count.matrix)), as.character(sample_table[["ID"]]))
  sample_table_sorted <- sample_table[ idx,]
  # Ensure order are correct
  if (identical(colnames(count.matrix), as.character(sample_table_sorted[["ID"]]))){
    # Rename columns
    colnames(count.matrix) <- paste0(sample_table_sorted[["ID"]], ".", sample_table_sorted[[voi]])
    count.matrix <- as.matrix(count.matrix)
    
    
    
    ## Perform GSVA
    result_GSVA <- gsva(expr = count.matrix, 
                        gset.idx.list = signature_list, 
                        method = method, 
                        mx.diff = mx.diff,
                        abs.ranking = abs.ranking,
                        kcdf = kcdf, 
                        verbose = F, 
                        parallel.sz = 1)
    result_GSVA_melted <- melt(result_GSVA)
    
    result_GSVA_melted[["ID"]] <- sapply(result_GSVA_melted$Var2, function(x){
      unlist(strsplit(as.character(x), split = "\\."))[1]
    })
    
    result_GSVA_melted[[voi]] <- sapply(result_GSVA_melted$Var2, function(x){
      unlist(strsplit(as.character(x), split = "\\."))[2]
    })
    
    stats.df <-result_GSVA_melted %>%
      group_by(Var1) %>%
      rstatix::wilcox_test(as.formula(paste0("value", "~", voi)), comparisons = comparisons) %>% 
      rstatix::adjust_pvalue(method = "BH") %>%
      rstatix::add_significance() %>%
      rstatix::add_xy_position(x = voi)
    stats.df$comp <- apply(stats.df[,c("group1","group2")], 1, paste, collapse = "_vs_")
    stats.df$comp_rev <- apply(stats.df[,c("group2","group1")], 1, paste, collapse = "_vs_")
    
    if (plot) {
      p <- ggplot(result_GSVA_melted, aes_string(x = voi, y = "value")) +
        geom_boxplot(aes_string(fill = voi),alpha = 0.5, outlier.colour = "white", outlier.size = 0.1) + 
        geom_jitter(aes_string(color = voi),width = 0.3, size = 2) + 
        scale_fill_manual(values = col_condition) + 
        scale_color_manual(values = col_condition) + 
        facet_wrap(~Var1, scales = "free") + 
        theme_light() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.line = element_line(colour = "black"), 
              strip.text = element_text(color = "black"), 
              axis.text.x=element_text(angle=90,  vjust=0.5, hjust = 1),
              text = element_text(size = 12, color = "black")
        ) +
        ylab("Enrichment score") + 
        xlab("") +
        ggpubr::stat_pvalue_manual(stats.df,  label = "p.adj", tip.length = 0.01)+
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
      
      print(p)
    }
    
    return(stats.df)
  }
}

#wrap stat function in loop of signature with different length
GSVA_signature_optimization <- function(input_gsva = norm_anno,
                                        sample_table_gsva = sample_table,
                                        voi,
                                        signatures,
                                        comparisons,
                                        start = 10,
                                        step.size = 1){
  
  signature_stats <- lapply(signatures,function(x){
    
    #compute the number of iterations (plus 1 for starting point)
    num_iterations <- ceiling((length(x)-start)/step.size) + 1
    
    #generate a data frame to save the results
    output <- matrix(ncol = length(comparisons)+1, nrow = num_iterations) %>% as.data.frame()
    colnames(output) <- c("gene_n",
                          unlist(lapply(comparisons, paste, collapse = "_vs_")))
    
    #compute the stats for each iteration
    for (i in 1:num_iterations){
      
      #compute the number of signature genes to use
      genes_to_use <- ifelse((start + i*step.size -1) <= length(x), (start + i*step.size -1), length(x))
      
      tmp <- GSVA_stats(input = input_gsva,
                        sample_table = sample_table_gsva,
                        signatures = list(x[1:genes_to_use]),
                        voi = voi,
                        method = "gsva",
                        kcdf = "Gaussian",
                        statistics = TRUE,
                        comparisons = comparisons
      )
      #add info to output
      output[i,"gene_n"] <- genes_to_use
      for (k in 1:length(comparisons)){
        output[i,1+k] <- tmp[tmp$comp == colnames(output)[1+k] | tmp$comp_rev == colnames(output)[1+k], "p.adj"]
      }
    }
    output
  })
  signature_stats
}

#wrap stat function in loop of random signatures to create signature specific probability distribution
GSVA_emperical_distribution <- function(input_list,
                                        sample_table_list,
                                        voi,
                                        comparisons,
                                        signature,
                                        iterations = 500,
                                        seed = 42){
  
  #set the seed
  set.seed(seed)
  
  #generate a data frame to save the results
  output <- lapply(1:length(input_list), function(x){
    tmp <- matrix(ncol = length(comparisons)+1, nrow = iterations) %>% as.data.frame()
    colnames(tmp) <- c("iteration",
                       unlist(lapply(comparisons, paste, collapse = "_vs_")))
    tmp
  })
  names(output) <- names(input_list)
  
  
  #get the intersection of all available genes across data sets
  gene_union <- input_list[[1]]$SYMBOL
  #subset the signature to the available to genes
  shared_signature <- signature[signature %in% gene_union]
  
  #compute the signature.size
  signature.size <- length(shared_signature)
  
  #compute a list of used signatures
  previous_signatures <- list(shared_signature)
  
  #compute the stats for each iteration
  for (i in 1:iterations){
    
    #compute the random signature
    genes_to_use <- sample(gene_union, size = signature.size)
    
    #check whether the signature was used before
    while(!all(is.na(unlist(lapply(previous_signatures, function(x) all(match(x, genes_to_use))))))){
      genes_to_use <- sample(gene_union, size = signature.size)
    }
    
    previous_signatures[[length(previous_signatures)+1]] <- genes_to_use
    
    for (n in 1:length(input_list)){
      
      tmp <- GSVA_stats(input = input_list[[n]],
                        sample_table = sample_table_list[[n]],
                        signatures = list(genes_to_use),
                        voi = voi,
                        method = "gsva",
                        kcdf = "Gaussian",
                        statistics = TRUE,
                        comparisons = comparisons
      )
      #add info to ouput
      output[[n]][i,"iteration"] <- i
      print(paste0("Iteration: ",i))
      for (k in 1:length(comparisons)){
        output[[n]][i,1+k] <- tmp[tmp$comp == colnames(output[[n]])[1+k] | tmp$comp_rev == colnames(output[[n]])[1+k], "p.adj"]
      }
    }
  }
  output
  
}

#investigate the signature edge
GSVA_signature_edge <- function(input_list,
                                sample_table_list,
                                signature){
  
  #assign input list to variable with intersection of genes from both data sets
  expr_data <- input_list
  
  #generate a rank matrices
  rank_list <- lapply(1:length(expr_data), function(x){
    
    gene_number <- nrow(expr_data[[x]])
    sample_number <- nrow(sample_table_list[[x]])
    
    #initiate C script from GSVA
    A = .C("matrix_density_R",
           as.double(t(as.matrix(expr_data[[x]][,as.character(sample_table_list[[x]]$ID), drop = F]))),
           as.double(t(as.matrix(expr_data[[x]][,as.character(sample_table_list[[x]]$ID)]))),
           R = double(sample_number * gene_number),
           as.integer(length(sample_table_list[[x]]$ID)),
           as.integer(sample_number),
           as.integer(gene_number),
           as.integer(F))$R
    
    gene.density <- t(matrix(A, 
                             sample_number, 
                             gene_number))
    
    rank.scores <- as.data.frame(gene.density)
    #compute rank
    colnames(rank.scores) <- colnames(expr_data[[x]][,as.character(sample_table_list[[x]]$ID)])
    rank.scores$SYMBOL <- expr_data[[x]]$SYMBOL
    #rank the mean scores
    rank.scores$mean <- apply(rank.scores[, colnames(rank.scores) != "SYMBOL" & colnames(rank.scores) %in% sample_table_list[[x]][sample_table_list[[x]]$condition == "severe_deceased",]$ID],1,mean)
    rank.scores <- rank.scores %>% arrange(desc(mean))
    rank.scores$rank <- 1:gene_number
    
    #compute a gsea of the mean rank
    tmp_rank <- -(rank.scores$rank-nrow(rank.scores)/2)
    names(tmp_rank) <- rank.scores$SYMBOL
    #visualize the result
    signature_list <- list(signature = signature)
    plot(fgsea::plotEnrichment(signature_list[["signature"]], tmp_rank) + theme_bw())
    fgsea <- fgsea::fgsea(signature_list[1], tmp_rank)
    #add leading edge info
    rank.scores$lead <- sapply(rank.scores$SYMBOL, function(y){
      if (y %in% unlist(fgsea$leadingEdge)) "edge"
      else NA
    })
    
    list(rank.scores[rank.scores$SYMBOL %in% signature,] %>% dplyr::select(SYMBOL, lead, rank, mean, everything()),
         fgsea)
    
  })
  
  
  return(rank_list)
  
}


### Heatmap colors
scaleColors <- function(data = input_scale, # data to use
                        maxvalue = NULL # value at which the color is fully red / blue
){
  if(is.null(maxvalue)){
    maxvalue <- floor(min(abs(min(data)), max(data)))
  }
  if(max(data) > abs(min(data))){
    if(ceiling(max(data)) == maxvalue){
      myBreaks <- c(floor(-max(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(max(data)))
    } else{
      myBreaks <- c(floor(-max(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(max(data)))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  } else {
    if(-floor(min(data)) == maxvalue){
      myBreaks <- c(floor(min(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(min(data)))
    } else{
      myBreaks <- c(floor(min(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(abs(min(data))))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  }
  return(list(breaks = myBreaks, color = myColor))
}

### Heatmap

plotHeatmap <- function(input = norm_anno,
                        geneset = "all",
                        title = "",
                        keyType = "Ensembl",
                        gene_type = "all",
                        show_rownames = FALSE,
                        cluster_cols = FALSE,
                        cluster_rows = TRUE,
                        method = "complete",
                        dist_method = "euclidean",
                        sample_annotation = sample_table,
                        plot_mean = FALSE)
{
  if (geneset[1] != "all") {
    if (keyType == "Ensembl") {
      input <- input[input$GENEID %in% geneset, ]
    }
    else if (keyType == "Symbol") {
      input <- input[input$SYMBOL %in% geneset, ]
    }
    else {
      print("Wrong keyType. Choose Ensembl or Symbol!")
    }
  }
  
  
  
  if (gene_type != "all") {
    input <- input[input$GENETYPE %in% gene_type, ]
  }
  
  rownames(input) <- paste(input$SYMBOL) 
  
  if(plot_mean == FALSE ){
    input <- input[ , colnames(input) %in% sample_annotation$ID]
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[, order(sample_annotation[[plot_order]], decreasing = FALSE)]
    
    column_annotation <- plot_annotation
    
  } else {
    input <- input[ , colnames(input) %in% plot_annotation_mean$condition]
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[, order(plot_annotation_mean$condition, decreasing = FALSE)]
    
    column_annotation <- plot_annotation_mean
    test <<- input_scale
  }
  
  plot <-  pheatmap(input_scale, main = title,
                    show_rownames = show_rownames,
                    show_colnames = TRUE,
                    cluster_cols = cluster_cols,
                    cluster_rows = cluster_rows,
                    fontsize = 7,
                    clustering_method = method,
                    clustering_distance_rows = dist_method,
                    clusering_distance_cols = dist_method,
                    annotation_col = column_annotation,
                    annotation_colors = ann_colors,
                    breaks = scaleColors(data = input_scale, maxvalue = 1.5)[["breaks"]],
                    color = scaleColors(data = input_scale, maxvalue = 1.5)[["color"]])
  plot
  return(plot)
}