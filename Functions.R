GOEA.of.list<- function(seurat_object,
                        condition.list,
                        GeneSets =c("GO","KEGG","Hallmark", "Reactome"),
                        GOntology = "BP",
                        pCorrection = "bonferroni",      # choose the p-value adjustment method (Steffi: bonferroni)
                        pvalueCutoff = 0.1,      # set the unadj. or adj. p-value cutoff (depending on correction method)
                        qvalueCutoff = 0.1       # set the q-value cutoff (FDR corrected)
){
  
  OrgDb = org.Hs.eg.db
  org = "hsa"
  
  
  #present genes
  present_genes<-rownames(seurat_object)
  #present_genes <- as.matrix(GetAssayData(object = seurat_object, slot = 'counts'))
  #present_genes <- rownames(present_genes[apply(present_genes, 1, function (x) {sum(x >= 1)})>0,])
  present_genes_entrez <- bitr(present_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  # output lists
  if("GO" %in% GeneSets){
    GO_list<-list()
  }
  if("KEGG" %in% GeneSets){
    KEGG_list<-list()
  }
  if("Hallmark" %in% GeneSets){
    HALLMARK_list<-list()
  }
  if("Reactome" %in% GeneSets){
    REACTOME_list<-list()
  }
  if("C_PATHWAYS" %in% GeneSets){
    C_PATHWAYS_list<-list()
  }
  
  
  #go through different Gene lists
  
  for(i in names(condition.list)){
    #print(i)
    #genes of interest
    genes_oi<-condition.list[[i]]
    if(i == "Plasmablasts_mod_Dex_vs_ctrl_down") entrez_genes_oi <- ""
    else   entrez_genes_oi   <- bitr(genes_oi, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
    
    
    # GO enrichment
    if("GO" %in% GeneSets){
      #print("Performing GO enrichment")
      GO <- as.data.frame(enrichGO(gene = entrez_genes_oi,
                                   universe = present_genes_entrez,
                                   OrgDb = OrgDb,
                                   ont = GOntology,
                                   pAdjustMethod = pCorrection,
                                   pvalueCutoff  = pvalueCutoff,
                                   qvalueCutoff  = qvalueCutoff,
                                   readable      = T))
      if(nrow(GO)>0){GO$Condition <- paste(i)}
    }
    
    # KEGG enrichment
    if("KEGG" %in% GeneSets){
      #print("Performing KEGG enrichment")
      
      org = "hsa"
      
      KEGG <- as.data.frame(enricher(entrez_genes_oi,
                                     TERM2GENE=KEGG_genes,
                                     universe = present_genes_entrez,
                                     pAdjustMethod = pCorrection,
                                     pvalueCutoff  = pvalueCutoff,
                                     qvalueCutoff = qvalueCutoff))
      if(nrow(KEGG)>0){KEGG$Condition <- paste(i)}
    }
    
    
    # Hallmark enrichment
    if("Hallmark" %in% GeneSets){
      #print("Performing Hallmark enrichment")
      
      HALLMARK <- as.data.frame(enricher(entrez_genes_oi,
                                         TERM2GENE=hallmark_genes,
                                         universe = present_genes_entrez,
                                         pAdjustMethod = pCorrection,
                                         pvalueCutoff  = pvalueCutoff,
                                         qvalueCutoff = qvalueCutoff))
      if(nrow(HALLMARK)>0){HALLMARK$Condition <- paste(i)}
    }
    
    
    # Reactome enrichment
    if("Reactome" %in% GeneSets){
      #print("Performing Reactome signature enrichment")
      
      REACTOME <- as.data.frame(enricher(entrez_genes_oi,
                                         TERM2GENE=reactome_genes,
                                         universe = present_genes_entrez,
                                         pAdjustMethod = pCorrection,
                                         pvalueCutoff  = pvalueCutoff,
                                         qvalueCutoff = qvalueCutoff))
      if(nrow(REACTOME)>0){REACTOME$Condition <- paste(i)}
    }
    
    # Canonical pathway enrichment
    if("C_PATHWAYS" %in% GeneSets){
      #print("Performing CANNONICAL PATHWAYS signature enrichment")
      
      C_PATHWAYS <- as.data.frame(enricher(entrez_genes_oi,
                                           TERM2GENE=cannonicalPathway_genes_symbols,
                                           universe = present_genes_entrez,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff = qvalueCutoff))
      if(nrow(C_PATHWAYS)>0){C_PATHWAYS$Condition <- paste(i)}
    }   
    
    #concatenate output
    if("GO" %in% GeneSets){
      GO_list[[paste(i)]]<-GO
    }
    if("KEGG" %in% GeneSets){
      KEGG_list[[paste(i)]]<-KEGG
    }
    if("Hallmark" %in% GeneSets){
      HALLMARK_list[[paste(i)]]<-HALLMARK
    }
    if("Reactome" %in% GeneSets){
      REACTOME_list[[paste(i)]]<-REACTOME
    }
    if("C_PATHWAYS" %in% GeneSets){
      C_PATHWAYS_list[[paste(i)]]<-C_PATHWAYS
    }
    
  }
  results<-list()
  results$GOresults <- do.call(rbind, GO_list)
  results$KEGGresults <- do.call(rbind, KEGG_list)
  results$HALLMARKresults <- do.call(rbind, HALLMARK_list)
  results$REACTOMEresults <- do.call(rbind, REACTOME_list)
  # results$C_PATHWAYSresults <- do.call(rbind, C_PATHWAYS_list)
  
  return(results)
}

dotplot.GOEA.list<-function(x,
                            show =10,
                            orderBy = "count", #padj or count
                            colorBy = c("padj", "regulation") #padj or regulation
){
  if(nrow(x)<1){
    print("No enrichment found.")
  }else{
    
    
    
    #ORDERING
    if(orderBy=="padj"){
      x %>% group_by(Condition) %>% arrange(desc(Count), .by_group = TRUE) ->x
      x %>% group_by(Condition) %>% arrange(desc(Count), .by_group = TRUE) %>% top_n(n= -show, wt = p.adjust) %>% pull(Description)->terms
      x<-x[x$Description %in% terms,]
      # x <- x[order(x$Count,decreasing=FALSE),]
      # x <- x[order(x$p.adjust,decreasing=TRUE),]
      x$Description <- ifelse(nchar(x$Description)>80,
                              paste(substr(x$Description, 1, 80),"[...]",sep=""),
                              x$Description)
      x$Description <- factor(x$Description, levels = rev(unique(x$Description)))
      #x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
      
      if(colorBy=="regulation"){
        #change color code of dots
        x$DE.regulation<-"no"
        if(nrow(x[grep("_up", x$Condition),])>0){
          x[grep("_up", x$Condition),]$DE.regulation<-"up"
        }
        if(nrow(x[grep("_down", x$Condition),])>0){
          x[grep("_down", x$Condition),]$DE.regulation<-"down"
        }
        
        p<-ggplot(data = x ,aes(x=Condition, y=Description, size=Count , fill=DE.regulation)) +
          geom_point(aes(fill=DE.regulation),pch=21)+theme_classic()+
          theme_linedraw()+
          theme(axis.text.y = element_text(color = "black"),
                axis.text.x = element_text(angle=90, vjust=0.5,hjust=1, color = "black"),
                text = element_text(size = 12))+
          scale_fill_manual(values= c(up = "firebrick3", down = "dodgerblue2"))+
          xlab("")+ylab("")+
          #scale_size_area()
          scale_radius()
        plot(p)
      }
      
      if(colorBy=="padj"){  
        p<- ggplot(data = x, aes(x = Condition, y = Description, color = p.adjust)) +
          geom_point(aes(size = Count)) +
          scale_colour_gradientn(colours=c('red',
                                           'orange',
                                           'darkblue',
                                           'darkblue'),
                                 limits=c(0,1),
                                 values   = c(0,0.05,0.2,0.5,1),
                                 breaks   = c(0.05,0.2,1),
                                 labels = format(c(0.05,0.2,1))) +
          ylab(NULL) +
          theme_bw() +
          theme(text = element_text(size=12), 
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        plot(p)
      }
    }
    
    if(orderBy=="count"){
      x %>% group_by(Condition) %>% arrange(p.adjust, .by_group = TRUE) ->x
      x %>% group_by(Condition) %>% arrange(p.adjust, .by_group = TRUE) %>% top_n(n= show, wt = Count) %>% pull(Description) ->terms
      x<-x[x$Description %in% terms,]
      x$Description <- ifelse(nchar(x$Description)>80,
                              paste(substr(x$Description, 1, 80),"[...]",sep=""),
                              x$Description)
      x$Description <- factor(x$Description, levels = rev(unique(x$Description)))
      #x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
      
      if(colorBy=="regulation"){
        #change color code of dots
        x$DE.regulation<-"no"
        if(nrow(x[grep("_up", x$Condition),])>0){
          x[grep("_up", x$Condition),]$DE.regulation<-"up"
        }
        if(nrow(x[grep("_down", x$Condition),])>0){
          x[grep("_down", x$Condition),]$DE.regulation<-"down"
        }
        
        p<-ggplot(data = x ,aes(x=Condition, y=Description, size=Count , fill=DE.regulation)) +
          geom_point(aes(fill=DE.regulation),pch=21)+theme_classic()+
          theme_linedraw()+
          theme(axis.text.y = element_text(color = "black"),
                axis.text.x = element_text(angle=90, vjust=0.5,hjust=1, color = "black"),
                text = element_text(size = 12))+
          scale_fill_manual(values= c(up = "firebrick3", down = "dodgerblue2"))+
          xlab("")+ylab("")+
          #scale_size_area()
          scale_radius()
        plot(p)
      }
      
      if(colorBy=="padj"){  
        p<- ggplot(data = x, aes(x = Condition, y = Description, color = p.adjust)) +
          geom_point(aes(size = Count)) +
          scale_colour_gradientn(colours=c('red',
                                           'orange',
                                           'darkblue',
                                           'darkblue'),
                                 limits=c(0,1),
                                 values   = c(0,0.05,0.2,0.5,1),
                                 breaks   = c(0.05,0.2,1),
                                 labels = format(c(0.05,0.2,1))) +
          ylab(NULL) +
          theme_bw() +
          theme(text = element_text(size=12), 
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        plot(p)
      }
    }
    
    
  }
}


confusionMatrix <- function (i = NULL, j = NULL)
{
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(i = match(i, ui), j = match(j,
                                                        uj), x = rep(1, length(i)), dims = c(length(ui), length(uj)))
  rownames(m) <- ui
  colnames(m) <- uj
  m
}


# load some functions that may e useful
multiplot<-function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

GSEA <- function(object,
                 DEG = DEG,
                 name = name,
                 GeneSets =c("GO","KEGG","DO","Hallmark","cannonicalPathways","Motifs","ImmunoSignatures"),
                 GOntology = "BP",
                 logfc.threshold = logfc.threshold,
                 pval.use = "p_val_adj", # or p_val for unadjusted
                 pCorrection = "bonferroni", # choose the p-value adjustment method
                 pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                 qvalueCutoff = 0.05 # set the q-value cutoff (FDR corrected)
){
  
  present_genes <- GetAssayData(object = object, slot = 'counts')
  present_genes <- names(which(rowSums(present_genes) > 3))
  present_genes_entrez <- bitr(present_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  up_reg_genes <- row.names(DEG[DEG[,pval.use]<.05&DEG$avg_log2FC>logfc.threshold,])
  entrez_up <- bitr(up_reg_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  down_reg_genes <- row.names(DEG[DEG[,pval.use]<.05&DEG$avg_log2FC<(logfc.threshold*-1),])
  entrez_down <- bitr(down_reg_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  OrgDb = org.Hs.eg.db
  
  results <- list()
  
  
  # GO enrichment
  if("GO" %in% GeneSets){
    #print("Performing GO enrichment")
    results$GOup <- as.data.frame(enrichGO(gene = entrez_up,
                                           universe = present_genes_entrez,
                                           OrgDb = OrgDb,
                                           ont = GOntology,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff  = qvalueCutoff,
                                           readable      = T))
    
    if(nrow(results$GOup)>0){results$GOup$Enrichment <- paste("GO up")}
    if(nrow(results$GOup)>0){results$GOup$comparison <- paste(name)}
    
    if(nrow(results$GOup)==0) results$GOup <- data.frame(Enrichment = "GO up", comparison = paste(name), result = "no enrichment found")
    
    results$GOdown <- as.data.frame(enrichGO(gene = entrez_down,
                                             universe = present_genes_entrez,
                                             OrgDb = OrgDb,
                                             ont = GOntology,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff  = qvalueCutoff,
                                             readable      = T))
    if(nrow(results$GOdown)>0){results$GOdown$Enrichment <- paste("GO down")}
    if(nrow(results$GOdown)>0){results$GOdown$comparison <- paste(name)}
    
    if(nrow(results$GOdown)==0) results$GOdown <- data.frame(Enrichment = "GO down", comparison = paste(name), result = "no enrichment found")
    
    
    
  }
  
  # Reactome enrichment
  if("Reactome" %in% GeneSets){
    #print("Performing Reactome enrichment")
    results$Reactomeup <- as.data.frame(enrichPathway(gene = entrez_up,
                                                      universe = present_genes_entrez,
                                                      OrgDb = OrgDb,
                                                      ont = GOntology,
                                                      pAdjustMethod = pCorrection,
                                                      pvalueCutoff  = pvalueCutoff,
                                                      qvalueCutoff  = qvalueCutoff,
                                                      readable      = T))
    
    if(nrow(results$Reactomeup)>0){results$Reactomeup$Enrichment <- paste("Reactome up",name,sep="")}
    
    results$Reactomedown <- as.data.frame(enrichGO(gene = entrez_down,
                                                   universe = present_genes_entrez,
                                                   OrgDb = OrgDb,
                                                   ont = GOntology,
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff  = qvalueCutoff,
                                                   readable      = T))
    if(nrow(results$Reactomedown)>0){results$Reactomedown$Enrichment <- paste("Reactome down",name,sep="")}
  }
  
  # KEGG enrichment
  if("KEGG" %in% GeneSets){
    #print("Performing KEGG enrichment")
    
    org = "hsa"
    
    results$KEGGup <- as.data.frame(enrichKEGG(gene = entrez_up,
                                               organism = org,
                                               universe = present_genes_entrez,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$KEGGup)>0){results$KEGGup$Enrichment <- paste("KEGG up")}
    if(nrow(results$KEGGup)>0){results$KEGGup$comparison <- paste(name)}
    
    if(nrow(results$KEGGup)==0) results$KEGGup <- data.frame(Enrichment = "KEGG up", comparison = paste(name), result = "no enrichment found")
    
    results$KEGGdown <- as.data.frame(enrichKEGG(gene = entrez_down,
                                                 organism = org,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$KEGGdown)>0){results$KEGGdown$Enrichment <- paste("KEGG down")}
    if(nrow(results$KEGGdown)>0){results$KEGGdown$comparison <- paste(name)}
    
    if(nrow(results$KEGGdown)==0) results$KEGGdown <- data.frame(Enrichment = "KEGG down", comparison = paste(name), result = "no enrichment found")
    
  }
  
  # DO enrichment
  if("DO" %in% GeneSets){
    #print("Performing Disease Ontology enrichment")
    
    results$DOup <- as.data.frame(enrichDO(gene = entrez_up,
                                           universe = present_genes_entrez,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff = qvalueCutoff,
                                           minGSSize     = 5,
                                           maxGSSize     = 500,
                                           readable=TRUE))
    if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in ",name,sep="")}
    if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in ",name,sep="")}
    
    results$DOdown <- as.data.frame(enrichDO(gene = entrez_down,
                                             universe = present_genes_entrez,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff = qvalueCutoff,
                                             minGSSize     = 5,
                                             maxGSSize     = 500,
                                             readable=TRUE))
    if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in ",name,sep="")}
    if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in ",name,sep="")}
    
  }
  
  # Hallmark enrichment
  if("Hallmark" %in% GeneSets){
    #print("Performing Hallmark enrichment")
    
    results$HALLMARKup <- as.data.frame(enricher(entrez_up,
                                                 TERM2GENE=hallmark_genes,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$HALLMARKup)>0){results$HALLMARKup$Enrichment <- paste("HALLMARK up")}
    if(nrow(results$HALLMARKup)>0){results$HALLMARKup$comparison <- paste(name)}
    
    if(nrow(results$HALLMARKup)==0) results$HALLMARKup <- data.frame(Enrichment = "HALLMARK up", comparison = paste(name), result = "no enrichment found")
    
    
    results$HALLMARKdown <- as.data.frame(enricher(entrez_down,
                                                   TERM2GENE=gmtfile_hallmarks,
                                                   universe = present_genes_entrez,
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
    if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste("HALLMARK down")}
    if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste(name)}
    
    if(nrow(results$HALLMARKdown)==0) results$HALLMARKdown <- data.frame(Enrichment = "HALLMARK down", comparison = paste(name), result = "no enrichment found")
    
    
  }
  
  # Cannonical Pathway enrichment
  if("cannonicalPathways" %in% GeneSets){
    #print("Performing Cannonical Pathway (C2) enrichment")
    
    results$cannonicalPathwaysup <- as.data.frame(enricher(entrez_up,
                                                           TERM2GENE=cannonicalPathway_genes,
                                                           universe = present_genes_entrez,
                                                           pAdjustMethod = pCorrection,
                                                           pvalueCutoff  = pvalueCutoff,
                                                           qvalueCutoff = qvalueCutoff))
    if(nrow(results$cannonicalPathwaysup)>0){results$cannonicalPathwaysup$Enrichment <- paste("Cannonical pathway enrichment for genes upregulated in ",name_cluster,sep="")}
    
    results$cannonicalPathwaysdown <- as.data.frame(enricher(entrez_down,
                                                             TERM2GENE=cannonicalPathway_genes,
                                                             universe = present_genes_entrez,
                                                             pAdjustMethod = pCorrection,
                                                             pvalueCutoff  = pvalueCutoff,
                                                             qvalueCutoff = qvalueCutoff))
    if(nrow(results$cannonicalPathwaysdown)>0){results$cannonicalPathwaysdown$Enrichment <- paste("Cannonical pathway enrichment for genes downregulated in ",name_cluster,sep="")}
  }
  
  # Motif enrichment
  if("Motifs" %in% GeneSets){
    #print("Performing Motif enrichment")
    
    results$Motifup <- as.data.frame(enricher(entrez_up,
                                              TERM2GENE=motifs,
                                              universe = present_genes_entrez,
                                              pAdjustMethod = pCorrection,
                                              pvalueCutoff  = pvalueCutoff,
                                              qvalueCutoff = qvalueCutoff))
    if(nrow(results$Motifup)>0){results$Motifup$Enrichment <- paste("TF binding motif enrichment for genes upregulated in cluster ",name,sep="")}
    
    results$Motifdown <- as.data.frame(enricher(entrez_down,
                                                TERM2GENE=motifs,
                                                universe = present_genes_entrez,
                                                pAdjustMethod = pCorrection,
                                                pvalueCutoff  = pvalueCutoff,
                                                qvalueCutoff = qvalueCutoff))
    if(nrow(results$Motifdown)>0){results$Motifdown$Enrichment <- paste("TF binding motif enrichment for genes downregulated in cluster",name,sep="")}
  }
  
  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    #print("Performing immunesignature enrichment")
    
    results$ImmSigup <- as.data.frame(enricher(entrez_up,
                                               TERM2GENE=gmtfile_immunosignatures,
                                               universe = present_genes_entrez,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$ImmSigup)>0){results$ImmSigup$Enrichment <- paste("Immunosignature enrichment for genes upregulated in cluster ",name,sep="")}
    
    results$ImmSigdown <- as.data.frame(enricher(entrez_down,
                                                 TERM2GENE=gmtfile_immunosignatures,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$ImmSigdown)>0){results$ImmSigdown$Enrichment <- paste("Immunosignature enrichment for genes downregulated in cluster ",name,sep="")}
  }
  return(results)
}



GSEA2 <- function(object,
                  up_reg_genes = up_reg_genes,
                  down_reg_genes = down_reg_genes,
                  name = name,
                  GeneSets =c("GO","KEGG","DO","Hallmark","cannonicalPathways","Motifs","ImmunoSignatures"),
                  GOntology = "BP",
                  logfc.threshold = logfc.threshold,
                  pval.use = "p_val_adj", # or p_val for unadjusted
                  pCorrection = "bonferroni", # choose the p-value adjustment method
                  pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                  qvalueCutoff = 0.05 # set the q-value cutoff (FDR corrected)
){
  
  present_genes <- GetAssayData(object = object, slot = 'counts')
  present_genes <- names(which(rowSums(present_genes) > 3))
  present_genes_entrez <- bitr(present_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  up_reg_genes <- up_reg_genes
  entrez_up <- bitr(up_reg_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  down_reg_genes <- down_reg_genes
  entrez_down <- bitr(down_reg_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  OrgDb = org.Hs.eg.db
  
  results <- list()
  
  
  # GO enrichment
  if("GO" %in% GeneSets){
    #print("Performing GO enrichment")
    results$GOup <- as.data.frame(enrichGO(gene = entrez_up,
                                           universe = present_genes_entrez,
                                           OrgDb = OrgDb,
                                           ont = GOntology,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff  = qvalueCutoff,
                                           readable      = T))
    
    if(nrow(results$GOup)>0){results$GOup$Enrichment <- paste("GO up")}
    if(nrow(results$GOup)>0){results$GOup$comparison <- paste(name)}
    
    if(nrow(results$GOup)==0) results$GOup <- data.frame(Enrichment = "GO up", comparison = paste(name), result = "no enrichment found")
    
    results$GOdown <- as.data.frame(enrichGO(gene = entrez_down,
                                             universe = present_genes_entrez,
                                             OrgDb = OrgDb,
                                             ont = GOntology,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff  = qvalueCutoff,
                                             readable      = T))
    if(nrow(results$GOdown)>0){results$GOdown$Enrichment <- paste("GO down")}
    if(nrow(results$GOdown)>0){results$GOdown$comparison <- paste(name)}
    
    if(nrow(results$GOdown)==0) results$GOdown <- data.frame(Enrichment = "GO down", comparison = paste(name), result = "no enrichment found")
    
    
    
  }
  
  # Reactome enrichment
  if("Reactome" %in% GeneSets){
    #print("Performing Reactome enrichment")
    results$Reactomeup <- as.data.frame(enrichPathway(gene = entrez_up,
                                                      universe = present_genes_entrez,
                                                      OrgDb = OrgDb,
                                                      ont = GOntology,
                                                      pAdjustMethod = pCorrection,
                                                      pvalueCutoff  = pvalueCutoff,
                                                      qvalueCutoff  = qvalueCutoff,
                                                      readable      = T))
    
    if(nrow(results$Reactomeup)>0){results$Reactomeup$Enrichment <- paste("Reactome up",name,sep="")}
    
    results$Reactomedown <- as.data.frame(enrichGO(gene = entrez_down,
                                                   universe = present_genes_entrez,
                                                   OrgDb = OrgDb,
                                                   ont = GOntology,
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff  = qvalueCutoff,
                                                   readable      = T))
    if(nrow(results$Reactomedown)>0){results$Reactomedown$Enrichment <- paste("Reactome down",name,sep="")}
  }
  
  # KEGG enrichment
  if("KEGG" %in% GeneSets){
    #print("Performing KEGG enrichment")
    
    org = "hsa"
    
    results$KEGGup <- as.data.frame(enrichKEGG(gene = entrez_up,
                                               organism = org,
                                               universe = present_genes_entrez,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$KEGGup)>0){results$KEGGup$Enrichment <- paste("KEGG up")}
    if(nrow(results$KEGGup)>0){results$KEGGup$comparison <- paste(name)}
    
    if(nrow(results$KEGGup)==0) results$KEGGup <- data.frame(Enrichment = "KEGG up", comparison = paste(name), result = "no enrichment found")
    
    results$KEGGdown <- as.data.frame(enrichKEGG(gene = entrez_down,
                                                 organism = org,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$KEGGdown)>0){results$KEGGdown$Enrichment <- paste("KEGG down")}
    if(nrow(results$KEGGdown)>0){results$KEGGdown$comparison <- paste(name)}
    
    if(nrow(results$KEGGdown)==0) results$KEGGdown <- data.frame(Enrichment = "KEGG down", comparison = paste(name), result = "no enrichment found")
    
  }
  
  # DO enrichment
  if("DO" %in% GeneSets){
    #print("Performing Disease Ontology enrichment")
    
    results$DOup <- as.data.frame(enrichDO(gene = entrez_up,
                                           universe = present_genes_entrez,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff = qvalueCutoff,
                                           minGSSize     = 5,
                                           maxGSSize     = 500,
                                           readable=TRUE))
    if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in ",name,sep="")}
    if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in ",name,sep="")}
    
    results$DOdown <- as.data.frame(enrichDO(gene = entrez_down,
                                             universe = present_genes_entrez,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff = qvalueCutoff,
                                             minGSSize     = 5,
                                             maxGSSize     = 500,
                                             readable=TRUE))
    if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in ",name,sep="")}
    if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in ",name,sep="")}
    
  }
  
  # Hallmark enrichment
  if("Hallmark" %in% GeneSets){
    #print("Performing Hallmark enrichment")
    
    results$HALLMARKup <- as.data.frame(enricher(entrez_up,
                                                 TERM2GENE=gmtfile_hallmarks,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$HALLMARKup)>0){results$HALLMARKup$Enrichment <- paste("HALLMARK up")}
    if(nrow(results$HALLMARKup)>0){results$HALLMARKup$comparison <- paste(name)}
    
    if(nrow(results$HALLMARKup)==0) results$HALLMARKup <- data.frame(Enrichment = "HALLMARK up", comparison = paste(name), result = "no enrichment found")
    
    
    results$HALLMARKdown <- as.data.frame(enricher(entrez_down,
                                                   TERM2GENE=gmtfile_hallmarks,
                                                   universe = present_genes_entrez,
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
    if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste("HALLMARK down")}
    if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste(name)}
    
    if(nrow(results$HALLMARKdown)==0) results$HALLMARKdown <- data.frame(Enrichment = "HALLMARK down", comparison = paste(name), result = "no enrichment found")
    
    
  }
  
  # Cannonical Pathway enrichment
  if("cannonicalPathways" %in% GeneSets){
    #print("Performing Cannonical Pathway (C2) enrichment")
    
    results$cannonicalPathwaysup <- as.data.frame(enricher(entrez_up,
                                                           TERM2GENE=cannonicalPathway_genes,
                                                           universe = present_genes_entrez,
                                                           pAdjustMethod = pCorrection,
                                                           pvalueCutoff  = pvalueCutoff,
                                                           qvalueCutoff = qvalueCutoff))
    if(nrow(results$cannonicalPathwaysup)>0){results$cannonicalPathwaysup$Enrichment <- paste("Cannonical pathway enrichment for genes upregulated in ",name_cluster,sep="")}
    
    results$cannonicalPathwaysdown <- as.data.frame(enricher(entrez_down,
                                                             TERM2GENE=cannonicalPathway_genes,
                                                             universe = present_genes_entrez,
                                                             pAdjustMethod = pCorrection,
                                                             pvalueCutoff  = pvalueCutoff,
                                                             qvalueCutoff = qvalueCutoff))
    if(nrow(results$cannonicalPathwaysdown)>0){results$cannonicalPathwaysdown$Enrichment <- paste("Cannonical pathway enrichment for genes downregulated in ",name_cluster,sep="")}
  }
  
  # Motif enrichment
  if("Motifs" %in% GeneSets){
    #print("Performing Motif enrichment")
    
    results$Motifup <- as.data.frame(enricher(entrez_up,
                                              TERM2GENE=motifs,
                                              universe = present_genes_entrez,
                                              pAdjustMethod = pCorrection,
                                              pvalueCutoff  = pvalueCutoff,
                                              qvalueCutoff = qvalueCutoff))
    if(nrow(results$Motifup)>0){results$Motifup$Enrichment <- paste("TF binding motif enrichment for genes upregulated in cluster ",name,sep="")}
    
    results$Motifdown <- as.data.frame(enricher(entrez_down,
                                                TERM2GENE=motifs,
                                                universe = present_genes_entrez,
                                                pAdjustMethod = pCorrection,
                                                pvalueCutoff  = pvalueCutoff,
                                                qvalueCutoff = qvalueCutoff))
    if(nrow(results$Motifdown)>0){results$Motifdown$Enrichment <- paste("TF binding motif enrichment for genes downregulated in cluster",name,sep="")}
  }
  
  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    #print("Performing immunesignature enrichment")
    
    results$ImmSigup <- as.data.frame(enricher(entrez_up,
                                               TERM2GENE=gmtfile_immunosignatures,
                                               universe = present_genes_entrez,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$ImmSigup)>0){results$ImmSigup$Enrichment <- paste("Immunosignature enrichment for genes upregulated in cluster ",name,sep="")}
    
    results$ImmSigdown <- as.data.frame(enricher(entrez_down,
                                                 TERM2GENE=gmtfile_immunosignatures,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$ImmSigdown)>0){results$ImmSigdown$Enrichment <- paste("Immunosignature enrichment for genes downregulated in cluster ",name,sep="")}
  }
  return(results)
}

