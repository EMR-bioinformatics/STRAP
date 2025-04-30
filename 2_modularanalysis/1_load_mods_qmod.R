######### loading modules

loadmodsforqusage <- function(load='LI'){ # tmod, wgcna
  if (load == 'LI'){
    library(tmod)
    #data(tmod)
    tmod <- readRDS('~/../../media/gcpeac/Anna/functions/modular_analysis/LI_modules')
    c2.indices <- tmod$MODULES2GENES
    c2.indices <- c2.indices[1:334]
    names(c2.indices) <- paste(names(c2.indices),tmod$MODULES$Title[1:334])
    names(c2.indices) <- gsub('LI.','',names(c2.indices))
    
  }else if (load == 'LI_reduced'){
    library(tmod)
    #data(tmod)
    tmod <- readRDS('~/../../media/gcpeac/Anna/functions/modular_analysis/LI_modules')
    c2.indices <- tmod$MODULES2GENES
    c2.indices <- c2.indices[1:334]
    names(c2.indices) <- paste(names(c2.indices),tmod$MODULES$Title[1:334])
    names(c2.indices) <- gsub('LI.','',names(c2.indices))
    TBA <- grep('TBA', names(c2.indices))
    c2.indices <- c2.indices[-TBA]
    
  }else if (load == 'wgcna'){
    # wgcna
    df <- read.csv('~/../../media/gcpeac/Chris_Desktop/PEAC/8_modularframework/synovium_baseline_modules_annotatedv11.csv')
    df <- df[,-1]
    final_list <- list()
    modules <- unique(df$Module)
    for (module in modules){
      tempdf <- subset(df, Module == module)
      new_list <- list(as.character(tempdf$Symbol))
      names(new_list)[1]<-module
      final_list <- c(final_list, new_list)
    }
    c2.indices <- final_list
  } else if (load == 'reduced_wgcna'){
    df <- read.csv('~/../../media/gcpeac/Anna/functions/modular_analysis/reduced_wgcna_modules_V1.csv')
    final_list <- list()
    modules <- unique(df$Module)
    for (module in modules){
      tempdf <- subset(df, Module == module)
      new_list <- list(as.character(tempdf$Symbol))
      names(new_list)[1]<-module
      final_list <- c(final_list, new_list)
    }
    c2.indices <- final_list
  }else if (load == 'singlecell'){
    # stephenson single cell
    df <- read.csv('~/../../media/gcpeac/Chris_Desktop/PEAC/8_modularframework/singlecelldata/SCreadyforannotation.csv')
    df <- subset(df, df$avg_logFC > 0)
    df <- df[,c('cluster','gene')]
    colnames(df) <- c('Module','Symbol')
    final_list <- list()
    modules <- unique(df$Module)
    for (module in modules){
      tempdf <- subset(df, Module == module)
      new_list <- list(as.character(tempdf$Symbol))
      names(new_list)[1]<-module
      final_list <- c(final_list, new_list)
    }
    c2.indices <- final_list
  }else if (load == 'KEGG'){
    df <- clusterProfiler::read.gmt('~/../../media/gcpeac/Chris_Desktop/Genesets/c2.cp.kegg.v6.1.symbols.gmt')
    colnames(df) <- c('Module','Symbol')
    final_list <- list()
    modules <- unique(df$Module)
    for (module in modules){
      tempdf <- subset(df, Module == module)
      new_list <- list(as.character(tempdf$Symbol))
      names(new_list)[1]<-module
      final_list <- c(final_list, new_list)
    }
    c2.indices <- final_list
  }else if (load == 'REACTOME'){
    df <- clusterProfiler::read.gmt("~/../../media/gcpeac/Chris_Desktop/Genesets/c2.cp.reactome.v6.1.symbols.gmt")
    colnames(df) <- c('Module','Symbol')
    final_list <- list()
    modules <- unique(df$Module)
    for (module in modules){
      tempdf <- subset(df, Module == module)
      new_list <- list(as.character(tempdf$Symbol))
      names(new_list)[1]<-module
      final_list <- c(final_list, new_list)
    }
    c2.indices <- final_list
  }
  return(c2.indices)
}
