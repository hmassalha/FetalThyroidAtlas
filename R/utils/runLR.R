runLR = function(REF.srat,ref_annot='annot',srat,LR_level='both',srat_annot='',minGeneMatch=0.99,model_fp,tissue='',scLR_TGTtype='',scLR_REFtype='',outDir=NULL,plot_prefix=NULL,out_prefix=NULL,seed=2397,maxCells=2000,skipIfExists=T,...){
  #source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/logisticRegression.R")
  
  ### Check outDir ##
  if(!is.null(outDir) & !dir.exists(outDir)){
    message('Generating output directory')
    dir.create(outDir,recursive = T)
  }
  plotPath = outDir
  REF.srat@meta.data$annot = as.character(REF.srat@meta.data[[ref_annot]])
  
  #### Logistic Regression ####
  # Get raw count matrix
  REF_1 = c()
  REF_2 = c()
  for(cluster in unique(REF.srat$annot)){
    tmp = REF.srat@meta.data[REF.srat@meta.data$annot == cluster,]
    set.seed(seed)
    cells1 = sample(rownames(tmp),size = nrow(tmp)*0.7)
    REF_1 = c(REF_1,cells1)
    set.seed(seed)
    cells2 = rownames(tmp)[!rownames(tmp) %in% cells1]
    REF_2 = c(REF_2,cells2)
  }
  
  REF_1_mtx = REF.srat@assays$RNA@counts[,colnames(REF.srat@assays$RNA@counts) %in% REF_1]
  REF_2_mtx = REF.srat@assays$RNA@counts[,colnames(REF.srat@assays$RNA@counts) %in% REF_2]
  m = match(colnames(REF_2_mtx), rownames(REF.srat@meta.data))
  sum(is.na(m))
  colnames(REF_2_mtx) = c(paste0('ref_',colnames(REF_2_mtx)))
  # Generate test matrix (made up of AK mtx and REF_2 mtx for controls)
  srat@meta.data$cellID = rownames(srat@meta.data)
  mtx = srat@assays$RNA@counts
  
  ## Perform LR Training and Prediction ##
  m=match(colnames(REF_1_mtx), rownames(REF.srat@meta.data))
  sum(is.na(m))
  refClasses = REF.srat@meta.data$annot[m]
  
  ## Remove refClasses with < 10 cells in training dataset
  w = which(table(refClasses) < 10)
  if(length(w) > 0){
    # remove these classes from REF_1_mtx and REF_2_mtx
    REF_1_mtx = REF_1_mtx[,colnames(REF_1_mtx) %in% rownames(REF.srat@meta.data[!REF.srat$annot %in% names(w),])]
    REF_2_mtx = REF_2_mtx[,gsub('ref_','',colnames(REF_2_mtx)) %in% rownames(REF.srat@meta.data[!REF.srat$annot %in% names(w),])]
    # re-compute refClasses
    m=match(colnames(REF_1_mtx), rownames(REF.srat@meta.data))
    sum(is.na(m))
    refClasses = REF.srat@meta.data$annot[m]
  }
  
  
  
  
  if((!is.null(model_fp)) & file.exists(model_fp) & skipIfExists==T){
    fit_REF_1 = readRDS(model_fp)
  }else if((is.null(model_fp)) | !file.exists(model_fp) | !skipIfExists){
    set.seed(seed) 
    fit_REF_1 = trainModel(REF_1_mtx, refClasses,maxCells = maxCells,...)
    if(!is.null(model_fp)){
      saveRDS(fit_REF_1,model_fp)
    }  
  }
  
  
  LR_outputs = list()
  if(LR_level == 'cluster'| LR_level == 'both'){
    # Prediction on test dataset (AK)
    m = match(gsub('ref_','',colnames(REF_2_mtx)), rownames(REF.srat@meta.data))
    sum(is.na(m))
    REF_clusterIDs = c(paste0('ref_',REF.srat@meta.data$annot[m]))
    
    if(srat_annot == ''){
      message('No annotation column for test_srat object provided - use seurat_clusters by default')
      srat_annot = 'seurat_clusters'
    }
    clusterIDs = as.character(srat@meta.data[[srat_annot]])
    
    set.seed(seed)
    tgt_logit_clusterIDs = predictSimilarity(fit_REF_1, mtx,logits = T,classes = clusterIDs,minGeneMatch = minGeneMatch)
    REF_logit_clusterIDs = predictSimilarity(fit_REF_1, REF_2_mtx,logits = T,classes = REF_clusterIDs)
    
    m = match(colnames(REF_logit_clusterIDs),colnames(tgt_logit_clusterIDs))
    sum(is.na(m))
    logit_clusterIDs = rbind(REF_logit_clusterIDs,tgt_logit_clusterIDs[,m])
    LR_outputs = append(LR_outputs,list(clusterLR = logit_clusterIDs))
  }
  
  if(LR_level == 'sc'| LR_level == 'both'){
    #===========================================
    # LR prediction at single cell levels     ##  
    set.seed(seed)
    tgt_logit_scLevel = predictSimilarity(fit_REF_1, mtx,logits = T,minGeneMatch = minGeneMatch)
    REF_logit_scLevel = predictSimilarity(fit_REF_1, REF_2_mtx,logits = T)
    REF_logit_scLevel = as.data.frame(REF_logit_scLevel)
    #REF_logit_scLevel$celltype = big.srat@meta.data[gsub('ref_','',rownames(REF_logit_scLevel)),]$finalAnn
    ### Now, what I want to do is to use the information from training REF_on_REF (validation) to set a threshold for what's supposed to be confident cell labels
    ### The threshold for each celltype should be at the point where 90% of the cells within that cluster is a PASS
    # Find the threshold
    # Slice the matrix by celltype cluster ID
    REF_logit_scLevel_byCT = split(REF_logit_scLevel,REF.srat@meta.data[gsub('ref_','',rownames(REF_logit_scLevel)),ref_annot])
    # calculate_10th_percentile = function(logit){
    #   # convert back to probs
    #   p = ((exp(-logit[!is.na(logit)])+1)^-1)
    #   t = quantile(p, probs = c(.1))
    #   return(t)
    # }
    
    percentile_10th_cutoff_byCT = lapply(seq_along(REF_logit_scLevel_byCT), function(i){
      if(!names(REF_logit_scLevel_byCT)[i] %in% (colnames(REF_logit_scLevel_byCT[[i]]))){
        return('NA')
      }else if(sum(is.na(REF_logit_scLevel_byCT[[i]][,names(REF_logit_scLevel_byCT)[i]]))>1){
        return('NA')
      }else{
        return(REF_logit_scLevel_byCT[[i]][,names(REF_logit_scLevel_byCT)[i]] %>%  quantile(probs=c(.1)))  
      }
    })
      
    percentile_5th_cutoff_byCT = lapply(seq_along(REF_logit_scLevel_byCT), function(i){
      if(!names(REF_logit_scLevel_byCT)[i] %in% (colnames(REF_logit_scLevel_byCT[[i]]))){
        return('NA')
      }else if(sum(is.na(REF_logit_scLevel_byCT[[i]][,names(REF_logit_scLevel_byCT)[i]]))>1){
        return('NA')
      }else{
        return(REF_logit_scLevel_byCT[[i]][,names(REF_logit_scLevel_byCT)[i]] %>%  quantile(probs=c(.05)))  
      }
    })
    
    names(percentile_10th_cutoff_byCT) = names(REF_logit_scLevel_byCT)
    names(percentile_5th_cutoff_byCT) = names(REF_logit_scLevel_byCT)
    
    # Now that we've got the threshold for each cell type, let's go back to the single cell target LR output and label cells
    # For each cell, choose the cell label with the highest match. 
    # Look at the logit score for this highest match, if it is >= the corresponding 10th percentile threshold for that cell type --> confident match
    # Look at the logit score for this highest match, if it is >= the corresponding 20th percentile threshold for that cell type --> low confident match
    # else - unknown
    
    x = apply(tgt_logit_scLevel,MARGIN = 1, FUN = function(x){
      x = x[!is.na(x)]
      ifelse(sum(x == max(x)) > 1,# ie. equivalent cells
             paste(c('ambiguous',names(x[x==max(x)])),sep = ':'),
             ifelse(max(x) >= percentile_10th_cutoff_byCT[[names(x)[x==max(x)]]],
                    names(x[x==max(x)]),
                    ifelse(max(x) >= percentile_5th_cutoff_byCT[[names(x[x==max(x)])]],
                           paste0('lowConf_',names(x[x==max(x)])),'unknown')))})
    
    if(length(x) != nrow(srat@meta.data)){
      stop(sprintf('Something went wrong with LR...'))
    }
    srat@meta.data$celltype_LR = x[rownames(srat@meta.data)]
    
    
    # TMP solution: too many REF cells --> subset to keep 10% for plotting only
    #set.seed(seed)
    #cells1 = sample(rownames(REF_logit_scLevel),size = nrow(REF_logit_scLevel)*0.1)
    #REF_logit_scLevel = REF_logit_scLevel[cells1,]
    
    m = match(colnames(REF_logit_scLevel),colnames(tgt_logit_scLevel))
    sum(is.na(m))
    logit_sc = rbind(REF_logit_scLevel,tgt_logit_scLevel[,m])
    #LR_outputs = append(LR_outputs,list(scLR = logit_sc))
    
    # tmp sol again
    LR_outputs = append(LR_outputs,list(scLR_all = logit_sc,
                                        scLR_tgt = tgt_logit_scLevel,
                                        scLR_ref = REF_logit_scLevel))
  }
  
  # Save output LR matrices 
  if(!is.null(out_prefix)){
    saveRDS(LR_outputs,file.path(outDir,paste0(out_prefix,'raw_LR_outputs.RDS')))
  }
  
  ## Add results to srat object and return ##
  for(i in 1:length(LR_outputs)){
    output_df = as.data.frame(LR_outputs[[i]])
    write.csv(output_df,file.path(outDir,paste0(out_prefix,names(LR_outputs[i]),'.csv')))
    if(names(LR_outputs[i]) == 'clusterLR'){
      # Add automated cluster prioritization
      
      # convert from logit back to probability between 0-1
      output_df$LR_celltype = apply(output_df,MARGIN = 1,FUN = function(x){ifelse(max(((exp(-x[!is.na(x)])+1)^-1))<0.7,'unknown',
                                                                                  ifelse(max((exp(-x[!is.na(x)])+1)^-1)<0.9,paste0('lowScore_',colnames(output_df)[which(x==max(x[!is.na(x)]))]),(colnames(output_df)[which(x==max(x[!is.na(x)]))])))})
      output_df$highest_LR_LLH = apply(output_df[,-ncol(output_df)],MARGIN = 1,FUN = function(x){max(x[!is.na(x)])})
      
      output_df$clusterID = rownames(output_df)
      
      m = match(as.character(srat@meta.data[[srat_annot]]), output_df$clusterID)
      sum(is.na(m))
      srat$LR_celltype.ClusterIDs = output_df$LR_celltype[m]
      srat$highest_LR_LLH.ClusterIDs = output_df$highest_LR_LLH[m]
    }else if(names(LR_outputs[i]) == 'scLR_tgt'){#else if(grepl('scLR',names(LR_outputs[i]))){
    
      # Convert the number back to probability
      output_df$LR_celltype.scLevel = apply(output_df,MARGIN = 1,FUN = function(x){ifelse(max(((exp(-x[!is.na(x)])+1)^-1))<0.7,'unknown',
                                                                                                        ifelse(max((exp(-x[!is.na(x)])+1)^-1)<0.9,paste0('lowScore_',colnames(output_df)[which(x==max(x[!is.na(x)]))]),(colnames(output_df)[which(x==max(x[!is.na(x)]))])))})
      output_df$highest_LR_LLH.scLevel = apply(output_df[,-ncol(output_df)],MARGIN = 1,FUN = function(x){max(x[!is.na(x)])})
      if('Row.names' %in% colnames(srat@meta.data)){
        srat@meta.data = srat@meta.data[,-which(colnames(srat@meta.data) == 'Row.names')]  
      }
      
      #m = match(rownames(srat@meta.data), rownames(output_df))
      #sum(is.na(m))
      #srat@meta.data = merge(srat@meta.data,output_df[m,],by=0,all = T)
      #rownames(srat@meta.data) = srat@meta.data$cellID
    }
  }  
    
  # Plot HeatMap
  if(!is.null(plot_prefix)){
    for(i in 1:length(LR_outputs)){
      message(paste0('Set #',i))
      output = LR_outputs[[i]]
      if(names(LR_outputs[i]) == 'clusterLR'){
        type = ifelse(grepl('ref_',rownames(output)),'REF','TGT')
        show_row_names = T
      }else if(names(LR_outputs)[i] == 'scLR_all'){
        next
      }else if(grepl('scLR',names(LR_outputs[i]))){
        type = ifelse(grepl('ref_',rownames(output)),'REF','TGT')
        show_row_names = F
        
        
        #type = rep('cells',nrow(output))
        #t=match(rownames(output),rownames(srat@meta.data))
        #output=output[!is.na(t),]
        #t=match(rownames(output),rownames(srat@meta.data))
        #type[!is.na(t)] = ifelse(rep((scLR_TGTtype == ''),sum(!is.na(t))),'TGT',
        #                         as.character(srat@meta.data[[scLR_TGTtype]])[t[!is.na(t)]])
        
        #r=match(rownames(output),rownames(REF.srat@meta.data))
        #type[!is.na(r)] = ifelse(rep(scLR_REFtype == '',sum(!is.na(r))),'REF',
        #                         as.character(REF.srat@meta.data[[scLR_REFtype]])[r[!is.na(r)]])
        
        
        if(sum(is.na(type))>0){
          stop()
        }
        
      }
      
      
      if(tissue == 'kidney'){
        column_order = c('End', 'MSC', "ICb","ICa","RVCSB","NPC","ErPrT","Pod","SSBm.d","DTLH","SSBpod","UBCD","CnT","SSBpr") 
      }else if(tissue == 'adrenal'){
        column_order = c('SCPs','Bridge','Sympathoblastic','Chromaffin','Mesenchyme - early','Mesenchyme','Cortex','Endothelium','Leukocytes','Erythroblasts','Other')
      }else if(tissue == 'Adrenal.Kidney'){
        column_order = c('SCPs','Bridge','Sympathoblastic','Chromaffin','Mesenchyme - early','Mesenchyme','Cortex','Endothelium','Leukocytes','Erythroblasts','Other','End', 'MSC', "ICb","ICa","RVCSB","NPC","ErPrT","Pod","SSBm.d","DTLH","SSBpod","UBCD","CnT","SSBpr")
      }else{
        column_order = colnames(output)
      }
      
      pdf(file.path(outDir,paste0(plot_prefix,names(LR_outputs[i]),'.pdf')),width = 20,height = 30)
      
      hm = similarityHeatmap(output,
                             #column_order = column_order,
                             row_title_rot = 0,
                        row_title_gp = gpar(fontsize=10),row_names_gp = gpar(fontsize=10),row_names_max_width = unit(6,'cm'),
                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
      draw(hm)
      dev.off()
      
      
      
      #type=type[grepl('^SETBP1_',rownames(output))]
      #pdf(file.path(outDir,paste0(plot_prefix,names(LR_outputs[i]),'_CT.pdf')),width = 30,height = 15)
      #output = output[grepl('^SETBP1_',rownames(output)),]
      
      #similarityHeatmap(output,column_order = column_order,row_title_rot = 0,
      #                  row_title_gp = gpar(fontsize=20),row_names_gp = gpar(fontsize=20),row_names_max_width = unit(80,'cm'),
      #                  column_names_gp = gpar(fontsize=20),column_names_max_height = unit(80,'cm'),
      #                  split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
      
      
      #dev.off()
    
    }
  }
  
  
  
  return(list(srat,LR_outputs))
  
}
