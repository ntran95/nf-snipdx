### set s4 object for snipdx

# set multiple types options
setClassUnion("SegmentsOrNull",c("data.frame","NULL"))
setClassUnion("Seg_VariantSelectedOrNull",c("data.frame","NULL"))
setClassUnion("Seg_CallsOrNull",c("data.frame","NULL"))
setClassUnion("Seg_dfInventoryOrNull",c("data.frame","NULL"))

setClass("snipdx", representation(dfIntegInventory = "data.frame",
                                  Stats = "data.frame",
                                  Variant = "data.frame",
                                  Facets = "data.frame",
                                  #CNVLong = "data.frame",
                                  Calls = "data.frame",
                                  Segments = "SegmentsOrNull",
                                  Seg_VariantSelected = "Seg_VariantSelectedOrNull",
                                  Seg_Calls = "Seg_CallsOrNull",
                                  Seg_dfInventory = "Seg_dfInventoryOrNull",
                                  ngs = "list",
                                  trial = "character",
                                  date = "Date"))



#### TODO:
createSNiPDxObj <- function(Stats,
                            Variant,
                            Facets,
                            Inventory,
                            OrigSegments,
                            trial){
  
  #### ==== GENE LVL ==== ####
  
  # modify Stats Quality params (updated 2022-07-29)
  Stats <- Stats %>%
    mutate(Facets.Quality = case_when(`Mean Library Insert Size` < 100 |
                                        `Percentage Duplicate Reads` > 50 |
                                        `Mean Coverage`  < 1000 ~ "Poor",
                                      TRUE ~ as.character("Good")))
  
  
  Variant$VAF <- Variant$FAF/100
  
  allIDs <- data.frame(`Sample Sequencing Name` =  unique(c(Inventory$`Sample Sequencing Name`)),
                       check.names = F)
  
  StatsSelected <-  Stats %>%
    select(`Sample Sequencing Name`,
           `Percentage Duplicate Reads`,
           `Mean Library Insert Size`,
           `Mean Coverage`,
           Facets.Quality) %>%
    inner_join(allIDs)
  
  FacetsSelected <- Facets %>%
    select(`Sample Sequencing Name` = full_sample_id,
           facets.ploidy,
           facets.purity,
           Gene,
           facets.tcn.em.median,
           facets.tcn.em.min,
           facets.tcn.em.max,
           facets.lcn.em.median,
           facets.lcn.em.min,
           facets.loh.pier) %>%
    inner_join(allIDs)
  
  # CNVLong <- CNV %>%
  #   inner_join(allIDs) %>%
  #   select(-F.Ploidy, -F.Purity) %>%
  #   melt(.,id.vars=c(1:2),variable.name = "Gene",value.name = "Call") %>%
  #   inner_join(allIDs)
  
  VariantSelected <- Variant %>%
    select(`Sample Sequencing Name`,
           Gene,
           `hGVS Protein`,
           `hGVS cDNA`,
           VAF,
           `Variant Effect`,
           `Variant Type`,
           `Rep. Class.`,
           `Germline Class.`,
           Rs.Clinical.Significance,
           Chr.,
           Pos.,
           Ref.,
           Alt.,
           FDP) %>%
    inner_join(allIDs)
  
  dfIntegSNiPDx <- allIDs %>%
    left_join(StatsSelected) %>%
    left_join(FacetsSelected, by = c("Sample Sequencing Name")) %>%
    #left_join(CNVLong,by=c("Sample Sequencing Name", "Gene")) %>%
    left_join(VariantSelected)
  
  # ### Calculate Expected biallelic VAF and confidenc eintervals ---OBSOLETE
  #
  # dfIntegSNiPDx <- dfIntegSNiPDx %>%
  #   mutate(SomaticBiallelicVAF = case_when(is.na(facets.purity) ~ 0,
  #                                          TRUE ~ CalcSBVAF(facets.purity,facets.tcn.em.median))) %>%
  #   mutate(GermlineBiallelicVAF = case_when(is.na(facets.purity) ~ 0,
  #                                           TRUE ~ CalcGBVAF(facets.purity,facets.tcn.em.median)))
  
  #### SGZ ####
  # iterate through multiple mutational changes within the same sample
  mut_change_df <- dfIntegSNiPDx %>%
    filter(!is.na(`hGVS Protein`)) %>%
    group_by(`Sample Sequencing Name`) %>%
    distinct(`hGVS Protein`, .keep_all = T)
  
  dfSGZList <- list()
  
  #iterate thru each sample and alteration level, subset, and perform SGZ
  for(currMutSample in seq(nrow(mut_change_df))){
    print(mut_change_df[currMutSample,])
    
    #filter on current mutational change in current sample ID
    dfMutChange <- mut_change_df[currMutSample,]
    
    # if lcn is a real value
    if(!is.na(dfMutChange$facets.lcn.em.min) & is.numeric(dfMutChange$facets.lcn.em.min)){
      Mi <- dfMutChange$facets.lcn.em.min
    }else if (is.na(dfMutChange$facets.lcn.em.min) & dfMutChange$facets.loh.pier  %in% 1){
      Mi <- 0
    }else{
      Mi <- 1
    }
    
    #### Call ####
    dfSGZ <- SGZ(sample.seq.name = dfMutChange$`Sample Sequencing Name`,
                 gene = dfMutChange$Gene,
                 AAChange = dfMutChange$`hGVS Protein`,
                 Ci = dfMutChange$facets.tcn.em.min,
                 Mi = Mi,
                 p = dfMutChange$facets.purity,
                 n = dfMutChange$FDP,
                 f = dfMutChange$VAF
    )
    
    dfSGZList[[currMutSample]] <- dfSGZ
  }
  
  ### compile mutational changes preds into a single df
  resSGZ <- rbindlist(dfSGZList)
  
  ### merge back with integrated df
  dfIntegSNiPDx <- left_join(dfIntegSNiPDx, resSGZ,
                             by = c("Sample Sequencing Name",
                                    "Gene",
                                    "hGVS Protein" = "AAChange"))
  
  ## define Copy Number
  dfIntegSNiPDx$TCN <- ceiling(dfIntegSNiPDx$facets.tcn.em.median)
  
  #Ploidy Adjusted Copy Number
  dfIntegSNiPDx$PACN <- (dfIntegSNiPDx$TCN / dfIntegSNiPDx$facets.ploidy) * 2
  dfIntegSNiPDx$CRP <- (dfIntegSNiPDx$TCN - dfIntegSNiPDx$facets.ploidy)
  
  
  ### Define CNV - HomDel, Het Loss, Neutral, Gain, Amp
  dfIntegSNiPDx <- dfIntegSNiPDx %>%
    mutate(CNV = case_when(is.na(facets.purity)~"no call",
                           CRP >= 4  ~ "amp",
                           CRP >= 2~ "gain",
                           facets.tcn.em.min == 1 ~"hetloss",
                           facets.tcn.em.min == 0 ~ "homdel",
                           TRUE ~ "neutral"
    )) %>%
    mutate(Zygosity = case_when(is.na(facets.purity)~ "no call",
                                # if facets.lcn.em.min == 0 & facets.tcn.em.min == 0 ~ "homdel" then LOH takes precedence
                                facets.tcn.em.min == 0 ~ "homdel",
                                facets.lcn.em.min == 0 ~ "LOH",
                                is.na(facets.lcn.em.min) & facets.loh.pier == 1 ~ "LOH",
                                is.na(facets.lcn.em.min) & facets.loh.pier == 0 ~ "no call",
                                TRUE ~ "No-LOH")) %>%
    mutate(`Allelic Status` = case_when(CNV == "homdel"~"biallelic",
                                        (SGZ_pred %in% "germline" | SGZ_pred %in% "somatic") & !is.na(VAF) & Zygosity == "LOH"~"biallelic",
                                        (SGZ_pred %in% "germline" | SGZ_pred %in% "somatic") & !is.na(VAF) & Zygosity == "No-LOH"~"monoallelic",
                                        SGZ_pred %in% "somatic subclonal" ~ "subclonal",
                                        SGZ_pred %in% "ambiguous -- poor fit" | SGZ_pred %in% "ambiguous" ~ "no call",
                                        !is.na(VAF) & Zygosity == "no call"~"indeterminate",
                                        TRUE ~ as.character(NA)))
  
  Calls <- dfIntegSNiPDx %>%
    select(any_of(c("Sample Sequencing Name",
                    "Gene",
                    "facets.tcn.em.min",
                    "facets.lcn.em.min",
                    "facets.purity",
                    "facets.ploidy" ,
                    "facets.loh.pier",
                    "CNV",
                    "TCN",
                    "PACN",
                    #copy number relative to ploidy
                    "CRP",
                    "Zygosity",
                    "Allelic Status",
                    "hGVS Protein",
                    "hGVS cDNA",
                    "Rep. Class.",
                    "VAF",
                    colnames(resSGZ)))) %>%
    #remove distinct, only display if protein change is not NA or if CNV is amp, gain, homdel
    filter(!is.na(`hGVS cDNA`) | CNV %in% c("amp", "gain", "homdel"))
  
  # not all sampleIDs are present in Inventory
  dfIntegInventory <- inner_join(dfIntegSNiPDx, Inventory)
  
  if(trial != "Non-Clinical"){
    # top level categorization of sample types
    # dfIntegInventory <- mutate(dfIntegInventory, SAMPTYPE2 =
    #                              case_when(SAMPTYPE == "Normalized ctDNA DNA Extraction" ~ "PBMC",
    #                                        SAMPTYPE == "PBMC" ~ "PBMC",
    #                                        is.na(SAMPTYPE) ~ "NA",
    #                                        TRUE ~ "Tumor"))
    
    ### categorize sample types using test codes
    dfIntegInventory <- mutate(dfIntegInventory, SAMPTYPE2 =
                                 case_when(`Test Code` == "AA099a" ~ "PBMC",
                                           `Test Code` == "AA091a" ~ "Tumor",
                                           is.na(`Test Code`) ~ "NA",
                                           TRUE ~ "Tumor"))
    
    dfIntegInventory <- dfIntegInventory %>%
      mutate(`Test Type` = case_when(SAMPTYPE2 %in% "PBMC" ~ "germline",
                                     SAMPTYPE2 %in% "Tumor" ~ "tissue NGS"))
  }
  
  #### revalue splice variants
  # revalue variants to "splice" if
  # Variant Effect is "." AND p.? or
  # if the c. annotation starts with a number, contains "-", or "+" and ends with number
  dfIntegInventory <- dfIntegInventory %>%
    mutate(`Variant Effect` = case_when(`Variant Effect` %in% "." &
                                          `hGVS Protein` %in% "p.?" |
                                          grepl("c\\.[0-9]+[-|\\+][0-9]",`hGVS cDNA`) ~
                                          "spice-site",
                                        TRUE ~ `Variant Effect`)) %>%
    mutate(`Variant Effect` = case_when(`Variant Effect` %in% "." &
                                          grepl("fs", `hGVS Protein`) ~"frameshift",
                                        TRUE ~ `Variant Effect`))
  
  ## convert to shorten mutation annotations:
  dfIntegInventory <- dfIntegInventory %>%
    mutate(AAChange = case_when(`hGVS Protein` == "p.?"~`hGVS cDNA`,
                                `hGVS Protein` == "p.Ser2018*"~"p.S2018fs",
                                `hGVS Protein` == "p.Gln754Cysfs*9"~"p.L752*",
                                `hGVS Protein` == "p.Lys750="~"p.K750K",
                                TRUE~ShortenMuts(`hGVS Protein`)))
  
  ### BAM LOCATION ON GCP
  trial_to_gcp <- list("RP3500-01" = "350001",
                       "RP3500-03" = "350003",
                       "RP6306-01" = "rp630601",
                       "RP6306-02" = "rp630602",
                       "Non-Clinical" = "nonclin")
  
  sampID_to_gcp <- list("RP3500-01" = "SAMPID2",
                        "RP3500-03" = "SAMPID",
                        "RP6306-01" = "SAMPID",
                        "RP6306-02" = "SAMPID",
                        "Non-Clinical" = "SAMPID2")
  
  trial_to_gcp <- trial_to_gcp[[trial]]
  
  sampID_to_gcp <- sampID_to_gcp[[trial]]
  
  dfIntegInventory$SEQ <- paste0("SEQ_",
                                 gsub(".*SEQ_","", dfIntegInventory$`Sample Sequencing Name`))
  
  dfIntegInventory$GCP_path_to_bam <-
    paste0("gs://quartz-bio/prd/repare_",
           trial_to_gcp,
           "_bmdb/data/source/GENOSITY/ARCHER/",
           dfIntegInventory$SEQ, "/",
           dfIntegInventory[,sampID_to_gcp], "/",
           "BAM/",
           dfIntegInventory$`Sample Sequencing Name`,
           ".final.bam"
    )
  
  
  #### ==== SEGMENT LVL ==== ####
  if(nrow(OrigSegments) >0 ){
    ### separate snipdx and non snipdx genes, expand genes into indiv. rows
    Segments_SNIPDX <- OrigSegments %>%
      mutate(Gene = strsplit(as.character(snipdx_genes_on_segment), split = ","),
             SNIPDX_Gene = "snipdx gene",
             chrom = as.numeric(gsub("chr", "", chrom))) %>%
      tidyr::unnest(Gene) %>%
      select(-c("snipdx_genes_on_segment",
                "other_genes_on_segment")) %>%
      filter(!is.na(Gene))
    
    Segments_NonSNIPDX <- OrigSegments %>%
      mutate(Gene = strsplit(as.character(other_genes_on_segment), split = ","),
             SNIPDX_Gene = "nonsnipdx gene",
             chrom = as.numeric(gsub("chr", "", chrom))) %>%
      tidyr::unnest(Gene) %>%
      select(-c("snipdx_genes_on_segment",
                "other_genes_on_segment")) #%>%
    #filter(!is.na(Gene))
    
    Segments <- rbind(Segments_SNIPDX, Segments_NonSNIPDX)
    
    #### CNV Calls ####
    Segments$Seg_TCN <- ceiling(Segments$total_cn)
    
    #Ploidy Adjusted Copy Number
    Segments$Seg_PACN <- (Segments$Seg_TCN / Segments$ploidy) * 2
    Segments$Seg_CRP <- (Segments$Seg_TCN - Segments$ploidy)
    
    Segments <- Segments %>%
      mutate(Seg_CNV = case_when(is.na(purity)~"no call",
                                 Seg_CRP >= 4  ~ "amp",
                                 Seg_CRP >= 2~ "gain",
                                 total_cn == 1 ~"hetloss",
                                 total_cn == 0 ~ "homdel",
                                 TRUE ~ "neutral")) %>%
      mutate(Zygosity = case_when(is.na(purity)~ "no call",
                                  # if minor_cn == 0 & total_cn == 0 ~ "homdel" then LOH takes precedence
                                  total_cn == 0 ~ "homdel",
                                  minor_cn == 0 ~ "LOH",
                                  is.na(minor_cn) & segments.loh.pier == 1 ~ "LOH",
                                  is.na(minor_cn) & segments.loh.pier == 0 ~ "no call",
                                  TRUE ~ "No-LOH"))
    
    
    #merge segments df with inventory
    Seg_dfInventory <- right_join(Inventory,
                                  Segments)
    
    #merge variants df with inventory
    Seg_dfInventory <- left_join(Seg_dfInventory,
                                 VariantSelected,
                                 by = c("chrom" = "Chr.",
                                        "Gene" = "Gene",
                                        "Sample Sequencing Name" = "Sample Sequencing Name"))
    
    # filter for variants only in Variant df
    Seg_VariantSelected <- Seg_dfInventory %>%
      filter(!is.na(`hGVS cDNA`) ) %>%
      relocate(Gene, chrom,  45:length(colnames(Seg_dfInventory)), 28:44)
    
    # if there are identiical genes on two separate segments, select variant on the correct genomic position
    # iterate through multiple mutational changes within the same sample
    seg_mut_change_df <- Seg_VariantSelected %>%
      group_by(`Sample Sequencing Name`, `hGVS Protein`) %>%
      filter((Pos. >= start & Pos. <= end) & !is.na(`hGVS Protein`)) %>%
      ungroup()
    
    #### ==== variant/segment lvl SGZ ==== ####
    seg_dfSGZList <- list()
    
    if(nrow(seg_mut_change_df) != 0){
      #iterate thru each sample and alteration level, subset, and perform SGZ
      for(currMutSample in seq(nrow(seg_mut_change_df))){
        print(seg_mut_change_df[currMutSample,])
        
        #filter on current mutational change in current sample ID
        dfMutChange <- seg_mut_change_df[currMutSample,]
        
        # if lcn is a real value
        if(!is.na(dfMutChange$minor_cn) & is.numeric(dfMutChange$minor_cn)){
          Mi <- dfMutChange$minor_cn
          ### TODO: segment level pier score
        }else if (is.na(dfMutChange$minor_cn)){
          Mi <- 0
        }else{
          Mi <- 1
        }
        
        #### Call ####
        dfSGZ <- SGZ(sample.seq.name = dfMutChange$`Sample Sequencing Name`,
                     gene = dfMutChange$Gene,
                     AAChange = dfMutChange$`hGVS Protein`,
                     Ci = dfMutChange$total_cn,
                     Mi = Mi,
                     p = dfMutChange$purity,
                     n = dfMutChange$FDP,
                     f = dfMutChange$VAF
        )
        
        seg_dfSGZList[[currMutSample]] <- dfSGZ
      }
      
      ### compile mutational changes preds into a single df
      resSGZ <- rbindlist(seg_dfSGZList)
      
      ### merge back with integrated df
      Seg_dfInventory <- right_join(resSGZ,Seg_dfInventory,
                                    by = c("Sample Sequencing Name",
                                           "Gene",
                                           "AAChange" = "hGVS Protein"))
      
      Seg_dfInventory <- Seg_dfInventory %>%
        relocate(36:length(colnames(Seg_dfInventory)),
                 .after = AAChange)
      
      
    }else{
      #build empty SGZ df
      resSGZ <- data.frame(`Sample Sequencing Name` = Seg_dfInventory$`Sample Sequencing Name`,
                           Gene = NA,
                           AAChange = NA,
                           pred_AF_somatic = NA,
                           pred_AF_germline = NA,
                           pvalue_pred_somatic = NA,
                           pvalue_pred_germline = NA,
                           SGZ_pred = NA,
                           check.names = F)
      
      ### merge back with integrated df
      Seg_dfInventory <- right_join(resSGZ,Seg_dfInventory,
                                    by = c("Sample Sequencing Name",
                                           "Gene",
                                           "AAChange" = "hGVS Protein"))
      
      Seg_dfInventory <- Seg_dfInventory%>%
        relocate(`hGVS cDNA`,
                 34:length(colnames(Seg_dfInventory)),
                 .after = AAChange)
      
    }
    
    #### ==== variant/segment lvl Calls ==== ####
    
    Seg_dfInventory <- Seg_dfInventory%>%
      ### filter only for segments that have variants, amps, and homdel (variants or non snipdx gene)
      filter(!is.na(`hGVS cDNA`) | Seg_CNV %in% c("amp", "homdel"))
    
    ### Define CNV - HomDel, Het Loss, Neutral, Gain, Amp
    Seg_dfInventory <- Seg_dfInventory %>%
      mutate(`Allelic Status` = case_when(Seg_CNV == "homdel"~"biallelic",
                                          (SGZ_pred %in% "germline" | SGZ_pred %in% "somatic") & !is.na(VAF) & Zygosity == "LOH"~"biallelic",
                                          (SGZ_pred %in% "germline" | SGZ_pred %in% "somatic") & !is.na(VAF) & Zygosity == "No-LOH"~"monoallelic",
                                          # if SGZ_pred == subclonal ~ "subclonal"
                                          # if SGZ_pred == "ambig" ~ "no call
                                          # if SGZ_pred =="ambig --poor fit" ~ "no call",
                                          SGZ_pred %in% "somatic subclonal" ~ "subclonal",
                                          SGZ_pred %in% "ambiguous -- poor fit" | SGZ_pred %in% "ambiguous" ~ "no call",
                                          !is.na(VAF) & Zygosity == "no call"~"indeterminate",
                                          TRUE ~ as.character(NA)))
    
    Seg_Calls <- Seg_dfInventory %>%
      select(any_of(c("Sample Sequencing Name",
                      "Gene",
                      "chrom",
                      "hGVS Protein",
                      "hGVS cDNA",
                      "total_cn",
                      "minor_cn",
                      "purity",
                      "ploidy" ,
                      "segments.loh.pier",
                      "Seg_CNV",
                      "Seg_TCN",
                      "Seg_PACN",
                      #copy number relative to ploidy
                      "Seg_CRP",
                      "Zygosity",
                      "Allelic Status",
                      "Rep. Class.",
                      "VAF",
                      colnames(resSGZ)))) #%>%
    #remove distinct, only display if protein change is not NA or if CNV is amp, gain, homdel
    #filter(!is.na(`hGVS cDNA`) | Seg_CNV %in% c("amp", "gain", "homdel"))
    
  }else{
    Segments <- NULL
    Seg_VariantSelected <- NULL
    Seg_Calls <- NULL
    Seg_dfInventory <- NULL
  }
  # #### ==== CONCORDANCE PREPROCESSING ==== ####
  #### ==== NGS ==== ####
  if(trial != "Non-Clinical"){
    # remove na in SAMPTYPE
    ngs <- dfIntegInventory %>%
      tidyr::drop_na(SAMPTYPE)
    
    #### ==== SNV.INDEL ==== ####
    # # snv: synonymous, missense, stop gain/nonsense
    # drop entries that have na in Variant effect
    ngs.snv.indel <- ngs %>%
      tidyr::drop_na(`Variant Effect`,
                     `hGVS Protein`,
                     `hGVS cDNA`)
    
    
    #### ==== SNV.INDEL.Test type ==== ####
    ngs.snv.indel.test.type <- split(ngs.snv.indel, ngs.snv.indel$`Test Type`)
    
    #### ==== CNV ==== ####
    ngs.cnv <- ngs %>%
      filter(CNV %notin% "no call")
    
    #### ==== CNV.Test type ==== ####
    ngs.cnv.test.type <- split(ngs.cnv, ngs.cnv$`Test Type`)
    
    #### ==== STORE ==== ####
    ngs.list <- list(ngs = ngs,
                     snv.indel = c(list(snv.indel = ngs.snv.indel),
                                   ngs.snv.indel.test.type),
                     cnv = c(list(cnv = ngs.cnv),
                             ngs.cnv.test.type)
    )
  }else{
    ngs.list <- list()
  }
  
  #### ==== set s4 ==== ####
  snipdx <- new("snipdx",
                dfIntegInventory = dfIntegInventory,
                Stats = StatsSelected,
                Variant = VariantSelected,
                Facets = FacetsSelected,
                #CNVLong = CNVLong,
                Calls = Calls,
                Segments = Segments,
                Seg_VariantSelected = Seg_VariantSelected,
                Seg_Calls = Seg_Calls,
                Seg_dfInventory = Seg_dfInventory,
                ngs = ngs.list,
                trial = trial,
                date = Sys.Date())
  
  return(snipdx)
  
}

#### SGZ ####

#### AF PRED FUNCTION ####
#### SOMATIC ####
predicted_AF_somatic <- function(Ci,Vi, p){
  #### formula: (pVi)/(pCi + 2(1-p)), where Vi == Ci
  
  #Vi = Mi/Ci or Ci-Mi, variant allele count
  #Vi <- Ci
  
  pred_AF_somatic <- as.numeric(p*Vi)/(p*Ci + 2*(1-p))
  
  return(pred_AF_somatic)
}


#### GERMLINE ####
predicted_AF_germline <- function(Ci,Vi, p){
  #### formula: (pVi + 1 - p)/(pCi + 2(1-p)), where Vi == Ci
  
  #Vi = Mi/Ci or Ci-Mi, variant allele count
  #Vi <- Ci
  
  pred_AF_germline <- (p*Vi + 1 - p)/(p*Ci + 2*(1-p))
  
  return(pred_AF_germline)
}


####  STAT SIGN FUNCTIONS  ####
#### SOMATIC ####
p_somatic <- function(Ci, Vi, p, n, f){
  
  # calculate predicted somatic AF
  pred_AF_somatic <- predicted_AF_somatic(Ci = Ci,Vi = Vi, p = p)
  
  # calculate stat sig of prediction
  ### formula: y ~ binomial(x,n,p) where n = coverage, p = predicated AF,
  # and x = (n*f), where f = VAF
  ### "y is obtained using the 2-tailed binomial test P(y|S; AFsomatic) = Bin (nf, n, AFsomatic)"
  y <- binom.test(x = round(n*f),
                  n = n,
                  p = pred_AF_somatic,
                  alternative = "two.sided")
  
  y <- y$p.value
  
  return(y)
}

#### GERMLINE ####
p_germline <- function(Ci,Vi, p, n,f){
  
  # calculate predicted germline AF
  pred_AF_germline <- predicted_AF_germline(Ci = Ci,Vi = Vi, p = p)
  
  # calculate stat sig of prediction
  ### formula: y ~ binomial(x,n,p) where n = coverage, p = predicated AF, and x = (n*f)
  ### "y is obtained using the 2-tailed binomial test P(y|G; AFgermline)"
  y <- binom.test(x = round(n*f),
                  n = n,
                  p = pred_AF_germline,
                  alternative = "two.sided")
  
  y <- y$p.value
  
  return(y)
}

##### SGZ function #####
### wrap all SGZ calculations into a single function call,
#return estimated somatic and germline VAFs, pvalue for somatic and germline
#p = facets.purity
#f = VAF, mutational allele freq
#Mi = facets.lcn.em.min, minor allele count, this will come up when labelling hetero/homo
#Ci = facets.tcn.em.min, total copy number
#n = variant level coverage

SGZ <- function(sample.seq.name, gene, AAChange,
                Ci, Mi, p, n, f ) {
  
  # a = alpha, constant in SGZ
  a <- 0.01
  
  # CN relative to variant segment, Ci-Mi
  Vi <- Ci-Mi
  
  # set conditions if purity is NA or 0, set to .01 (pseudovalue)
  if(is.na(p)){
    p <- 0.01
  }
  
  pred_AF_somatic <- predicted_AF_somatic(Ci = Ci,Vi = Vi, p = p)
  
  pred_AF_germline <- predicted_AF_germline(Ci = Ci,Vi = Vi, p = p)
  
  y_somatic <- p_somatic(Ci = Ci,Vi = Vi, p = p, n = n, f = f )
  
  y_germline <- p_germline(Ci = Ci, Vi = Vi, p = p, n = n, f = f )
  
  ### set SGZ calls based on conditions -- obsolete as of 2022-08-05
  # if(y_somatic > a & y_germline < a){
  #   outcome <- "somatic"
  # } else if (y_somatic <= a & y_germline > a){
  #   outcome <- "germline"
  # } else if(y_somatic <= a & y_germline <= a & f < (pred_AF_somatic/1.5) & p > .2){
  #   outcome <- "subclonal somatic"
  # }else if(y_somatic <= a & y_germline <= a & f >= (pred_AF_somatic/1.5) | p <= .2){
  #   outcome <- "ambiguous -- poor fit"
  # }else if(y_somatic > a & y_germline > a){
  #   outcome <- "ambiguous"
  # }
  
  ### new implementation, recovery step for germline/somatic
  pvaldifflog10  <- abs(-log(y_somatic,10) - log(y_germline,10))
  
  delta_somatic_vaf  <- ((abs(f - pred_AF_somatic)) / abs(pred_AF_somatic - pred_AF_germline)) -
    ((abs(f - pred_AF_germline) / abs(pred_AF_somatic - pred_AF_germline)))
  
  # update purity from .2 to .3,
  # subclonal somatic --> somatic subclonal
  if(delta_somatic_vaf < -.5 & pvaldifflog10 > 10 & f < (pred_AF_somatic/1.5) & p > .3){
    outcome <- "somatic subclonal"
    # original -- else if(delta_somatic_vaf < -.6 & pvaldifflog10 > 10 )
  }else if((delta_somatic_vaf < -.5 & pvaldifflog10 > 10) | (y_somatic > a & y_germline < 0.001) ){
    outcome <- "somatic"
    #delta_somatic_vaf > .6 & pvaldifflog10 > 10
  }else if((delta_somatic_vaf > .5 & pvaldifflog10 > 10) | (y_somatic <= 0.001 & y_germline > a)){
    outcome <- "germline"
  }else{
    outcome <- "ambiguous -- poor fit"
  }
  
  
  dfSGZ <- data.frame(`Sample Sequencing Name` = sample.seq.name,
                      Gene = gene,
                      AAChange = AAChange,
                      pred_AF_somatic = pred_AF_somatic,
                      pred_AF_germline = pred_AF_germline,
                      pvalue_pred_somatic = y_somatic,
                      pvalue_pred_germline = y_germline,
                      pvaldifflog10 = pvaldifflog10,
                      delta_somatic_vaf = delta_somatic_vaf,
                      SGZ_pred = outcome,
                      check.names = F)
  
  return(dfSGZ)
}



A2AAA <- function(x){
  require(stringr)
  str_replace_all(x,
                  c(
                    "A"="Ala", "R"="Arg", "N"="Asn", "D"="Asp",
                    "C"="Cys", "E"="Glu", "Q"="Gln", "G"="Gly",
                    "H"="His", "I"="Ile", "L"="Leu", "K"="Lys",
                    "M"="Met", "F"="Phe", "P"="Pro", "S"="Ser",
                    "T"="Thr", "W"="Trp", "Y"="Tyr", "V"="Val"))
}

AAA2A <- function(x){
  require(stringr)
  str_replace_all(x,
                  c(
                    "Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D",
                    "Cys"="C", "Glu"="E", "Gln"="Q", "Gly"="G",
                    "His"="H", "Ile"="I", "Leu"="L", "Lys"="K",
                    "Met"="M", "Phe"="F", "Pro"="P", "Ser"="S",
                    "Thr"="T", "Trp"="W", "Tyr"="Y", "Val"="V","Ter"="*"))
}


ShortenMuts <- function(x){
  # Shorten from 3 letter to 1 letter
  x <- AAA2A(x)
  #remove p.
  x <- gsub("^p.","",x)
  
  ### remove AA before fs
  x <- gsub("[A-Z]fs","fs",x)
  
  x <- gsub("[*][0-9]+","",x)
  # x <- gsub("fs[0-9]+","fs",x)
  
  # preserve "p." prefix
  # if alterations do not have a "p." and is not NA --> append prefix, else leave as is
  x <- ifelse( !grepl("^p.",x) & !is.na(x) &  x != "",
               paste0("p.", x),
               x)
  return(x)
}

