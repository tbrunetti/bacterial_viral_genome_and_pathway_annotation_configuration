library(clusterProfiler)
library(readxl)

mappings = read.csv("staph_aureus_kegg_mappings.tsv", sep = "\t")

excel_sheets <- readxl::excel_sheets("202303_TnSeq_Significant_Hits_for_KEGG.xlsx")

valid_excel_sheets <- excel_sheets[1:length(excel_sheets)]
for (i in valid_excel_sheets){
  de_results <- read_excel("202303_TnSeq_Significant_Hits_for_KEGG.xlsx", col_names =FALSE, sheet = i)
  names(de_results)[1]  <- "locus_tag"
  merged_kegg <- merge(de_results, mappings, by = "locus_tag", all.x = T)
  kk <- enrichKEGG(gene = merged_kegg$ko_ids, organism = 'ko', pvalueCutoff = 0.05, keyType = "kegg")
  plot <- dotplot(kk, title = i)
  ggsave(filename = paste(i, "_kegg_pathway_ora.pdf", sep = ""), plot = plot, device = "pdf")
  hits_df <- kk@result
  hits_df$kegg_id <- hits_df$geneID
  hits_df$locus_tag <- hits_df$geneID
  for (h in seq(1:nrow(hits_df))){
    geneids <- c()
    locusids <- c()
    orig <- unlist(str_split(hits_df[h, 'geneID'], '/'))
    for (j in orig){
      if (j == ""){
        next
      }else{
        print(j)
        map_index <- which(mappings$ko_ids %in% j)
        translations <- mappings[map_index, c("gene", "locus_tag", "ko_ids")]
        geneids <- c(geneids, translations$gene)
        locusids <- c(locusids, translations$locus_tag)      
      }
    }
    hits_df[h, "geneID"] <- paste(geneids, collapse = "/")
    hits_df[h, "locus_tag"] <- paste(locusids, collapse = "/")
  }
  write.csv(hits_df, file = paste(i, "_kegg_pathway_ora_pathway_analysis_results.csv", sep = ""))
}

# csv <- "~/7v4L_UniqueHits2.csv"
# csv_data <- read.csv(csv, skip = 1)
# de_results <- read.csv(csv, skip = 1)
# locustagloc <- which(colnames(de_results) == "old.locus.tag")
# names(de_results)[locustagloc]  <- "locus_tag"
# merged_kegg <- merge(de_results, mappings, by = "locus_tag", all.x = T)
# kk <- enrichKEGG(gene = merged_kegg$ko_ids, organism = 'ko', pvalueCutoff = 0.05, keyType = "kegg")
# plot <- dotplot(kk, title = "7v4L Unique Hits Underrepresented")
# ggsave(filename = paste(i, "7v4L_Unique_Hits_Underrepresented_kegg_pathway_ora.pdf", sep = ""), plot = plot, device = "pdf")

# hits_df <- kk@result
# hits_df$kegg_id <- hits_df$geneID
# hits_df$locus_tag <- hits_df$geneID
# 
# for (i in seq(1:nrow(hits_df))){
#   geneids <- c()
#   locusids <- c()
#   orig <- unlist(str_split(hits_df[i, 'geneID'], '/'))
#   for (j in orig){
#     if (j == ""){
#       continue
#     }else{
#       map_index <- which(mappings$ko_ids %in% j)
#       translations <- mappings[map_index, c("gene", "locus_tag", "ko_ids")]
#       geneids <- c(geneids, translations$gene)
#       locusids <- c(locusids, translations$locus_tag)      
#     }
#   }
#   print(geneids)
#   print(locusids)
#   hits_df[i, "geneID"] <- paste(geneids, collapse = "/")
#   hits_df[i, "locus_tag"] <- paste(locusids, collapse = "/")
# }
# 
# 
# write.csv(hits_df, "7v4L_UniqueHits2_ora_pathway_analysis_results.csv", sep = ",")
# 
# 
# 

# gsea <- merged_kegg$log2FC
# names(gsea) <- merged_kegg$ko_ids
# gsea[which(names(gsea) != "")] -> drop_empty
# sorted_gsea <- sort(gsea, decreasing = T)
# 
# kk2 <- gseKEGG(geneList     = sorted_gsea,
#                organism     = 'ko',
#                pvalueCutoff = 0.05,
#                verbose      = FALSE)
# head(kk2)
# 
# browseKEGG(kk2, 'map01100')
# 
# library("pathview")
# hsa04110 <- pathview(gene.data  = sorted_gsea,
#                      pathway.id = "map01100",
#                      species    = "ko",
#                      limit      = list(gene=max(abs(sorted_gsea)), cpd=1))
