library(dplyr)
library(tidyr)
library(snpsettest)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list = ls())
# merge 3 wt samples(shGFP)
wt_edit<-list()
for (i in 1:3){
  wt_edit[[i]]<-read.table(paste("./shGFP/shGFP",i,"CT_filter.edit",sep = "_"),
                           col.names = c("chr","position","ref","alt","ref_num","alt_num"))%>% 
    as_tibble() %>% unite(id,chr,position,sep = ":") %>% dplyr::select(id,ref_num,alt_num)
}
wt_edit_merge<-merge(merge(wt_edit[[1]],wt_edit[[2]],by="id",all = T),wt_edit[[3]],all = T)

wt_edit_merge[is.na(wt_edit_merge)]<-0
wt_edit_merge_2<-as_tibble(wt_edit_merge) %>% mutate(ref_num_all=ref_num.x+ref_num.y+ref_num,alt_num_all=alt_num.x+alt_num.y+alt_num) %>%
  dplyr::select(id,ref_num_all,alt_num_all)
# calculate specific site(subtract wt site)
apo_specific<-list()
for (i in c("A","B","C","D","E")){
  tmp_apo<-read.table(paste("./",i,"/APO_",i,"_CT_filter.edit",sep = ""),
                      col.names = c("chr","position","ref","alt","ref_num","alt_num"))%>% 
    as_tibble() %>% unite(id,chr,position,sep = ":") %>% dplyr::select(id,ref_num,alt_num)
  apo_specific[[i]]<-dplyr::setdiff(tmp_apo[,1],wt_edit_merge_2[,1])
}
# merge different linker length FAM46C and GFP sites
apo_total<-rbind(apo_specific[["A"]],apo_specific[["B"]],apo_specific[["C"]],apo_specific[["D"]],
                 apo_specific[["E"]]) %>% arrange(id) %>% unique()
write.table(file = "apo_total_2.txt",x=apo_total,col.names = FALSE, row.names = FALSE,quote = FALSE) 



apo_total<-list()
apo_merge<-list()
apo_merge_fisher<-list()
apo_merge_map<-list()
ego_apo<-list()
apo_merge_map_genes<-list()
gene_description<-readRDS("D:/NGS/gene_description_grch38.rds")
apo_total[["A"]]<-read.table("./APO_A_total.site",col.names = c("chr","pos","ref","A_ref_num","A_alt_num","A_dup"))
apo_total[["B"]]<-read.table("./APO_B_total.site",col.names = c("chr","pos","ref","B_ref_num","B_alt_num","B_dup"))
apo_total[["C"]]<-read.table("./APO_C_total.site",col.names = c("chr","pos","ref","C_ref_num","C_alt_num","C_dup"))
apo_total[["D"]]<-read.table("./APO_D_total.site",col.names = c("chr","pos","ref","D_ref_num","D_alt_num","D_dup"))
apo_total[["E"]]<-read.table("./APO_E_total.site",col.names = c("chr","pos","ref","E_ref_num","E_alt_num","E_dup"))




dev.off()
pdf("ego_apobec.pdf")
for(i in (c("B","C","D","E"))){
 
  apo_merge[[i]]<-merge(apo_total[["A"]][,c(1,2,4,5,6)],apo_total[[i]][,c(1,2,4,5,6)],by=c("chr","pos")) %>% 
    unite(id,chr,pos,sep = ":")
  
  
  # 创建一个空的数据框，用于存储每个位点的Fisher's Exact Test结果
  fisher_results <- data.frame(site = character(), p_value = numeric())
  # 
  # # 针对每个位点进行Fisher's Exact Test
  for (j in 1:nrow(apo_merge[[i]])) {
    # 提取两个样本在当前位点上的变异和未变异计数
    ref_counts <- c(apo_merge[[i]][j, "A_ref_num"], apo_merge[[i]][j, paste(i,"_ref_num",sep = "")])
    alt_counts <- c(apo_merge[[i]][j, "A_alt_num"], apo_merge[[i]][j,  paste(i,"_alt_num",sep = "")])
    
    # 进行Fisher's Exact Test
    fisher_result <- fisher.test(matrix(c(ref_counts, alt_counts), nrow = 2), alternative = "greater")
    
    # 将结果存储到数据框中
    fisher_results <- rbind(fisher_results, data.frame(apo_merge[[i]][j, 1:7],p_value = fisher_result$p.value))
  }
  
  apo_merge_fisher[[i]]<-fisher_results %>% filter(p_value<0.05&id!="1:117623069") %>% mutate(ID=id)%>%
    separate(id,into = c("chr","pos"),sep = ":")%>% mutate(chr2=paste("chr",chr,sep = ""))
 #  
  write.table(file = paste("apo_merge_AB_",i,".bed",sep=""),x =apo_merge_fisher[[i]][,c("chr2","pos","pos","ID")] ,
              col.names = FALSE, row.names = FALSE,quote = FALSE,sep = "\t")

  
  apo_merge_map[[i]]<-apo_merge_fisher[[i]][,c("ID","chr","pos")]
  colnames( apo_merge_map[[i]])<-c("id","chr","pos")
  tmp_apo_merge_map_genes<-map_snp_to_gene( apo_merge_map[[i]],gene.curated.GRCh38,extend_start = 0,extend_end = 0)[[2]]%>%
    separate(gene.id,"ensembl_gene_id",sep = "\\.")
  apo_merge_map_genes[[i]]<-merge(tmp_apo_merge_map_genes,gene_description,by="ensembl_gene_id",all.x=T)
 # go enrichment 
  for (k in c("CC","BP","MF")){
    tmp_ego<-enrichGO(gene =  apo_merge_map_genes[[i]]$ensembl_gene_id,
                          OrgDb = org.Hs.eg.db,
                          keyType='ENSEMBL',
                          ont = k,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
    try(print(dotplot(tmp_ego,title=paste("ego",i,k,sep = "_"),showCategory=10,font=8)),silent=TRUE)
    ego_apo[[i]][[k]]<-tmp_ego %>% data.frame()
    #try(print(dotplot(ego_apo[[i]][[k]],title=paste("ego",i,k,sep = "_"),showCategory=10,font=8)),silent=TRUE)
  }
}
dev.off()

