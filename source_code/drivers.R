###combine driver mutations of cancers
library(dplyr)
files <- list.files("data/Onco_mutations_OncoVar_TCGA/",full.names = T)
dt <- do.call(rbind,lapply(files,data.table::fread,data.table=F))

files_gene <- list.files("data/Onco_genes_OncoVar_TCGA/",full.names = T)
dt_genes <- do.call(rbind,lapply(files_gene,data.table::fread,data.table=F))
###
dt <- dt %>% filter(Gene_symbol %in% dt_genes$Gene_symbol)

dt %>% group_by(Cancer,Gene_symbol) %>% summarise(affect=sum(TCGA_affected)) -> a
saveRDS(dt,file = "data/TCGA_driver_mutations.rds")
##看driver 基因中的 driver mutations 的新抗原
unique(dt$Hg38_Location) %>% length() ##[1] 3068

##HLA common and well-documented alleles in China
hla <- read.csv("~/common_driver/data/common_hla_chinese.csv")
hla_1_common <- hla %>%
  filter(Category=="common") %>%
  filter(!grepl("DRB",Allele)) %>%
  filter(!grepl("DQB",Allele))
hla_1_common <- hla_1_common %>%
  mutate(Allele=paste0("HLA-",Allele)) 
hla_1_common$fre <- hla_1_common$Number/812211
hla_1_common <- hla_1_common %>% select(-Locus)
saveRDS(hla_1_common,file = "data/HLA_I_alleles.rds")

###预测
library(MHCbinding)
library(purrr)
library(dplyr)
library(parallel)
TCGA_driver_mutations <- readRDS("~/MHCbindshiny/data/TCGA_driver_mutations.rds")
HLA_I_alleles <- readRDS("~/MHCbindshiny/data/HLA_I_alleles.rds")
dt <- TCGA_driver_mutations %>% 
  select(Hg38_Location,Ref,Alt) %>% 
  distinct(Hg38_Location,.keep_all = T) %>% 
  tidyr::separate(col = Hg38_Location,into = c("chr","position"),sep = ":") %>% 
  tidyr::separate(col = position,into = c("start","end"),sep = "-")
write.table(dt,file = "data/TCGA_driver_mutations_annovar_input.txt",sep = "\t",col.names = F,row.names = F,quote = F)

available_alleles("MHC-I","netmhcpan_ba") -> a
HLA_I_alleles <- HLA_I_alleles %>% 
  filter(Allele %in% a$alleles)
res <- vector("list",nrow(HLA_I_alleles))
names(res) <- HLA_I_alleles$Allele
# for (i in seq_along(res)){
#   res[[i]] <- txt2binding(annovar_path = "~/software/annovar/",
#                           txt_path = "data/TCGA_driver_mutations_annovar_input.txt",
#                           genome_version = "hg38",mhc_type = "MHC-I",pep_length = c(8,9,10,11,12,13,14),
#                           allele = names(res)[i],pre_method = "ann",get_method = "client",client_path = "~/software/mhc_i/src/")
#   
# }
library(parallel)
cl <- makeCluster(getOption("cl.cores", 35), type = "FORK")
res <- parLapply(cl = cl, names(res), function(x) {
  dt <- txt2binding(annovar_path = "~/software/annovar/",
                    txt_path = "data/TCGA_driver_mutations_annovar_input.txt",
                    genome_version = "hg38",mhc_type = "MHC-I",pep_length = c(8,9,10,11,12,13,14),
                    allele = x,pre_method = "netmhcpan_ba",
                    get_method = "client",client_path = "~/software/mhc_i/src/",
                    tmp_dir = paste0("~/tmp/mhcbinding/",x,"/"),num_thread = 1)
  return(dt)
})
stopCluster(cl)
saveRDS(res,file = "~/tmp/TCGA_driver_res2.rds")

TCGA_driver_res <- readRDS("~/tmp/TCGA_driver_res.rds")
drivers <- bind_rows(TCGA_driver_res)
saveRDS(drivers,file = "../tmp/TCGA_drivers_all.rds")

###按照基因-突变来预测
HLA_I_alleles <- readRDS("~/MHCbindshiny/data/HLA_I_alleles.rds")
available_alleles("MHC-I","netmhcpan_ba") -> a
HLA_I_alleles <- HLA_I_alleles %>% 
  filter(Allele %in% a$alleles)
dt <- TCGA_driver_mutations %>% 
  select(Hg38_Location,Ref,Alt,Gene_symbol) %>% 
  distinct_all(.keep_all = T) %>% 
  tidyr::separate(col = Hg38_Location,into = c("chr","position"),sep = ":") %>% 
  tidyr::separate(col = position,into = c("start","end"),sep = "-")
genes <- unique(dt$Gene_symbol)
for (i in 1:length(genes)){
  dir.create(paste0("data/driver/",genes[i]))
  gene_dt <- dt %>% 
    filter(Gene_symbol == genes[i]) %>% 
    select(-Gene_symbol)
  write.table(gene_dt,file = paste0("data/driver/",genes[i],"/anno_input.txt"),
              sep = "\t",col.names = F,row.names = F,quote = F)
}
alleles <- HLA_I_alleles$Allele
library(parallel)

a <- list.files("data/driver/",recursive = T) %>% `[`(.,grepl("neo.rds",.)) %>% gsub("/neo.rds","",.)
genes <- genes[which(!(genes %in% a))]
cl <- makeCluster(getOption("cl.cores", 18), type = "FORK")
for (i in 1:length(genes)){
  res <- parLapply(cl = cl, alleles, function(x,y) {
    dt <- txt2binding(annovar_path = "~/software/annovar/annovar/",
                      txt_path = paste0("data/driver/",y,"/anno_input.txt"),
                      genome_version = "hg38",mhc_type = "MHC-I",pep_length = c(8,9,10,11,12,13,14),
                      allele = x,pre_method = "netmhcpan_ba",
                      get_method = "client",client_path = "~/software/mhc_i/src/",
                      tmp_dir = paste0("~/tmp/mhcbinding/",y,"/",x,"/"),num_thread = 1)
    return(dt)
  },genes[i])
  res <- bind_rows(res)
  saveRDS(res,file = paste0("data/driver/",genes[i],"/neo.rds"))
}
stopCluster(cl)

dt <- txt2binding(annovar_path = "~/software/annovar/annovar/",
                  txt_path = paste0("~/MHCbindshiny/data/driver/","B2M","/anno_input.txt"),
                  genome_version = "hg38",mhc_type = "MHC-I",pep_length = c(8,9,10,11,12,13,14),
                  allele = "HLA-C*14:02",pre_method = "netmhcpan_ba",
                  get_method = "client",client_path = "~/software/mhc_i/src/",
                  tmp_dir = paste0("~/tmp/mhcbinding/","B2M","/","HLA-C*14:02","/"),num_thread = 1)
which(genes == "CNOT3")


###anchar
genes <- list.files("data/driver",full.names = T)
res <- vector("list",476)
for (i in 1:length(genes)){
  t <- readRDS(paste0(genes[i],"/neo.rds"))
  if(nrow(t)==0){
    next
  }
  t <- t %>% 
    select(pos_alter,transcript) %>% 
    distinct_all(.keep_all = T)
  tmp <- vector("list",length = nrow(t))
  for (j in seq_along(tmp)){
    tmp_t <- data.frame(aa = strsplit(t$pos_alter[j],split = ",")[[1]],
                        trans = strsplit(t$transcript[j],split = ",")[[1]])
    tmp[[j]] <- tmp_t
  }
  tmp <- bind_rows(tmp)
  tmp <- tmp %>% 
    group_by(aa) %>% 
    summarise(transcript=paste(trans,collapse = ",")) %>% 
    ungroup() %>% 
    mutate(gene=gsub("data/driver/","",genes[i]))
  res[[i]] <- tmp
  
}
res <- bind_rows(res)
saveRDS(res,file = "data/all_gene_aa_trans.rds")

library(vioplot)
dt_l <- data.frame(x=rnorm(100),y=c(rep("a",100)))
dt_r <- data.frame(x=rnorm(100,mean = 1),y=c(rep("b",100)))
vioplot(x~y, data=dt_l, col = "palevioletred", plotCentre = "line", side = "right")
vioplot(x~y, data=dt_r, col = "lightblue", plotCentre = "line", side = "left", add = T)
