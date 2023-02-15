library(dplyr)
library(ggplot2)
library(ggprism)
neodb <- readxl::read_xlsx("data/neodb_data.xlsx")
unique(neodb$pmid) %>% length()
neodb$`aa change` <- gsub("p[.]","",neodb$`aa change`)
neodb$`HLA allele` <- gsub("HLA-","",neodb$`HLA allele`)
neodb$`HLA allele` <- gsub("[*]","",neodb$`HLA allele`)
neodb$`HLA allele` <- gsub("[ ]","",neodb$`HLA allele`)
table(neodb$`HLA allele`) %>% as.data.frame() -> hla
neodb$pep_len <- nchar(neodb$`Mut Epitope`)
table(neodb$pep_len)
library(xlsx)
write.xlsx(neodb,"data/clean_neodb.xlsx",showNA = FALSE)
##把 HLA allele 人工整理一下
neodb <- readxl::read_xlsx("~/MHCbindshiny/data/clean_neodb.xlsx")
neodb <- neodb %>% select(-`T cell type`,-c(14:20)) %>% 
  distinct_all(.keep_all = T)

# hla$Var1 <- as.character(hla$Var1)
# strsplit(hla$Var1[1],"")
# strsplit(hla$Var1[2],"")[[1]][7]
# library(rebus)
# literal(strsplit(hla$Var1[2],"")[[1]][7])

neodb$HLA <- gsub("\\Q​\\E","",neodb$HLA)
##pep len HLA
neodb$`Mut Epitope` <- gsub(" ","",neodb$`Mut Epitope`) %>% gsub("\\Q​\\E","",.)
neodb$`WT peptide` <- gsub(" ","",neodb$`WT peptide`) %>% gsub("\\Q​\\E","",.)
pep <- neodb %>% select(`Mut Epitope`,`WT peptide`,HLA) %>% 
  distinct_all(.keep_all = T) %>% 
  mutate(pep_len=nchar(`Mut Epitope`)) %>% 
  arrange(pep_len) %>% 
  mutate(hla_type = case_when(
    grepl("^A",HLA) ~ "A",
    grepl("^B",HLA) ~ "B",
    grepl("^C",HLA) ~ "C",
    grepl("^DP",HLA) ~ "DP",
    grepl("^DQ",HLA) ~ "DQ",
    grepl("^DR",HLA) ~ "DR"
  )) %>% 
  mutate(super_type=ifelse(hla_type %in% c("A","B","C"),"HLA-I","HLA-II"))
pep$pep_len <- factor(pep$pep_len,levels = unique(pep$pep_len))
p2 <- ggplot(data=pep,aes(x=pep_len))+
  geom_bar()+
  theme_prism()+
  scale_y_continuous(expand = c(0,0))+
  labs(x="Peptide length",y="Counts")+
  facet_grid(~super_type,scales = "free",space = "free")

p2
table(pep$HLA) %>% as.data.frame() -> hla

hla <- hla %>% 
  mutate(hla_type = case_when(
    grepl("^A",Var1) ~ "A",
    grepl("^B",Var1) ~ "B",
    grepl("^C",Var1) ~ "C",
    grepl("^DP",Var1) ~ "DP",
    grepl("^DQ",Var1) ~ "DQ",
    grepl("^DR",Var1) ~ "DR"
  )) %>% 
  mutate(super_type=ifelse(hla_type %in% c("A","B","C"),"HLA-I","HLA-II"))
hla %>% 
  group_by(super_type,hla_type) %>% 
  summarise(counts=sum(Freq)) -> summ
library(ggplot2)
library(ggprism)
p1 <- ggplot(data=summ,aes(x=hla_type,y=counts))+
  geom_bar(stat = "identity")+
  facet_grid(~super_type,scales = "free_x")+
  scale_y_continuous(expand = c(0,0))+
  theme_prism()+labs(y="Counts",x="HLA alleles")

library(patchwork)
p1+p2+plot_layout(widths = c(1,2))
###IC50 预测
pep <- neodb %>% 
  select(`Mut Epitope`,`WT peptide`,HLA) %>% 
  filter(!is.na(`WT peptide`)) %>% 
  distinct_all(.keep_all = T)%>% 
  mutate(pep_len=nchar(`Mut Epitope`))
library(MHCbinding)

pep <- pep %>% mutate(hla_type = case_when(
  grepl("^A",HLA) ~ "A",
  grepl("^B",HLA) ~ "B",
  grepl("^C",HLA) ~ "C",
  grepl("^DP",HLA) ~ "DP",
  grepl("^DQ",HLA) ~ "DQ",
  grepl("^DR",HLA) ~ "DR"
))

pep %>% 
  mutate(type = case_when(
    hla_type %in% c("A","B","C") & pep_len <= 15 & pep_len >=8 ~ "yes_i",
    hla_type %in% c("DP","DR","DQ") & pep_len <= 30 & pep_len >= 11 ~ "yes_ii",
    TRUE ~ "no"
  )) -> pep
pep <- pep %>% filter(type != "no")
pep <- pep %>% 
  filter((type == "yes_ii" & nchar(HLA)==9) | (type == "yes_i" & nchar(HLA)==6))
which(grepl("DQA",pep$HLA))
pep$HLA[309] <- paste0(pep$HLA[309],"/",pep$HLA[310])
pep <- pep[-310,]
pep <- pep %>% 
  rowwise() %>% 
  mutate(HLA_1 = ifelse(type=="yes_i",
                      sub("(.{1})(.*)", "\\1*\\2", HLA) %>% paste0("HLA-",.),
                      sub("(.{4})(.*)", "\\1*\\2", HLA)))
pep$HLA_1[309] <- "DQA1*03:01/DQB1*03:02"

pep$len_wt <- nchar(pep$`WT peptide`)
which(pep$len_wt != pep$pep_len)
pep$`WT peptide`[135] <- "GIINFYTALL"
pep$`WT peptide`[164] <- "KMLLTEILLI"
pep$`WT peptide`[237] <- "ARDIYRASYY"
pep$`WT peptide`[400] <- "CISSCNPNP"
pep$`WT peptide`[490] <- "DPPALASKNAEVT"
pep$`WT peptide`[491] <- "TYDTVHRRL"
pep$`WT peptide`[494] <- "GYNSYSVSNSEKDIM"
pep$ic50_wt <- 1
pep$rank_wt <- 1
pep$ic50_mt <- 1
pep$rank_mt <- 1

pep$pep_len <- nchar(pep$`Mut Epitope`)
res <- vector("list",nrow(pep))
for (i in seq_along(res)){
  if (pep$type[i] == "yes_i"){
    t <- try(mhcIbinding_client(peptide = pep$`Mut Epitope`[i],
                                allele = pep$HLA_1[i],length = pep$pep_len[i],
                                pre_method = "netmhcpan_ba",client_path="~/software/mhc_i/src/",
                                tmp_dir=tempdir()),silent = TRUE)
    t2 <- try(mhcIbinding_client(peptide = pep$`WT peptide`[i],
                                 allele = pep$HLA_1[i],length = pep$pep_len[i],
                                 pre_method = "netmhcpan_ba",client_path="~/software/mhc_i/src/",
                                 tmp_dir=tempdir()),silent = TRUE)
    if ('try-error' %in% class(t)){
      next
    }else{
      pep$ic50_mt[i] <- t$ic50
      pep$rank_mt[i] <- t$rank
      pep$ic50_wt[i] <- t2$ic50
      pep$rank_wt[i] <- t2$rank
    }
  }else{
    t <- try(mhcIIbinding_client(peptide = pep$`Mut Epitope`[i],
                                 allele = pep$HLA_1[i],length = pep$pep_len[i],
                                 pre_method = "netmhciipan_ba",client_path="~/software/mhc_ii/",
                                 tmp_dir=tempdir()),silent = TRUE)
    t2 <- try(mhcIIbinding_client(peptide = pep$`WT peptide`[i],
                                  allele = pep$HLA_1[i],length = pep$pep_len[i],
                                  pre_method = "netmhciipan_ba",client_path="~/software/mhc_ii/",
                                  tmp_dir=tempdir()),silent = TRUE)
    if ('try-error' %in% class(t)){
      next
    }else{
      pep$ic50_mt[i] <- t$ic50
      pep$rank_mt[i] <- t$percentile_rank
      pep$ic50_wt[i] <- t2$ic50
      pep$rank_wt[i] <- t2$percentile_rank
    }
  }
  
}
saveRDS(pep,file = "data/neodb_pre.rds")
library(ggpubr)
pep <- readRDS("~/MHCbindshiny/data/neodb_pre.rds")
dt <- pep %>% 
  tidyr::pivot_longer(cols = c("ic50_wt","ic50_mt"),names_to = "ic50_type",values_to = "ic50")
p3 <- ggplot(data=dt,aes(x=ic50_type,y=log(ic50+1)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = c("MT","WT"))+
  labs(y="Log(IC50 + 1)")

dt <- pep %>% 
  tidyr::pivot_longer(cols = c("rank_wt","rank_mt"),names_to = "rank_type",values_to = "rank")
p4 <- ggplot(data=dt,aes(x=rank_type,y=log(rank+1)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = c("MT","WT"))+
  labs(y="Log(Rank + 1)")

p3 + p4 

#看看dbpepNeo2.0里面的数据
files <- list.files("data/other",full.names = T)
res <- lapply(files,read.csv)
res <- bind_rows(res)
##没有的拿过来
res <- res %>% 
  tidyr::separate(col = "PMID",into = c("ss","pmid"),sep = ":")
which(!(res$pmid %in% neodb$pmid)) %>% length()
##保存一下
which(grepl("DQA",neodb$HLA))
neodb$HLA[405] <- paste0(neodb$HLA[405],"/",neodb$HLA[406])
neodb <- neodb[-406,]

neodb$len_wt <- nchar(neodb$`WT peptide`)
neodb$pep_len <- nchar(neodb$`Mut Epitope`)
which(neodb$len_wt != neodb$pep_len)
neodb$`WT peptide`[160] <- "GIINFYTALL"
neodb$`WT peptide`[189] <- "KMLLTEILLI"
neodb$`WT peptide`[285] <- "ARDIYRASYY"
neodb$`WT peptide`[288] <- "HHGGWTTKM"
neodb$`WT peptide`[289] <- "MKQMNDAHH"
neodb$`WT peptide`[290] <- "RDIKNDSNYV"
neodb$`WT peptide`[291] <- "RDIKNDSNYV"
neodb$`WT peptide`[514] <- "CISSCNPNP"
neodb$`WT peptide`[641] <- "DPPALASKNAEVT"
neodb$`WT peptide`[642] <- "TYDTVHRRL"
neodb$`WT peptide`[646] <- "GYNSYSVSNSEKDIM"
res$pmid[which(!(res$pmid %in% neodb$pmid))] %>% unique()

neodb$pmid <- gsub(" ","",neodb$pmid)
neodb$pmid <- gsub("\\Q​\\E","",neodb$pmid)
res$pmid[which(!(res$pmid %in% neodb$pmid))] %>% unique() -> a
dt <- res %>% filter(pmid %in% a)
dt <- dt %>% filter(nchar(WT_peptide)>1)
dt$pmid %>% unique()

neodb %>% 
  select(-`HLA allele`,-`transcript id`) %>% 
  select(pmid,ID,"Mut Epitope","Immunizing peptide","aa change",HLA,everything()) -> neodb
saveRDS(neodb,file = "data/clean_neodb_part1.rds")
xlsx::write.xlsx(neodb,file = "data/clean_neodb_part1.xlsx",showNA = FALSE)
##II类的
res <- vector("list",9)
files <- list.files("data/other/MHCII",full.names = T)
res <- lapply(files,function(x){
  df <- read.csv(x)
  df$Mut_Affinity.nM. <- as.numeric(df$Mut_Affinity.nM.)
  df$Rank. <- as.numeric(df$Rank.)
  return(df)
  }
)
res <- bind_rows(res)
res <- res %>% 
  tidyr::separate(col = "PMID",into = c("ss","pmid"),sep = ":")
res$pmid[which(!(res$pmid %in% neodb$pmid))] %>% unique()

##
neodb_2 <- readxl::read_xlsx("data/neodb_part2.xlsx")
neodb_2$pmid <- gsub(" ","",neodb_2$pmid)
neodb_2$pmid <- gsub("\\Q​\\E","",neodb_2$pmid)
neodb_2$`aa change` <- gsub("p[.]","",neodb_2$`aa change`)
neodb_2$HLA <- gsub("HLA-","",neodb_2$HLA)
neodb_2$HLA <- gsub("[*]","",neodb_2$HLA)
neodb_2$HLA <- gsub("[ ]","",neodb_2$HLA)
neodb_2$HLA <- gsub("\\Q​\\E","",neodb_2$HLA)
neodb_2$`Mut Epitope` <- gsub(" ","",neodb_2$`Mut Epitope`) %>% gsub("\\Q​\\E","",.)
neodb_2$`WT peptide` <- gsub(" ","",neodb_2$`WT peptide`) %>% gsub("\\Q​\\E","",.)

neodb_2$pep_len <- nchar(neodb_2$`Mut Epitope`)
neodb_2$len_wt <- nchar(neodb_2$`WT peptide`)
neodb_2$`WT peptide`[52] <- "NLLGRNSFEN"

##将没有WT的补全
library(ensembldb)
library(EnsDb.Hsapiens.v105)
edb <- EnsDb.Hsapiens.v105

neodb_all <- bind_rows(neodb,neodb_2)
neodb_all$pep_len <- nchar(neodb_all$`Mut Epitope`)
neodb_all$len_wt <- nchar(neodb_all$`WT peptide`)
which(neodb_all$pep_len != neodb_all$len_wt)
neodb_all$gene <- gsub("\\Q​\\E","",neodb_all$gene)
neodb_all$gene <- gsub("[ ]","",neodb_all$gene)
neodb_all$gene[which(neodb_all$gene == "GPR133")] <- "ADGRD1"
neodb_all$gene[which(neodb_all$gene == "WDR63")] <- "DNAI3"
neodb_all$gene[which(neodb_all$gene == "PNMAL1")] <- "PNMA8A"
neodb_all$gene[which(neodb_all$gene == "C1orf170")] <- "PERM1"
neodb_all$gene[which(neodb_all$gene == "KIAA1683")] <- "IQCN"
neodb_all <- neodb_all %>% 
  mutate(index=row_number())
na_dt <- neodb_all %>% 
  dplyr::filter(is.na(`WT peptide`)) %>% 
  dplyr::filter(grepl("[A-Z][0-9]+[A-Z]",`aa change`)) %>% 
  dplyr::filter(`aa change` != "R119Gfs*38​")
for (i in 1:nrow(na_dt)){
  mut_pep <- na_dt$`Mut Epitope`[i]
  aa_change <- na_dt$`aa change`[i] %>% gsub("[0-9]","",.)
  gene <- na_dt$gene[i]
  prts <- proteins(edb, filter = GeneNameFilter(gene),
                   return.type = "data.frame")
  index_mut <- which(strsplit(mut_pep,"")[[1]] == substr(aa_change,2,2))
  for (j in index_mut){
    tmp <- mut_pep
    substr(tmp,j,j) <- substr(aa_change,1,1)
    a <- str_extract(string = prts$protein_sequence,pattern = tmp)
    if (all(is.na(a))){
      next
    }else{
      na_dt$`WT peptide`[i] <- tmp
    }
  }
}
neodb_all <- neodb_all %>% dplyr::filter(!(index %in% na_dt$index))
neodb_all <- bind_rows(neodb_all,na_dt)
which(neodb_all$gene == "CALR")
which(neodb_all$gene == "CALRp7")
neodb_all$`aa change`[198] <- "K385Nfs*47"
neodb_all$`aa change`[199] <- "K385Nfs*47"
which(neodb_all$gene=="ARMT1")
neodb_all$`aa change`[247] <- "P286L"
neodb_all$`WT peptide`[247] <- "FYGKTIPWF"
which(neodb_all$gene=="Caspase-5")
neodb_all$`aa change`[309] <- "FS"
saveRDS(neodb_all,file = "data/neodb_all.rds")
xlsx::write.xlsx(neodb_all,file = "data/neodb_all.xlsx",showNA = FALSE)
###plot
library(dplyr)
library(ggplot2)
library(ggprism)
neodb_all <- readRDS("data/neodb_all.rds")
pep <- neodb_all %>% select(`Mut Epitope`,`WT peptide`,HLA) %>% 
  distinct_all(.keep_all = T) %>% 
  mutate(pep_len=nchar(`Mut Epitope`)) %>% 
  arrange(pep_len) %>% 
  mutate(hla_type = case_when(
    grepl("^A",HLA) ~ "A",
    grepl("^B",HLA) ~ "B",
    grepl("^C",HLA) ~ "C",
    grepl("^DP",HLA) ~ "DP",
    grepl("^DQ",HLA) ~ "DQ",
    grepl("^DR",HLA) ~ "DR"
  )) %>% 
  mutate(super_type=ifelse(hla_type %in% c("A","B","C"),"HLA-I","HLA-II"))
pep$pep_len <- factor(pep$pep_len,levels = unique(pep$pep_len))
p2 <- ggplot(data=pep,aes(x=pep_len))+
  geom_bar()+
  theme_prism()+
  scale_y_continuous(expand = c(0,0))+
  labs(x="Peptide length",y="Counts")+
  facet_grid(~super_type,scales = "free",space = "free")

p2
table(pep$HLA) %>% as.data.frame() -> hla

hla <- hla %>% 
  mutate(hla_type = case_when(
    grepl("^A",Var1) ~ "A",
    grepl("^B",Var1) ~ "B",
    grepl("^C",Var1) ~ "C",
    grepl("^DP",Var1) ~ "DP",
    grepl("^DQ",Var1) ~ "DQ",
    grepl("^DR",Var1) ~ "DR"
  )) %>% 
  mutate(super_type=ifelse(hla_type %in% c("A","B","C"),"HLA-I","HLA-II"))
hla %>% 
  group_by(super_type,hla_type) %>% 
  summarise(counts=sum(Freq)) -> summ
library(ggplot2)
library(ggprism)
p1 <- ggplot(data=summ,aes(x=hla_type,y=counts))+
  geom_bar(stat = "identity")+
  facet_grid(~super_type,scales = "free_x")+
  scale_y_continuous(expand = c(0,0))+
  theme_prism()+labs(y="Counts",x="HLA alleles")
p1

##预测
pep <- neodb_all %>% 
  select(`Mut Epitope`,`WT peptide`,HLA) %>% 
  filter(!is.na(`WT peptide`)) %>% 
  distinct_all(.keep_all = T)%>% 
  mutate(pep_len=nchar(`Mut Epitope`))
library(MHCbinding)

pep <- pep %>% mutate(hla_type = case_when(
  grepl("^A",HLA) ~ "A",
  grepl("^B",HLA) ~ "B",
  grepl("^C",HLA) ~ "C",
  grepl("^DP",HLA) ~ "DP",
  grepl("^DQ",HLA) ~ "DQ",
  grepl("^DR",HLA) ~ "DR"
))

pep %>% 
  mutate(type = case_when(
    hla_type %in% c("A","B","C") & pep_len <= 15 & pep_len >=8 ~ "yes_i",
    hla_type %in% c("DP","DR","DQ") & pep_len <= 30 & pep_len >= 11 ~ "yes_ii",
    TRUE ~ "no"
  )) -> pep
pep <- pep %>% filter(type != "no")
pep <- pep %>% 
  filter((type == "yes_ii" & nchar(HLA)==9) | (type == "yes_i" & nchar(HLA)==6) | (type == "yes_ii" & nchar(HLA)==19))
pep <- pep %>% 
  rowwise() %>% 
  mutate(HLA_1 = ifelse(type=="yes_i",
                        sub("(.{1})(.*)", "\\1*\\2", HLA) %>% paste0("HLA-",.),
                        sub("(.{4})(.*)", "\\1*\\2", HLA)))

which(pep$HLA_1 == "DQA1*:0301/DQB1:0302")
pep$HLA_1[310] <- "DQA1*03:01/DQB1*03:02"
which(pep$HLA_1 == "DPA1*02:02/DPB103:01")
pep$HLA_1[c(566,594)] <- "DPA1*02:02/DPB1*03:01"
which(pep$HLA_1 == "DQA1*05:01/DQB102:01")
pep$HLA_1[565] <- "DQA1*05:01/DQB1*02:01"
which(pep$HLA_1 == "DPA1*03:01/DPB103:01")
pep$HLA_1[595] <- "DPA1*03:01/DPB1*03:01"
which(pep$HLA_1 == "DPA1*02:01/DPB109:01")
pep$HLA_1[590] <- "DPA1*02:01/DPB1*09:01"

pep$len_wt <- nchar(pep$`WT peptide`)
which(pep$len_wt != pep$pep_len)
pep$ic50_wt <- 1
pep$rank_wt <- 1
pep$ic50_mt <- 1
pep$rank_mt <- 1

pep$pep_len <- nchar(pep$`Mut Epitope`)
res <- vector("list",nrow(pep))
for (i in seq_along(res)){
  if (pep$type[i] == "yes_i"){
    t <- try(mhcIbinding_client(peptide = pep$`Mut Epitope`[i],
                                allele = pep$HLA_1[i],length = pep$pep_len[i],
                                pre_method = "netmhcpan_ba",client_path="~/software/mhc_i/src/",
                                tmp_dir=tempdir()),silent = TRUE)
    t2 <- try(mhcIbinding_client(peptide = pep$`WT peptide`[i],
                                 allele = pep$HLA_1[i],length = pep$pep_len[i],
                                 pre_method = "netmhcpan_ba",client_path="~/software/mhc_i/src/",
                                 tmp_dir=tempdir()),silent = TRUE)
    if ('try-error' %in% class(t)){
      next
    }else{
      pep$ic50_mt[i] <- t$ic50
      pep$rank_mt[i] <- t$rank
      pep$ic50_wt[i] <- t2$ic50
      pep$rank_wt[i] <- t2$rank
    }
  }else{
    t <- try(mhcIIbinding_client(peptide = pep$`Mut Epitope`[i],
                                 allele = pep$HLA_1[i],length = pep$pep_len[i],
                                 pre_method = "netmhciipan_ba",client_path="~/software/mhc_ii/",
                                 tmp_dir=tempdir()),silent = TRUE)
    t2 <- try(mhcIIbinding_client(peptide = pep$`WT peptide`[i],
                                  allele = pep$HLA_1[i],length = pep$pep_len[i],
                                  pre_method = "netmhciipan_ba",client_path="~/software/mhc_ii/",
                                  tmp_dir=tempdir()),silent = TRUE)
    if ('try-error' %in% class(t)){
      next
    }else{
      pep$ic50_mt[i] <- t$ic50
      pep$rank_mt[i] <- t$percentile_rank
      pep$ic50_wt[i] <- t2$ic50
      pep$rank_wt[i] <- t2$percentile_rank
    }
  }
  
}
saveRDS(pep,file = "data/neodb_all_pre.rds")
library(ggpubr)
pep <- readRDS("~/MHCbindshiny/data/neodb_all_pre.rds")
dt <- pep %>% 
  tidyr::pivot_longer(cols = c("ic50_wt","ic50_mt"),names_to = "ic50_type",values_to = "ic50")
p3 <- ggplot(data=dt,aes(x=ic50_type,y=log(ic50+1)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = c("MT","WT"))+
  labs(y="Log(IC50 + 1)")

dt <- pep %>% 
  tidyr::pivot_longer(cols = c("rank_wt","rank_mt"),names_to = "rank_type",values_to = "rank")
p4 <- ggplot(data=dt,aes(x=rank_type,y=log(rank+1)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = c("MT","WT"))+
  labs(y="Log(Rank + 1)")
library(patchwork)
p3 + p4 

##饼图
plot_pie <- function(vector, expression, label,main){
  per <- c(mean(eval(parse(text = expression))), 1 - mean(eval(parse(text = expression))))
  a <- pie(per, labels = paste0(label, " (", round(per, digits = 3) * 
                                  100, "%)"), border = "white", col = c("#66C2A5", "#FC8D62"),main = main)
  return(a)
}
a <- hla$super_type
plot_pie(a,"a=='HLA-I'",label = c("HLA-I","HLA-II"),main = "Total 102 HLA alleles")


###blast
library(seqinr)
neodb_all <- readRDS("~/MHCbindshiny/data/neodb_all.rds")
dt <- neodb_all %>% 
  select(pmid,`Mut Epitope`,`aa change`,HLA,`WT peptide`,gene) %>% 
  distinct_all() %>% 
  mutate(name=paste(pmid,gene,sep = "_"))
dt$name <- gsub("​","",dt$name)
dt$name <- gsub(" ","",dt$name)
aa <- as.list(dt$`Mut Epitope`)
names(aa) <- dt$name
write.fasta(aa,names = names(aa),file.out = "data/neodb.fasta")

library(Biostrings)
data(BLOSUM50)
BLOSUM50
seq1 <- "TRAATGRMV"
seq2 <- dt$`Mut Epitope`
scores <- vector("numeric",length = length(seq2))
for (i in 1:length(scores)){
  globalAligns <- pairwiseAlignment(seq1, seq2[i], substitutionMatrix = "BLOSUM50", gapOpening = -2,
                                    gapExtension = -8, scoreOnly = FALSE)
  scores[i] <- globalAligns@score
}
globalAligns <- pairwiseAlignment(seq1, seq2[which.max(scores)], substitutionMatrix = "BLOSUM50", gapOpening = -2,
                                  gapExtension = -8, scoreOnly = FALSE)

globalAligns
library(ggmsa)
library(Biostrings)
aa <- list(`Input seq`=globalAligns@pattern,`Subject seq`=globalAligns@subject)
seqinr::write.fasta(aa,names = names(aa),file.out = "data/test_msa.fasta")
readAAStringSet("data/test_msa.fasta") -> test_msa
ggmsa(test_msa, char_width = 0.5, seq_name = T)+
  labs(title=paste0("Score = ",globalAligns@score))

##
neodb_all <- readRDS("data/neodb_all.rds")
neodb_all <- neodb_all %>% 
  mutate(hla_type = case_when(
    grepl("^A",HLA) ~ "A",
    grepl("^B",HLA) ~ "B",
    grepl("^C",HLA) ~ "C",
    grepl("^DP",HLA) ~ "DP",
    grepl("^DQ",HLA) ~ "DQ",
    grepl("^DR",HLA) ~ "DR"
  )) %>% 
  mutate(super_type=ifelse(hla_type %in% c("A","B","C"),"HLA-I","HLA-II"))
saveRDS(neodb_all,file = "data/neodb_all.rds")

xlsx::write.xlsx(neodb_all,file = "filter_db.xlsx",showNA = FALSE)

##pie plot
install.packages("ggpie")
library(ggpie)
hlas <- neodb_all %>% select(HLA,hla_type,super_type) %>% distinct_all()

ggdonut(data = hlas, group_key = "hla_type", count_type = "full",
        label_info = "all", label_type = "horizon", label_split = NULL,
        label_size = 4, label_pos = "in")
##画个 seqlog
pep_9 <- neodb_all %>% filter(pep_len==9)

wt_pep <- pep_9 %>% filter(!is.na(`WT peptide`))
install.packages("ggseqlogo")
library(ggseqlogo)
data(ggseqlogo_sample)
ggseqlogo( pep_9$`Mut Epitope`, seq_type='aa')
ggseqlogo( wt_pep$`WT peptide`, seq_type='aa')
