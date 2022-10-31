####
library(dplyr)
iedb <- read.csv("~/NeoTCR/data/tcell_table_export_1655876131.csv",header = T,skip = 1)
iedb <- iedb %>%
  select(Description,Method.Technique,Allele.Name,Qualitative.Measure,
         Number.of.Subjects.Tested,Number.of.Subjects.Responded) %>% 
  filter(grepl(":",Allele.Name)) %>% 
  filter(!grepl("mutant",Allele.Name)) %>% 
  filter(grepl("HLA",Allele.Name)) %>% 
  mutate(HLA=Allele.Name) %>% 
  filter(!is.na(Number.of.Subjects.Tested)) %>% 
  filter(!is.na(Number.of.Subjects.Responded))
iedb <- iedb %>% 
  distinct_all() %>% 
  rename(type=Qualitative.Measure) %>% 
  mutate(index=paste(Description,HLA,type,sep = "-"))
iedb <- iedb %>% 
  group_by(index) %>% 
  mutate(test=sum(Number.of.Subjects.Tested),response=sum(Number.of.Subjects.Responded)) %>% 
  ungroup() %>% 
  distinct(index,.keep_all = TRUE)
iedb <- iedb %>% 
  mutate(tt=ifelse((type =="Negative" & test >=4 | 
                      grepl("Positive",type)),"yes","no")) %>% 
  filter(tt=="yes") %>% select(-tt)
##把 neg 和 pos 中都有的肽删除
neg <- iedb %>% filter(type == "Negative") %>% 
  mutate(index = paste(Description,HLA,sep = "-"))
pos <- iedb %>% filter(type != "Negative") %>% 
  mutate(index = paste(Description,HLA,sep = "-"))
pos_filter <- pos %>% filter(!(index %in% neg$index))
neg_filter <- neg %>% filter(!(index %in% pos$index))

iedb <- bind_rows(pos_filter,neg_filter)
iedb <- iedb %>% select(Description,Method.Technique,Allele.Name,type,HLA,test,response)
colnames(iedb) <- c("pep","method","allele","immunogenicity","hla","test","respond")
iedb <- iedb %>% 
  select(pep,immunogenicity,hla,test,respond) %>% distinct_all()
all_dt <- iedb
all_dt <- all_dt %>% distinct_all()

pseudu  <- data.table::fread("~/software/neo_pre/netMHCpan-4.1/data/MHC_pseudo.dat",header = F,fill=TRUE,data.table = F)
pseudu <- pseudu %>% filter(grepl("HLA",V1)) %>% filter(grepl(":",V1))

all_dt <- all_dt %>% 
  mutate(hla=gsub("[*]","",hla)) %>% 
  filter(hla %in% pseudu$V1) %>% 
  rowwise() %>% 
  mutate(hla_seq=pseudu$V2[which(pseudu$V1==hla)]) %>% 
  ungroup()
aa <- c('A',
        'R',
        'N',
        'D',
        'C',
        'Q',
        'E',
        'G',
        'H',
        'I',
        'L',
        'K',
        'M',
        'F',
        'P',
        'S',
        'T',
        'W',
        'Y',
        'V')
all_dt <- all_dt %>% 
  rowwise() %>% 
  mutate(isinaa1=ifelse(length(which(!((strsplit(pep,"") %>% unlist()) %in% aa)))==0,"yes","no"),
         isinaa2=ifelse(length(which(!((strsplit(hla_seq,"") %>% unlist()) %in% aa)))==0,"yes","no")) %>% 
  ungroup()
all_dt <- all_dt %>% 
  filter(isinaa1 == "yes" & isinaa2 == "yes") %>% 
  select(-isinaa1,-isinaa2)
all_dt <- all_dt %>% 
  mutate(type = ifelse(immunogenicity == "Negative",0,1))
write.csv(all_dt,file = "data/train_data_iedb_2.csv",quote = F,row.names = F)

###验证数据
###在我们的数据库中进行验证
neodb_all <- readRDS("~/MHCbindshiny/data/neodb_all.rds")
val_dt <- neodb_all %>% 
  filter(super_type == "HLA-I") %>% 
  select(`Mut Epitope`,HLA) %>% 
  distinct_all()
all_dt <- read.csv("data/train_data_iedb_2.csv")
pseudu  <- data.table::fread("~/software/neo_pre/netMHCpan-4.1/data/MHC_pseudo.dat",header = F,fill=TRUE,data.table = F)
pseudu <- pseudu %>% filter(grepl("HLA",V1)) %>% filter(grepl(":",V1))
pseudu <- pseudu %>% filter(V2 %in% all_dt$hla_seq)

val_dt <- val_dt %>% filter(paste0("HLA-",HLA) %in% pseudu$V1)
val_dt <- val_dt %>% 
  filter(nchar(`Mut Epitope`)>=8 & nchar(`Mut Epitope`)<=11)
val_dt <- left_join(
  val_dt %>% mutate(HLA=paste0("HLA-",HLA)),
  pseudu %>% rename(HLA=V1)
)
colnames(val_dt)[1] <- "pep"
aa <- c('A',
        'R',
        'N',
        'D',
        'C',
        'Q',
        'E',
        'G',
        'H',
        'I',
        'L',
        'K',
        'M',
        'F',
        'P',
        'S',
        'T',
        'W',
        'Y',
        'V')
val_dt <- val_dt %>% 
  rowwise() %>% 
  mutate(isinaa1=ifelse(length(which(!((strsplit(pep,"") %>% unlist()) %in% aa)))==0,"yes","no"),
         isinaa2=ifelse(length(which(!((strsplit(V2,"") %>% unlist()) %in% aa)))==0,"yes","no")) %>% 
  ungroup()
val_dt <- val_dt %>% 
  filter(isinaa1 == "yes" & isinaa2 == "yes") %>% 
  select(-isinaa1,-isinaa2)
val_dt <- val_dt %>% rename(hla_seq=V2) %>%
  mutate(type=1) %>%
  select(pep,HLA,type,hla_seq)
val_dt <- val_dt %>% mutate(index=paste(pep,HLA,sep = "-"))

##add tesla neg
tesla_dt <- readRDS("~/NeoTCR/data/neo_model/tesla_dt.rds")
val_dt <- bind_rows(val_dt,
                    tesla_dt %>% mutate(index=paste(pep,HLA,sep = "-"))) %>% 
  distinct_all()
pos <- val_dt %>% filter(type==1)
neg <- val_dt %>% filter(type==0)
a <- intersect(neg$index,pos$index)
val_dt <- bind_rows(
  pos %>% filter(!(index %in% a)),
  neg %>% filter(!(index %in% a))
)
val_dt <- val_dt %>% 
  filter(nchar(pep)>=8 & nchar(pep)<=11)

train_dt <- read.csv("data/train_data_iedb_2.csv") %>% 
  mutate(index=paste(hla_seq,pep,sep = "-"))

val_dt <- val_dt %>% 
  mutate(index=paste(hla_seq,pep,sep = "-")) %>% 
  filter(hla_seq %in% train_dt$hla_seq) %>% 
  filter(!(index %in% train_dt$index))
table(val_dt$type)
saveRDS(val_dt,file = "data/val_dt.rds")
write.csv(val_dt, file = "data/db_pre_dt.csv",row.names = F,quote = F)
###预测
val_dt_res <- read.csv("data/db_pred_gnn.csv")
val_dt_res <- bind_cols(
  val_dt,val_dt_res %>% select(all_preds)
)
val_dt_res <- val_dt_res %>% select(-index)
saveRDS(val_dt_res,file = "data/own_immuno_res.rds")
###其他方法,不需要去掉训练数据
neodb_all <- readRDS("~/MHCbindshiny/data/neodb_all.rds")
val_dt <- neodb_all %>% 
  filter(super_type == "HLA-I") %>% 
  select(`Mut Epitope`,HLA) %>% 
  distinct_all()
colnames(val_dt)[1] <- "pep"
aa <- c('A',
        'R',
        'N',
        'D',
        'C',
        'Q',
        'E',
        'G',
        'H',
        'I',
        'L',
        'K',
        'M',
        'F',
        'P',
        'S',
        'T',
        'W',
        'Y',
        'V')
val_dt <- val_dt %>% 
  rowwise() %>% 
  mutate(isinaa1=ifelse(length(which(!((strsplit(pep,"") %>% unlist()) %in% aa)))==0,"yes","no")) %>% 
  ungroup()
val_dt <- val_dt %>% 
  filter(isinaa1 == "yes") %>% 
  select(-isinaa1)
val_dt <- val_dt %>% 
  mutate(type=1) %>% 
  select(pep,HLA,type) 
val_dt <- val_dt %>% mutate(index=paste(pep,paste0("HLA-",HLA),sep = "-"))

tesla_dt <- readRDS("~/NeoTCR/data/neo_model/tesla_dt.rds")
val_dt <- bind_rows(val_dt %>% select(-index) %>% mutate(HLA=paste0("HLA-",HLA)),
                    tesla_dt %>% select(-hla_seq)) %>% 
  distinct_all() %>%
  mutate(index=paste(HLA,pep,sep = "-"))
pos <- val_dt %>% filter(type==1)
neg <- val_dt %>% filter(type==0)
a <- intersect(neg$index,pos$index)
val_dt <- bind_rows(
  pos %>% filter(!(index %in% a)),
  neg %>% filter(!(index %in% a))
)

saveRDS(val_dt,file = "data/val_dt_other.rds")
##IEDB 
pseudu  <- data.table::fread("~/software/neo_pre/netMHCpan-4.1/data/MHC_pseudo.dat",header = F,fill=TRUE,data.table = F)
pseudu <- pseudu %>% filter(grepl("HLA",V1)) %>% filter(grepl(":",V1))
val_dt <- readRDS("data/val_dt_other.rds")
###remove train data
train_dt <- read.csv("data/train_data_iedb_2.csv") %>%
  left_join(.,pseudu %>% rename(hla_seq=V2)) %>% 
  mutate(index=paste(V1,pep,sep = "-"))
val_dt <- val_dt %>% 
  filter(!(index %in% train_dt$index))
iedb <- val_dt %>% 
  filter(nchar(pep) == 9)
iedb$HLA <- gsub("[:]","",iedb$HLA)
iedb <- iedb %>% 
  filter(HLA %in% available_alleles("Immuno","I","IEDB"))
iedb$immuno <- NA
for (i in 1:nrow(iedb)){
  tt <- mhcbinding_client(
    peptide = iedb$pep[i],
    allele = iedb$HLA[i],
    length=nchar(iedb$pep[i]),
    pre_method="IEDB",
    tmp_dir=tempdir(),
    hla_type="I",
    method_type="Immuno",
    Immuno_IEDB_path="~/software/immunogenicity/"
  )
  iedb$immuno[i] <- tt$score
}
saveRDS(iedb,file = "data/iedb_immuno_res.rds")
##prime
val_dt <- readRDS("data/val_dt_other.rds")
prime2 <- val_dt 
prime2$HLA <- gsub("HLA-","",prime2$HLA) %>% gsub("[:]","",.)
prime2 <- prime2 %>%  
  filter(nchar(pep)>=8 & nchar(pep)<=14) %>% 
  filter(HLA %in% available_alleles("Immuno","I","PRIME2.0"))
prime_train <- readxl::read_xlsx("~/NeoTCR/data/neo_model/mmc2.xlsx") %>% 
  mutate(index =paste(Allele,Mutant,sep = "-"))
prime2 <- prime2 %>% 
  mutate(index = paste(HLA,pep,sep="-")) %>% 
  filter(!(index %in% prime_train$index))
prime2$immuno <- NA
for (i in 1:nrow(prime2)){
  tt <- mhcbinding_client(
    peptide = prime2$pep[i],
    allele = prime2$HLA[i],
    length=nchar(prime2$pep[i]),
    pre_method="PRIME2.0",
    tmp_dir=tempdir(),
    hla_type="I",
    method_type="Immuno",
    PRIME_path="~/software/PRIME-2.0/",
    MixMHCpred_path="~/software/MixMHCpred/MixMHCpred"
  )
  prime2$immuno[i] <- tt$PRIME_score
}

saveRDS(prime2,file = "data/prime2_immuno_res.rds")
###deepimmuno 需要用网页
deepimmuno <- val_dt %>%
  rowwise() %>% 
  mutate(HLA=sub("(.{5})(.*)", "\\1*\\2", HLA)) %>% 
  filter(nchar(pep) %in% c(9,10)) %>% 
  mutate(HLA=gsub("[:]","",HLA)) %>% 
  mutate(index=paste(pep,HLA,sep = ":")) %>% 
  ungroup() %>% 
  filter(HLA %in% available_alleles("Immuno","I","DeepImmuno"))
deepimmu_train <- read.csv("~/NeoTCR/data/neo_model/deepimmun0_train.csv")
deepimmu_train <- deepimmu_train %>% mutate(index = paste(peptide,HLA,sep = "-"))
deepimmuno <- deepimmuno %>% filter(!(index %in% deepimmu_train$index))
write.csv(deepimmuno %>% select(pep,HLA),quote = F,row.names = F,file = "data/deepimmuno_input.csv")

deepimmuno_res <- read.table("data/deepimmuno_res.txt",header = TRUE)
deepimmuno_res$index <- paste(deepimmuno_res$peptide,deepimmuno_res$HLA,sep = ":")
deepimmuno_res <- left_join(
  deepimmuno,
  deepimmuno_res %>% select(index,immunogenicity)
) %>% ungroup() %>% filter(!is.na(immunogenicity))
saveRDS(deepimmuno_res,file = "data/deepimmuno_res.rds")
##seq2neo
seq2neo <- val_dt %>%
  filter(nchar(pep)>=8 & nchar(pep)<=11) %>% 
  rowwise() %>% 
  mutate(HLA=sub("(.{5})(.*)", "\\1*\\2", HLA)) %>% 
  ungroup() %>% 
  filter(HLA %in% available_alleles("Immuno","I","Seq2Neo-CNN"))
seq2neo_train <- read.table("data/add_need_info_total_data.tsv",sep = "\t",header = T) %>% 
  mutate(index=paste(pep,hla,sep = "-"))
seq2neo <- seq2neo %>% 
  mutate(index=paste(pep,HLA,sep = "-")) %>% 
  filter(!(index %in% seq2neo_train$index))
seq2neo$immuno <- NA
library(parallel)
library(doParallel)
library(foreach)
myCluster <- makeCluster(20,type = "PSOCK") # type of cluster
registerDoParallel(myCluster)
neo <- foreach(i=1:nrow(seq2neo),.export="seq2neo",.combine = "rbind") %dopar%  {
  library(dplyr)
  library(MHCbinding)
  tt <- mhcbinding_client(peptide = seq2neo$pep[i],
                          allele = seq2neo$HLA[i],
                          length=nchar(seq2neo$pep[i]),
                          pre_method="Seq2Neo-CNN",
                          tmp_dir=tempdir(),
                          hla_type="I",
                          method_type="Immuno",seq2neo_env = "DeepImmuno",
                          seq2neo_path = "~/software/Seq2Neo/seq2neo/function/immuno_Prediction/",
                          netchop_path = "~/software/netchop/",client_path = "~/software/mhc_i/src/")
  seq2neo$immuno[i] <- tt$immunogenicity
  tt <- seq2neo[i,]
  return(tt)
}
stopCluster(myCluster)
saveRDS(neo,file = "data/seq2neo_immuno_res.rds")
###比较
own_immuno_res <- readRDS("~/Neodb_model/data/own_immuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))
prime2_immuno_res <- readRDS("~/Neodb_model/data/prime2_immuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))
iedb_immuno_res <- readRDS("~/Neodb_model/data/iedb_immuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))
deepimmuno_res <- readRDS("~/Neodb_model/data/deepimmuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))
seq2neo_immuno_res <- readRDS("~/Neodb_model/data/seq2neo_immuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))

##sensitivity = recall 实际是新抗原中预测出来是新抗原的比例
re_own <- sum(own_immuno_res$type==1 & own_immuno_res$all_preds>0.5)/sum(own_immuno_res$type==1)
re_deepimm <- sum(deepimmuno_res$type==1 & deepimmuno_res$immunogenicity>0.5)/sum(deepimmuno_res$type==1)
re_iedb <- sum(iedb_immuno_res$type==1 & iedb_immuno_res$immuno>0)/sum(iedb_immuno_res$type==1)
re_prime <- sum(prime2_immuno_res$type==1 & prime2_immuno_res$immuno<=0.5)/sum(prime2_immuno_res$type==1)
re_seq2neo <- sum(seq2neo_immuno_res$type==1 & seq2neo_immuno_res$immuno>0.5)/sum(seq2neo_immuno_res$type==1)
###Precision 预测的是新抗原中实际是新抗原的比例
pre_own <- sum(own_immuno_res$type==1 & own_immuno_res$all_preds>0.5)/sum(own_immuno_res$all_preds>0.5)
pre_deepimm <- sum(deepimmuno_res$type==1 & deepimmuno_res$immunogenicity>0.5)/sum(deepimmuno_res$immunogenicity>0.5)
pre_iedb <- sum(iedb_immuno_res$type==1 & iedb_immuno_res$immuno>0)/sum(iedb_immuno_res$immuno>0)
pre_prime <- sum(prime2_immuno_res$type==1 & prime2_immuno_res$immuno<=0.5)/sum(prime2_immuno_res$immuno<=0.5)
pre_seq2neo <- sum(seq2neo_immuno_res$type==1 & seq2neo_immuno_res$immuno>0.5)/sum(seq2neo_immuno_res$immuno>0.5)

f1_own <- 2* ((pre_own*re_own)/(pre_own+re_own))
f1_deepimm <- 2* ((pre_deepimm*re_deepimm)/(pre_deepimm+re_deepimm))
f1_iedb <- 2* ((pre_iedb*re_iedb)/(pre_iedb+re_iedb))
f1_prime <- 2* ((pre_prime*re_prime)/(pre_prime+re_prime))
f1_seq2neo <- 2* ((pre_seq2neo*re_seq2neo)/(pre_seq2neo+re_seq2neo))

##画图
dt <- data.frame(
  `re_Immuno-GNN` = re_own,
  re_Deepimmuno = re_deepimm,
  re_IEDB = re_iedb,
  re_PRIME2 = re_prime,
  re_Seq2Neo = re_seq2neo,
  `pre_Immuno-GNN` = pre_own,
  pre_Deepimmuno = pre_deepimm,
  pre_IEDB = pre_iedb,
  pre_PRIME2 = pre_prime,
  pre_Seq2Neo = pre_seq2neo,
  `f1_Immuno-GNN` = f1_own,
  f1_Deepimmuno = f1_deepimm,
  f1_IEDB = f1_iedb, 
  f1_PRIME2 = f1_prime,
  f1_Seq2Neo = f1_seq2neo,check.names = F
) %>% t() %>% as.data.frame()
dt$type <- rownames(dt)
colnames(dt)[1] <- "score"
dt <- dt %>% tidyr::separate(col = type, into = c("metric","method"),sep = "_")
dt <- dt %>% 
  mutate(Metric=case_when(
    metric == "re" ~ "Sensitivity",
    metric == "pre" ~ "Precision",
    metric == "f1" ~ "F1"
  )) %>% filter(Metric != "Precision")
dt$score <- round(dt$score,3)
dt$method <- factor(dt$method,levels = c("Immuno-GNN","Seq2Neo",
                                         "Deepimmuno","PRIME2","IEDB"))
p1 <- ggbarplot(dt, "Metric", "score",
                fill = "method", color = "method", palette = "Paired",
                label = TRUE,
                position = position_dodge(0.9))
###topk
top_mat <- data.frame(method=c("Immuno-GNN","Seq2Neo","Deepimmuno","IEDB","PRIME2"),
                      `Top 20`=NA,
                      `Top 50`=NA,check.names = FALSE)
get_top <- function(own,deepimme,cnn,iedb,prime,top_num){
  tt <- own_immuno_res %>% slice_max(order_by = all_preds,n=top_num) %>% `[`(1:top_num,)
  tt1 <- seq2neo_immuno_res %>% slice_max(order_by = immuno,n=top_num) %>% `[`(1:top_num,) 
  tt2 <- deepimmuno_res %>% slice_max(order_by = immunogenicity,n=top_num) %>% `[`(1:top_num,)
  tt3 <- iedb_immuno_res %>% slice_max(order_by = immuno,n=top_num) %>% `[`(1:top_num,)
  tt4 <- prime2_immuno_res %>% slice_min(order_by = immuno,n=top_num) %>% `[`(1:top_num,)
  return(c(sum(tt$type==1),sum(tt1$type==1),sum(tt2$type==1),sum(tt3$type==1),sum(tt4$type==1)))
}

top_mat$`Top 20` <- get_top(own_immuno_res,
                            deepimmuno_res,
                            seq2neo_immuno_res,
                            iedb_immuno_res,
                            prime2_immuno_res,top_num = 20)
top_mat$`Top 50` <- get_top(own_immuno_res,
                            deepimmuno_res,
                            seq2neo_immuno_res,
                            iedb_immuno_res,
                            prime2_immuno_res,top_num = 50)

top_mat <- top_mat %>% 
  tidyr::pivot_longer(cols = c("Top 20","Top 50"),values_to = "counts",
                      names_to = "Top_n")
top_mat$method <- factor(top_mat$method,
                         levels = c("Immuno-GNN","Seq2Neo",
                                    "Deepimmuno","PRIME2","IEDB"))
p2 <- ggbarplot(top_mat, "Top_n", "counts",
                fill = "method", color = "method", palette = "Paired",
                label = TRUE,
                position = position_dodge(0.9))
library(patchwork)
p1 + p2 + plot_layout(guides = 'collect') & theme(legend.position='top')
##roc
library(plotROC)
p1 <- ggplot(own_immuno_res, aes(d = type, m = all_preds)) + 
  geom_roc(labels=FALSE,pointsize=0)
calc_auc(p1)$AUC
p1 <- ggplot(deepimmuno_res, aes(d = type, m = immunogenicity)) + 
  geom_roc(labels=FALSE,pointsize=0)
calc_auc(p1)$AUC
p1 <- ggplot(seq2neo_immuno_res, aes(d = type, m = immuno)) + 
  geom_roc(labels=FALSE,pointsize=0)
calc_auc(p1)$AUC

###ROC
library(dplyr)
library(ggprism)
files <- list.files("data/models/",pattern = "train",full.names = T)
train_dt <- lapply(files, function(x){
  dt <- read.csv(x) %>% select(-X)
  dt$fold <- gsub("data/models//gnn_train_res_","",x) %>% gsub(".csv","",.)
  dt
}) %>% bind_rows()
train_dt$fold <- paste0("Fold ",as.numeric(train_dt$fold)+1)
library(plotROC)
p1 <- ggplot(train_dt, aes(d = all_labels, m = all_preds_raw, color = fold)) + 
  geom_roc(labels=FALSE,pointsize=0)
train_auc <- calc_auc(p1)$AUC %>% round(.,2)
train_dt <- train_dt %>% 
  mutate(fold = case_when(
    fold == "Fold 1" ~ paste0("Fold 1: AUC=",train_auc[1]),
    fold == "Fold 2" ~ paste0("Fold 2: AUC=",train_auc[2]),
    fold == "Fold 3" ~ paste0("Fold 3: AUC=",train_auc[3]),
    fold == "Fold 4" ~ paste0("Fold 4: AUC=",train_auc[4]),
    fold == "Fold 5" ~ paste0("Fold 5: AUC=",train_auc[5]),
    fold == "Fold 6" ~ paste0("Fold 6: AUC=",train_auc[6]),
    fold == "Fold 7" ~ paste0("Fold 7: AUC=",train_auc[7]),
    fold == "Fold 8" ~ paste0("Fold 8: AUC=",train_auc[8]),
    fold == "Fold 9" ~ paste0("Fold 9: AUC=",train_auc[9]),
    fold == "Fold 10" ~ paste0("Fold 10: AUC=",train_auc[10])
  ))
train_dt$fold <- factor(train_dt$fold,levels = unique(train_dt$fold))
p1 <- ggplot(train_dt, aes(d = all_labels, m = all_preds_raw, color = fold)) + 
  geom_roc(labels=FALSE,pointsize=0) +
  theme_prism()+
  labs(x="FPR",y="TPR",title="Cross Validation-Train")+
  theme(legend.position = c(0.7, 0.4))
p1
files <- list.files("data/models/",pattern = "test",full.names = T)
test_dt <- lapply(files, function(x){
  dt <- read.csv(x) %>% select(-X)
  dt$fold <- gsub("data/models//gnn_test_res_","",x) %>% gsub(".csv","",.)
  dt
}) %>% bind_rows()
test_dt$fold <- paste0("Fold ",as.numeric(test_dt$fold)+1)
p2 <- ggplot(test_dt, aes(d = all_labels, m = all_preds_raw, color = fold)) + 
  geom_roc(labels=FALSE,pointsize=0)
test_auc <- calc_auc(p2)$AUC %>% round(.,2)
test_dt <- test_dt %>% 
  mutate(fold = case_when(
    fold == "Fold 1" ~ paste0("Fold 1: AUC=",test_auc[1]),
    fold == "Fold 2" ~ paste0("Fold 2: AUC=",test_auc[2]),
    fold == "Fold 3" ~ paste0("Fold 3: AUC=",test_auc[3]),
    fold == "Fold 4" ~ paste0("Fold 4: AUC=",test_auc[4]),
    fold == "Fold 5" ~ paste0("Fold 5: AUC=",test_auc[5]),
    fold == "Fold 6" ~ paste0("Fold 6: AUC=",test_auc[6]),
    fold == "Fold 7" ~ paste0("Fold 7: AUC=",test_auc[7]),
    fold == "Fold 8" ~ paste0("Fold 8: AUC=",test_auc[8]),
    fold == "Fold 9" ~ paste0("Fold 9: AUC=",test_auc[9]),
    fold == "Fold 10" ~ paste0("Fold 10: AUC=",test_auc[10])
  ))
test_dt$fold <- factor(test_dt$fold,levels = unique(test_dt$fold))
p2 <- ggplot(test_dt, aes(d = all_labels, m = all_preds_raw, color = fold)) + 
  geom_roc(labels=FALSE,pointsize=0) +
  theme_prism()+
  labs(x="FPR",y="TPR")+
  labs(x="FPR",y="TPR",title="Cross Validation-Test")+
  theme(legend.position = c(0.7, 0.4))
p2
library(patchwork)
p1 + p2


###save
own_immuno_res <- readRDS("~/Neodb_model/data/own_immuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))
prime2_immuno_res <- readRDS("~/Neodb_model/data/prime2_immuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))
iedb_immuno_res <- readRDS("~/Neodb_model/data/iedb_immuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))
deepimmuno_res <- readRDS("~/Neodb_model/data/deepimmuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))
seq2neo_immuno_res <- readRDS("~/Neodb_model/data/seq2neo_immuno_res.rds") %>% 
  filter(nchar(pep) %in% c(9,10))

xlsx::write.xlsx(own_immuno_res,file = "~/Neodb_model/data/val_res.xlsx",
                 sheetName = "Immuno-GNN")
xlsx::write.xlsx(prime2_immuno_res,file = "~/Neodb_model/data/val_res.xlsx",
                 sheetName = "PRIME2",append = TRUE)
xlsx::write.xlsx(iedb_immuno_res,file = "~/Neodb_model/data/val_res.xlsx",
                 sheetName = "IEDB",append = TRUE)
xlsx::write.xlsx(deepimmuno_res,file = "~/Neodb_model/data/val_res.xlsx",
                 sheetName = "Deepimmuno",append = TRUE)
xlsx::write.xlsx(seq2neo_immuno_res,file = "~/Neodb_model/data/val_res.xlsx",
                 sheetName = "Seq2Neo",append = TRUE)
