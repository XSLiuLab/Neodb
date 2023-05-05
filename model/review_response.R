library(dplyr)
train_dt <- read.csv("/home/data/sda/wt/Neodb_model/raw/train_data_iedb_2.csv")
all_weight <- data.table::fread("model/all_weight_mean.csv",data.table = F)
all_edge <- data.table::fread("model/all_edge_index.csv",data.table = F)
train_dt <- train_dt %>% 
  mutate(pep_len = nchar(pep),
         hla_len = nchar(hla_seq))
###交叉合并
pep_hla_len <- as.vector(matrix(c(train_dt$pep_len, train_dt$hla_len), 
                                nrow = 2, byrow = TRUE))
pep_hla_index <- c(0,cumsum(pep_hla_len))

index <- which(train_dt$pep_len == 9)
pep_index <- seq(1,length(pep_hla_index),by=2)[index]
hla_index <- seq(2,length(pep_hla_index),by=2)[index]

##得到所有 pep 氨基酸的编号
pep_res <- lapply(pep_index[-length(pep_index)],function(x){
  seq(pep_hla_index[x],pep_hla_index[x + 1]-1)
})
pep_res <- unlist(pep_res)

###所有 HLA 氨基酸的编号
hla_res <- lapply(hla_index,function(x){
  seq(pep_hla_index[x],pep_hla_index[x + 1]-1)
})
hla_res <- unlist(hla_res)

all_edge_need <- all_edge %>% 
  filter((index1 %in% pep_res) & (index2 %in% hla_res))
all_weight_need <- all_weight %>% filter(V1 %in% all_edge_need$V1)

###
all_weight_need$pep_index <- rep(c(1:9), each=34)
all_weight_need$hla_index <- rep(c(1:34), nrow(all_weight_need) / 34)
all_weight_need$item_index <- rep(c(1:(nrow(all_weight_need)/(34*9))), 
                                  each = 34*9)

# all_weight_need$pep_index <- rep(c(1:10), each=34)
# all_weight_need$hla_index <- rep(c(1:34), nrow(all_edge_need) / 34)
# all_weight_need$item_index <- rep(c(1:2281), each = 34*10)

all_weight_need_summ <- all_weight_need %>% 
  group_by(pep_index,hla_index) %>% 
  summarise(median_weight = median(weight_mean)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(names_from = hla_index, values_from = median_weight) %>% as.data.frame()
rownames(all_weight_need_summ) <- all_weight_need_summ$pep_index
all_weight_need_summ <- all_weight_need_summ %>% select(-pep_index) %>% as.matrix()
all_weight_need_summ <- all_weight_need_summ * 100

pos <- apply(all_weight_need_summ,1,sum) %>% as.data.frame()
colnames(pos) <- "weight_sum"
pos$pos <- rownames(pos)

library(ComplexHeatmap)
pdf(file="model/weight_heatmap.pdf",width=16,height = 4)
row_ha <- rowAnnotation(foo = anno_text(as.character(round(pos$weight_sum)), 
                                        location = 0.5, just = "center",
                                        gp = gpar(border = "black"),
                                        width = max_text_width(as.character(round(pos$weight_sum)))*1.2))
ht <- Heatmap(all_weight_need_summ,
              cluster_rows = F,cluster_columns = F,
              heatmap_legend_param = list(title = "weight"),
              right_annotation = row_ha)
ComplexHeatmap::draw(ht)
dev.off()

library(ggpubr)

ggbarplot(data = pos, x = "pos", y = "weight_sum")


