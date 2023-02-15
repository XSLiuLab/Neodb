library(dplyr)

example_driver <- readRDS("~/tmp/mhcbinding/example_driver.rds")
example_driver <- example_driver %>% 
  select(-pos_alter,-transcript,-exon,-index) %>% 
  select(trans,everything()) %>% 
  group_by(across(2:aa)) %>% 
  summarise(contain_trans = paste(sort(trans),collapse  = ","))
##文件太大，先把 gene 及其 突变 氨基酸改变 保存下来
anchar <- example_driver %>% 
  ungroup() %>% 
  select(genes,aa) %>%
  distinct_all(.keep_all = T)
saveRDS(anchar,file = "data/example_anchar.rds")
example_driver <- example_driver %>% ungroup()
###分突变保存
gene <- unique(anchar$genes)
for (i in gene){
  dir.create(paste0("data/test/",i))
  dt <- anchar %>% 
    filter(genes == i)
  for (j in dt$aa){
    dt2 <- example_driver %>% 
      filter(genes == i & aa == j)
    message(paste0(i," ",j,": ",paste(unique(dt2$contain_trans),collapse = ":"),"\n"))
    saveRDS(dt2,file = paste0("data/test/",i,"/",j,".rds"))
  }
}
###绘图
library(latex2exp)
p.E224D <- readRDS("~/MHCbindshiny/data/test/TP53/p.E224D.rds")
p.E224D <- p.E224D %>% 
  arrange(length)
#p.E224D$ic50 <- 1 -(log(p.E224D$ic50)/log(50000))
p.E224D$ic50 <- ifelse(p.E224D$ic50 < 500,1,0)
dt <- p.E224D %>% 
  select(length,allele,peptide,ic50) %>% 
  filter(grepl("HLA-A",allele))
library(ComplexHeatmap)
library(circlize)
res <- vector("list",7)
names(res) <- unique(dt$length)
for (i in seq_along(res)){
  df <- dt %>% 
    filter(length == names(res)[i]) %>% 
    select(-length) %>% 
    tidyr::pivot_wider(names_from = "allele",values_from = "ic50")
  df <- as.data.frame(df)
  rownames(df) <- df$peptide
  df <- df %>% select(-peptide)
  col_fun = colorRamp2(c(0, max(dt$ic50)), c("white", "blue"))
  if(i ==1 ){
    h1 <- Heatmap(df,col = col_fun,cluster_rows = F,cluster_columns = F,
                  rect_gp = gpar(col = "grey", lwd = 2),row_names_side = "left",
                  show_heatmap_legend=T,row_names_gp = gpar(fontsize = 8))
  }else{
    h1 <- Heatmap(df,col = col_fun,cluster_rows = F,cluster_columns = F,
                  rect_gp = gpar(col = "grey", lwd = 2),row_names_side = "left",
                  show_heatmap_legend=F,row_names_gp = gpar(fontsize = 8))
  }
  res[[i]] <- h1
}

ht_list1 <- res[[1]] %v%  res[[2]] %v%  res[[3]] %v%  res[[4]] %v%  res[[5]] %v%  res[[6]] %v%  res[[7]] 
draw(ht_list,merge_legends =T)

library(cowplot)
p1 <- grid::grid.grabExpr(draw(ht_list,merge_legends =T,heatmap_legend_side = "top"))
p2 <- grid::grid.grabExpr(draw(ht_list1,merge_legends =T,heatmap_legend_side = "top"))
plot_grid(p2,p1,ncol = 1,align = "v")
plot_grid(plotlist = list(ht_list,ht_list))
