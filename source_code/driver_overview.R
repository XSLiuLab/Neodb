
### 
genome_version <- "hg38"
temp_dir <- "~/tmp/mhcbinding/"
txt_path <- "data/TCGA_driver_mutations_annovar_input.txt"
annovar_path = "~/software/annovar/"
comm_annotate <- paste0(annovar_path,"/annotate_variation.pl -out ",temp_dir,"/test -build ",genome_version," ",txt_path,
                        " ",annovar_path,"/humandb/ -hgvs -dbtype ensGene")
system(comm_annotate)
comm_coding <- paste0(annovar_path,"/coding_change.pl ",temp_dir,"/test.exonic_variant_function ",annovar_path,"/humandb/",
                      genome_version,"_ensGene.txt ",annovar_path,"/humandb/",genome_version,"_ensGeneMrna.fa ",
                      "-includesnp -onlyAltering --alltranscript --outfile ",temp_dir,"/coding.txt")
system(comm_coding)

dt <- seqinr::read.fasta(file = paste0(temp_dir,"/coding.txt"),seqtype = "AA",as.string = TRUE,whole.header = T)
anno_mut <- data.table::fread("~/tmp/mhcbinding/test.exonic_variant_function")

### 一个基因可以有三个层面的变异：基因组的 hotspot + CDNA change + aa change
##先把每个转录本对应的基因 cDNA change 和 aa change 弄出来
anno_mut <- anno_mut %>% as.data.frame()  %>% 
  mutate(index=paste(V4,V5,V6,V7,V8,sep = ":")) %>% as.data.frame()
res <- vector("list",length = nrow(anno_mut))
for (i in seq_along(res)){
  trans <- sapply(anno_mut$V3[i] %>% strsplit(.,split = ",") %>% `[[`(.,1) %>% strsplit(.,split = ":"),`[[`,2)
  genes <- sapply(anno_mut$V3[i] %>% strsplit(.,split = ",") %>% `[[`(.,1) %>% strsplit(.,split = ":"),`[[`,1)
  exon <- sapply(anno_mut$V3[i] %>% strsplit(.,split = ",") %>% `[[`(.,1) %>% strsplit(.,split = ":"),`[[`,3)
  cdna <- sapply(anno_mut$V3[i] %>% strsplit(.,split = ",") %>% `[[`(.,1) %>% strsplit(.,split = ":"),`[[`,4)
  aa <- sapply(anno_mut$V3[i] %>% strsplit(.,split = ",") %>% `[[`(.,1) %>% strsplit(.,split = ":"),`[[`,5)
  df <- data.frame(trans,genes,exon,cdna,aa,index=anno_mut$index[i])
  res[[i]] <- df
}
res <- bind_rows(res)

###看每个基因的三个层面
res_summ <- res %>% 
  group_by(genes) %>% 
  summarise(G_mut = length(unique(index)),
            cdna_mut = length(unique(cdna)),
            aa_change = length(unique(aa))) %>% 
  arrange(desc(G_mut)) %>% filter(G_mut>15)
res_summ$genes <- factor(res_summ$genes,levels = res_summ$genes)
res_summ_t <- res_summ %>% 
  tidyr::pivot_longer(cols = G_mut:aa_change,names_to = "Type",values_to = "counts")

library(tidytext)
library(ggplot2)
library(ggprism)
ggplot(data=res_summ_t,aes(x=genes,y=log(counts)))+
  geom_bar(stat = "identity")+
  facet_wrap(~Type,ncol = 1)+
  theme_prism()+
  theme(axis.text.x = element_text(angle = 30))

####
drivers <- readRDS("../tmp/TCGA_drivers_all.rds")
res <- res %>% 
  mutate(trans = gsub("[.].+","",trans)) %>% 
  mutate(index=paste(index,trans,sep = ":"))
drivers <- drivers %>% 
  mutate(index=paste(chr,start,end,ref,alt,transcript,sep = ":"))
drivers <- left_join(
  drivers,
  res,by="index"
)
saveRDS(drivers,file = "../tmp/TCGA_drivers_all.rds")

##先拿一部分出来做 APP
test_gene <- unique(drivers$genes)[1:10]

test_dt <- drivers %>% filter(genes %in% test_gene) 
saveRDS(test_dt,file = "~/tmp/mhcbinding/example_driver.rds")




