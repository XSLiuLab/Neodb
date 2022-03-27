install.packages("easyPubMed")
devtools::install_github("dami82/easyPubMed")
library(easyPubMed)
new_query <- 'neoantigen' 
out.A <- batch_pubmed_download(pubmed_query_string = new_query, 
                               format = "txt", 
                               batch_size = 20,dest_dir = "data/",
                               dest_file_prefix = "easyPM_example")
