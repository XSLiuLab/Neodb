library(shiny)
library(ggplot2)
library(purrr)
library(shinythemes)
library(shinyWidgets)
library(MHCbinding)
library(DT)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(shinyFeedback)
library(waiter)
library(heatmaply)
library(slickR)
library(shinymaterial)
library(Biostrings)
#dyn.load('/home/wt/miniconda3/lib/libproj.so.15')
#dyn.load("/home/wt/miniconda3/lib/libpng16.so.16")
#dyn.load("/home/wt/miniconda3/lib/libjpeg.so.9")
#dyn.load("/home/wt/miniconda3/lib/libharfbuzz-icu.so.0")
library(Cairo)
library(ggmsa)
library(seqinr)
library(rintrojs)
system(". ~/.bashrc")
#system("conda")
message(system("echo $USER"))
message("a")
message(Sys.getenv("PATH"))
Sys.setenv(PATH="/home/ubuntu/miniconda3/bin:/home/ubuntu/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin")
message(Sys.getenv("PATH"))
options(shiny.maxRequestSize=10*1024^2)
TP53 <- readRDS("data/driver/TP53/neo.rds")
neodb_all <- readRDS("data/neodb_all.rds")
neodb_all <- neodb_all %>% filter(Gene != "NA") %>% filter(Gene != "ERBB2IP<U+00A0") %>%
        filter(!(HLA %in% c("Class I","ClassII","Class II","Class I and II")))
candidate_dirver_neo <- readRDS("data/candidate_dirver_neo.rds")
#colnames(neodb_all)[1:10] <- c("PMID","Sample ID","Mut Epitope","AA change","HLA","WT","Gene","T cell","Assay",
#                               "Cancer type")
anchar <- readRDS("data/all_gene_aa_trans.rds")
gene <- unique(anchar$gene)
plot_heatmap <- function(dt,hla_type,need_value,legend_title,need_thre=FALSE,origin_dt=NULL,thre){
  if (need_thre){
    dt$ic50 <- ifelse(dt$ic50 < as.numeric(thre),1,0)
    dt$rank <- ifelse(dt$rank < as.numeric(thre),1,0)
  }else{
    dt$ic50 <- 1 -(log(dt$ic50)/log(50000))
    dt$rank <- 100 - dt$rank
  }
  dt <- dt %>% 
    select(length,allele,peptide,need_value) %>% 
    filter(grepl(paste0("HLA-",hla_type),allele))
  df <- dt %>% 
    select(-length) %>% 
    tidyr::pivot_wider(names_from = "allele",values_from = need_value)
  df <- as.data.frame(df)
  rownames(df) <- df$peptide
  df <- df %>% select(-peptide)
  
  origin_dt <- origin_dt %>% 
    select(length,allele,peptide,need_value) %>% 
    filter(grepl(paste0("HLA-",hla_type),allele))
  origin_dt <- origin_dt %>% 
    select(-length) %>% 
    tidyr::pivot_wider(names_from = "allele",values_from = need_value)
  origin_dt <- as.data.frame(origin_dt)
  rownames(origin_dt) <- origin_dt$peptide
  origin_dt <- origin_dt %>% select(-peptide)
  
  hide_colorbar <- FALSE
  if (need_thre){
    hide_colorbar <- TRUE
  }
  
  p <- heatmaply(
    df,grid_gap = 1,grid_color = "grey",Rowv=F,Colv=F,fontsize_row=6.5,fontsize_col = 7,
    key.title=legend_title,colorbar_len=0.2,
    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
      low = "white", 
      high = "red", 
      limits = c(0, max(df))
    ),label_names=c("Peptide","HLA",toupper(need_value)),custom_hovertext=origin_dt,hide_colorbar = hide_colorbar
  )
  return(p)
}

source('tools_tab.R', local = TRUE)
source('driver_tab.R', local = TRUE)
source('neodb_tab.R', local = TRUE)
source('our_method_tab.R', local = TRUE)
ui <- bootstrapPage(
  navbarPage(id = "landing", title= "Neodb", theme = shinytheme("flatly"), collapsible = TRUE,
             windowTitle = "Neodb",
             tags$head(includeCSS("styles.css")),
             tabPanel(
               title = "Home", icon = icon("home"),
               fluidRow(h3("Welcome to the Neodb", style="text-align: center;")),
               div(class="landing-page-box",
                   div(img(src="neodb-04.png",height="80%", width="80%"),style="text-align: center;"),
                   fluidRow(
                     align = "center",
                     br(),
                     column(
                       12,
                       splitLayout(cellWidths = c("20%","20%","20%","20%"),
                                   actionBttn(
                                     inputId = "landing1",
                                     label = "Prediction Tools",
                                     color = "primary",
                                     style = "gradient",
                                     block = FALSE,size="lg"
                                   ),
                                   actionBttn(
                                     inputId = "landing2",
                                     label = "Driver mutations",
                                     color = "primary",
                                     style = "gradient",
                                     block = FALSE,size="lg"
                                   ),
                                   actionBttn(
                                     inputId = "landing3",
                                     label = "Validate neoantigens",
                                     color = "primary",
                                     style = "gradient",
                                     block = FALSE,size="lg"
                                   ),
                                   actionBttn(
                                     inputId = "landing4",
                                     label = "Immunogenicity-GNN",
                                     color = "primary",
                                     style = "gradient",
                                     block = FALSE,size="lg"
                                   )
                       )
                     )
                   ) # fluidRow
               )),
             tab_tools,
             tab_driver,
             tab_neodb,
             immuno_gnn_tab,
             tabPanel("About", br(),
                      fluidRow(
                        align = "left",
                        br(),
                        column(
                          12,
                          splitLayout(cellWidths = c("20%","20%","20%","20%"),
                                      
                                      downloadBttn(
                                        outputId = "download_maf_exp",
                                        label = "Example maf file",
                                        color = "primary",
                                        style = "gradient",
                                        block = FALSE,size="md"
                                      ),
                                      downloadBttn(
                                        outputId = "download_vcf_exp",
                                        label = "Example vcf file",
                                        color = "primary",
                                        style = "gradient",
                                        block = FALSE,size="lg"
                                      ),
                                      downloadBttn(
                                        outputId = "download_txt_exp",
                                        label = "Example txt file",
                                        color = "primary",
                                        style = "gradient",
                                        block = FALSE,size="lg"
                                      ),
                                      downloadBttn(
                                        outputId = "download_pep_exp",
                                        label = "Example pep file",
                                        color = "primary",
                                        style = "gradient",
                                        block = FALSE,size="lg"
                                      ),
				      downloadBttn(
                                        outputId = "download_neodb_exp",
                                        label = "All Neodb data",
                                        color = "primary",
                                        style = "gradient",
                                        block = FALSE,size="lg"
                                      )
                          )
                        )
                      ),br(),
                      includeMarkdown("About.Rmd"))
             )
             
)

server <- function(input, output, session) {
  
  observeEvent(input$landing1, {
    showModal(modalDialog(easyClose = TRUE,
      h3("What type of tool do you want to use?"),
      fluidRow(
        align = "center",
        br(),
        column(
          12,
          splitLayout(cellWidths = c("20%","30%","30%"),
                      actionBttn(
                        inputId = "jump_to_binding",
                        label = "Binding",
                        color = "primary",
                        style = "fill",
                        block = FALSE,size="md"
                      ),
                      actionBttn(
                        inputId = "jump_to_processing",
                        label = "Processing",
                        color = "primary",
                        style = "fill",
                        block = FALSE,size="md"
                      ),
                      actionBttn(
                        inputId = "jump_to_immuno",
                        label = "Immunogenicity",
                        color = "primary",
                        style = "fill",
                        block = FALSE,size="md"
                      )
          )
        )
      )
    )
  )
  })
  
  # observeEvent(c(input$jump_to_binding,input$jump_to_processing,input$jump_to_immuno), {
  #   removeModal()
  # })
  
  observeEvent(input$jump_to_binding, {
    updateTabsetPanel(session, "landing", selected = "Binding")
  })
  observeEvent(input$jump_to_processing, {
    updateTabsetPanel(session, "landing", selected = "Processing")
  })
  observeEvent(input$jump_to_immuno, {
    updateTabsetPanel(session, "landing", selected = "Immuno")
  })
  
  observeEvent(input$landing2, {
    updateTabsetPanel(session, "landing", selected = "Driver search")
  })
  
  observeEvent(input$landing3, {
    updateTabsetPanel(session, "landing", selected = "Neodb")
  })
  
  observeEvent(input$landing4, {
    updateTabsetPanel(session, "landing", selected = "Immuno-GNN")
  })
  
  ##tools
  tools_server("Binding",module_id = "Binding")
  tools_server("Processing",module_id = "Processing")
  tools_server("Immuno",module_id = "Immuno")
  
  ##our method 
  immuno_gnn_server("Immuno-GNN")
  
  ##driver
  aa <- reactive({
    dt <- anchar %>% 
      filter(gene == input$gene_name)
    return(c("User input",unique(dt$aa)))
  })
  
  observeEvent(input$gene_name,{
    ##发生改变前把冻住
    freezeReactiveValue(input,"aa_change")
    updatePickerInput(session = session,inputId = "aa_change",
                      choices = aa(),
                      selected = aa()[1])
  })
  
  observeEvent(input$gene_name,{
    updateTextInput(session, "aa_change_text",
                    value = aa()[2])
  })
  
  observeEvent(input$aa_change_text,{
    currentaa <- anchar[anchar$gene == input$gene_name,"aa"]
    if (input$aa_change_text %in% currentaa$aa) {
      hideFeedback("aa_change_text")
    } else {
      showFeedbackWarning(
        inputId = "aa_change_text",
        text = "This mutation does not occour in the gene"
      )
    }
    #feedbackWarning("aa_change_text", !(input$aa_change_text %in% currentaa$aa), "This mutation does not occour in the gene")
  })
  
  observeEvent(input$input_thre,{
    range <- ifelse(input$ic50 == "IC50",Inf,100)
    if (is.na(as.numeric(input$input_thre))){
      NULL
    }else{
      if(as.numeric(input$input_thre) >0 & as.numeric(input$input_thre) < range){
        hideFeedback("input_thre")
      }else{
        showFeedbackWarning(
          inputId = "input_thre",
          text = "Please input valid range: > 0 for IC50, 1-100 for Rank"
        ) 
      }
    }
  })
  gene_name <- reactive({
    input$gene_name
  })
  aa_change <- reactive({
    ifelse(input$aa_change == "User input",
           input$aa_change_text,input$aa_change)
  })
  
  data <- eventReactive(input$search,{
    currentaa <- anchar[anchar$gene == input$gene_name,"aa"]
    req((aa_change() %in% currentaa$aa))
    if (input$gene_name == "TP53"){
      dt <- TP53
    }else {
      dt <- readRDS(paste0("data/driver/",gene_name(),"/neo.rds"))
    }
    dt <- dt %>% 
      filter(grepl(aa_change(),pos_alter))
    return(dt)
  })
  data_candidate <- eventReactive(input$search,{
    currentaa <- anchar[anchar$gene == input$gene_name,"aa"]
    req((aa_change() %in% currentaa$aa))
    dt <- candidate_dirver_neo %>% 
      filter(gene == input$gene_name) %>% 
      filter(grepl(aa_change(),pos_alter)) %>% 
      select(-cdna,-core,-icore)
    return(dt)
  })
  
  #observeEvent(c(input$gene_name,input$aa_change,input$aa_change_text),{
  #  
  #  #dt <- data() %>% as.data.frame()
  #  enst <- anchar %>% 
  #    filter(gene == gene_name()) %>% 
  #    filter(aa == aa_change())
  #  enst <- enst$transcript
  #  print(length(enst))
  #  enst <- base::strsplit(enst,split = ",")[[1]]
  #  updatePickerInput(session = session,inputId = "enst",
  #                    choices = enst,
  #                    selected = enst[1])
  #})
  
  # data_dt <- reactive({
  #   dt <- readRDS(paste0("data/test/",gene_name(),"/",aa_change(),".rds"))
  #   dt <- dt %>% filter(contain_trans == input$enst)
  #   return(dt)
  # })
  
  out_plot <- eventReactive(c(input$search,input$ok),ignoreInit = T,{
    gene_name <- reactive({
      input$gene_name
    })
    aa_change <- reactive({
      ifelse(input$aa_change == "User input",
             input$aa_change_text,input$aa_change)
    })
    if (gene_name() == "TP53"){
      dt <- TP53
    }else {
      dt <- readRDS(paste0("data/driver/",gene_name(),"/neo.rds"))
    }
    dt <- dt %>% 
      filter(grepl(aa_change(),pos_alter))
    dt <- dt %>% 
      arrange(length)
    thre <- input$input_thre
    if(is.null(thre)){thre <- 500}
    if(input$ic50 == "IC50"){
      P <- plot_heatmap(dt,hla_type = input$hla,need_value = "ic50",
                        legend_title = "1-log(IC50)/log(50000)",need_thre = input$need_thre,
                        origin_dt = dt,thre = thre)
    }else{
      P <- plot_heatmap(dt,hla_type = input$hla,need_value = "rank",
                        legend_title = "1 - %Rank",need_thre = input$need_thre,
                        origin_dt = dt,thre = thre)
    }
    return(P)
  })
  
  output$threshold <- renderUI({
    if (input$need_thre){
      textAreaInput("input_thre","Input the threshold for IC50 or %Rank:",value = "500")
    }else{
      tags$div("   ", style="display:inline-block")
    }
  })
  
  observeEvent(input$search,{
    output$download_menu_image <- renderUI({
      dropdownButton(
        tags$div(selectInput("download_image", "Download image: ",   
                             choices = c("PNG","PDF"),
                             selected = c("PDF"),
                             multiple = FALSE,width = "100%"), style="display:inline-block"),
        downloadBttn(
          outputId = "ok_download_image",
          label = "OK",
          color = "primary",
          style = "pill",icon = icon("ok",lib = "glyphicon"),size = "sm"),
        
        #hr(style = "border-top: 1px solid #000000;"),
        circle = FALSE,label = "Download",icon = icon("save",lib = "glyphicon"),status = "primary")
    })
    
    output$download_menu_table_driver <- renderUI({
      dropdownButton(
        
        tags$div(selectInput("download_data", "Download data: ",   
                             choices = c("CSV","TXT"),selected = "CSV",
                             multiple = FALSE,width = "100%"), style="display:inline-block"),
        downloadBttn(
          outputId = "ok_download_data",
          label = "OK",
          color = "primary",
          style = "pill",icon = icon("ok",lib = "glyphicon"),size = "sm"),
        circle = FALSE,label = "Download",icon = icon("save",lib = "glyphicon"),status = "primary")
    })
  })
  
  width <- reactive({
    ifelse(input$hla == "B", 860,800)
  })
  # height <- reactive({
  #   if (input$need_thre){
  #     1600
  #   }else{
  #     800
  #   }
  # })
  
  output$heatmap <- renderPlotly({
    out_plot()
  })
  
  output$ok_download_image <- downloadHandler(
    filename = function(){
      paste(paste0(gene_name(),"_",aa_change()), tolower(input$download_image), sep=".")
    },
    content = function(file){
      save_image(p1,file = file)
    }
  )
  
  output$ok_download_data <- downloadHandler(
    filename = function(){
      paste(paste0(gene_name(),"_",aa_change()), tolower(input$download_data), sep=".")
    },
    
    content = function(file){
      if(input$download_data == "TXT"){
        write.table(data(),file = file, sep = "\t",col.names = T,row.names = F,quote = F)
      }else{
        write.table(data(),file = file, sep = ",",col.names = T,row.names = F,quote = F)
      }
    }
  )

  output$ok_download_data2 <- downloadHandler(
    filename = function(){
      paste(paste0(gene_name(),"_",aa_change(),"_candidate"), tolower(input$download_data2), sep=".")
    },

    content = function(file){
      if(input$download_data2 == "TXT"){
        write.table(data_candidate(),file = file, sep = "\t",col.names = T,row.names = F,quote = F)
      }else{
        write.table(data_candidate(),file = file, sep = ",",col.names = T,row.names = F,quote = F)
      }
    }
  )
  
  output$tbl2 = renderDT({
    
    data() %>% 
      select(-c("ext_seqs_mt","pep_start","pep_end","core","icore","pos_alter","cdna","transcript","ext_seqs_wt"))
  }, 
  options = list(lengthChange = TRUE,pageLength = 10), selection = 'none',width=400
  )
  
  output$tbl3 = renderDT({

    dt <- data_candidate() %>%
      select(-c("ext_seqs_mt","pep_start","pep_end","ext_seqs_wt",
                "chr","start","end","length","ic50","wt_ic50","transcript","pos_alter")) %>%
      dplyr::rename(ic50_rank=rank,wt_ic50_rank=wt_rank)
    colnames(dt) <- c("Gene", "Ref", "Alt", "HLA allele", "Peptide (MT)",
                      "IC50 Rank (MT)", "Peptide (WT)", "IC50 Rank (WT)",
                      "Stab", "TAP & Cleav", "DeepImmuno", "PRIME", "Seq2Neo-CNN", "Immuno-GNN")
    dt
  },
  options = list(lengthChange = TRUE,pageLength = 10,autoWidth = TRUE), selection = 'none'
  )

  observeEvent(input$switch_db, {
    updateTabsetPanel(session =session ,inputId = "neo_switch", selected = "search_neo")
  })
  
  output$OV <- renderSlickR({
    imgs <- list.files("Fig/ov", pattern=".png", full.names = TRUE)
    slickR(imgs,height = "400px",width = "800px")
  })
  
  ##neodb
  observeEvent(c(input$search_HLA,input$search_gene),{
    output$dt_or_plot <- renderUI({
      shinycssloaders::withSpinner(DT::DTOutput('tbl_db'))
    })
    output$plot_table <- renderUI({tags$br()})
  })
  
  observeEvent(c(input$search_seq),{
    output$dt_or_plot <- renderUI({
      shinycssloaders::withSpinner(plotOutput("ggmsa_plot"))
    })
    output$plot_table <- renderUI({
      DT::DTOutput('tbl_db2')
    })
  })
  
  observeEvent(input$switch_db, {
    updateTabsetPanel(session =session ,inputId = "neo_switch", selected = "search_neo")
  })
  
  output$OV <- renderSlickR({
    imgs <- list.files("Fig/ov", pattern=".png", full.names = TRUE)
    slickR(imgs,height = "400px",width = "800px")
  })
  
  hla_genes <- reactive({
    dt <- neodb_all %>% 
      filter(super_type == input$hla_suptyper) %>%
      filter(hla_type != "H2")
    unique(dt$hla_type)
  })
  
  hla_alleles <- reactive({
    dt <- neodb_all %>% 
      filter(super_type == input$hla_suptyper) %>% 
      filter(hla_type == input$hla_gene)
    unique(dt$HLA)
  })
  
  pmid <- reactive({
    dt <- neodb_all %>% 
      filter(super_type == input$hla_suptyper) %>% 
      filter(hla_type == input$hla_gene) %>% 
      filter(HLA %in% input$hla_allels)
    c("ALL",unique(dt$PMID))
  })
  
  get_seqs <- eventReactive(input$search_seq,{
    dt <- neodb_all %>% 
      select(PMID,`Mut Epitope`,`AA change`,HLA,`WT peptide`,Gene) %>% 
      distinct_all() 
    seq1 <- input$pep_seq
    seq2 <- dt$`Mut Epitope`
    scores <- vector("numeric",length = length(seq2))
    for (i in 1:length(scores)){
      globalAligns <- pairwiseAlignment(seq1, seq2[i], substitutionMatrix = "BLOSUM50", gapOpening = -2,
                                        gapExtension = -8, scoreOnly = FALSE)
      scores[i] <- globalAligns@score
    }
    globalAligns <- pairwiseAlignment(seq1, seq2[which.max(scores)], substitutionMatrix = "BLOSUM50", gapOpening = -2,
                                      gapExtension = -8, scoreOnly = FALSE)
    res <- list(seq1=globalAligns@pattern  %>% as.character(),
                subject= globalAligns@subject %>% as.character(),
                score=scores[which.max(scores)])
    return(res)
  })
  
  select_neo <- reactiveVal()
  observeEvent(input$search_HLA, { 
    dt <- neodb_all %>% 
      filter(super_type == input$hla_suptyper) %>% 
      filter(hla_type == input$hla_gene) %>% 
      filter(HLA %in% input$hla_allels) %>% 
      select(-hla_type,-super_type)
    if (input$pmid != "ALL"){
      dt <- dt %>% 
        filter(PMID %in% input$pmid)
    }
    select_neo(dt)
  }
  )
  
  observeEvent(input$search_gene, { 
    dt <- neodb_all %>% 
      filter(Gene %in% input$gene) %>% 
      select(-hla_type,-super_type)
    select_neo(dt)
  }
  )
  
  observeEvent(input$hla_suptyper,{
    updatePickerInput(session = session,inputId = "hla_gene",
                      choices = hla_genes(),
                      selected = hla_genes()[1])
  })
  
  observeEvent(c(input$hla_suptyper,input$hla_gene),{
    updatePickerInput(session = session,inputId = "hla_allels",
                      choices = hla_alleles(),
                      selected = hla_alleles()[1])
  })
  
  observeEvent(c(input$hla_suptyper,input$hla_gene,input$hla_allels),{
    updatePickerInput(session = session,inputId = "pmid",
                      choices = pmid(),
                      selected = pmid()[1])
  })
  
  observeEvent(c(input$search_HLA,input$search_gene),{
    output$tbl_db <-  renderDT({
      select_neo()
    }, 
    options = list(lengthChange = TRUE,pageLength = 20), selection = 'none',width=400
    )
  })
  
  output$ggmsa_plot <- renderPlot({
    globalAligns <- get_seqs()
    aa <- list(`Input seq`=globalAligns$seq1,`Subject seq`=globalAligns$subject)
    seqinr::write.fasta(aa,names = names(aa),file.out = "~/ShinyApps/Neodb/data/test_msa.fasta")
    test_msa <- readAAStringSet("~/ShinyApps/Neodb/data/test_msa.fasta")
    ggmsa(test_msa, char_width = 0.5, seq_name = T)+
      labs(title=paste0("Score = ",globalAligns$score))
    
  })
  
  observeEvent(input$search_seq,{
    output$tbl_db2 <- renderDT({
      globalAligns <- get_seqs()
      dt <- neodb_all %>% 
        filter(`Mut Epitope` == globalAligns$subject)
      dt
    },
    options = list(lengthChange = TRUE,pageLength = 20), selection = 'none',width=400)
  })

  ##download example
  output$download_maf_exp <- downloadHandler(
    filename = function(){
      "test.maf"
    },
    content = function(file){
      file.copy("./data/test.maf", file)
    }
  )
  
  output$download_vcf_exp <- downloadHandler(
    filename = function(){
      "test.vcf"
    },
    content = function(file){
      file.copy("./data/test.vcf", file)
    }
  )
  
  output$download_txt_exp <- downloadHandler(
    filename = function(){
      "test.txt"
    },
    content = function(file){
      file.copy("./data/test.txt", file)
    }
  )
  
  output$download_pep_exp <- downloadHandler(
    filename = function(){
      "test.pep"
    },
    content = function(file){
      file.copy("./data/random.pep", file)
    }
  )
  output$download_neodb_exp <- downloadHandler(
    filename = function(){
      "neodb_all.csv"
    },
    content = function(file){
      file.copy("./data/neodb_all.csv", file)
    }
  )
}
shinyApp(ui, server)
