library(shiny)
library(ggplot2)
library(purrr)
library(shinythemes)
library(shinyWidgets)
library(MHCbinding)
library(DT)
library(dplyr)
library(shinyFeedback)
library(waiter)
library(slickR)
library(shinymaterial)
library(Biostrings)
library(ggmsa)
library(seqinr)
neodb_all <- readxl::read_xlsx("data/filter_db.xlsx",sheet = 1) %>% select(-`Immunizing peptide`)
colnames(neodb_all)[1:10] <- c("PMID","Sample ID","Mut Epitope","AA change","HLA","WT","Gene","T cell resource","Assay",
                         "Cancer type")
neodb_all <- neodb_all %>% filter(Gene != "NA") %>% filter(Gene != "ERBB2IP<U+00A0") %>% 
	filter(!(HLA %in% c("Class I","Class II","Class I and II")))
neodb_switch <- tabsetPanel(
  id = "neo_switch",
  type = "hidden",
  selected = "overview",
  tabPanel("overview",
           br(),
           br(),
           br(),
           fluidRow(
             column(12,
                    div(slickROutput("OV",height = "400px",width = "800px"), style="text-align: center;"))
           ),
           br(),
           br(),
           br(),
           fluidRow(
             column(5),
             column(4,tags$div(actionBttn(
               inputId = "switch_db",
               label = "Search for validated neoantigens",
               color = "primary",
               style = "fill",size = "lg"
             ),  style="display:inline-block"))
           )
  ),
  tabPanel("search_neo", 
           sidebarLayout(
             sidebarPanel(
               width = 3,
               span(tags$i(h4("Search by HLA")), style="color:#045a8d"),
               pickerInput("hla_type", "HLA Type: ",   
                           choices = c("HLA-I","HLA-II"), 
                           selected = c("HLA-I"),
                           multiple = FALSE),
               pickerInput("hla_gene", "HLA Gene: ",   
                           choices = c("A","B","C"),
                           selected = "A",
                           multiple = FALSE),
               pickerInput("hla_allels", "HLA Allele: ",   
                           choices = c("A01:01"),
                           selected = "A01:01",
                           options = list(
                             `actions-box` = TRUE, 
                             size = 10,
                             `selected-text-format` = "count > 3",
                             pickerOptions(liveSearch = TRUE)
                           ),
                           multiple = TRUE),
               pickerInput("pmid", "PMID: ",   
                           choices = "a",
                           multiple = TRUE,
                           options = list(
                             `actions-box` = TRUE, 
                             size = 10,
                             `selected-text-format` = "count > 3",
                             pickerOptions(liveSearch = TRUE)
                           )),
               span(tags$i(h4("Search by Genes")), style="color:#045a8d"),
               pickerInput("gene", "Gene: ",   
                           choices = unique(neodb_all$Gene),selected = unique(neodb_all$Gene)[1],
                           multiple = TRUE,
                           options = list(
                             `actions-box` = TRUE, 
                             size = 10,
                             `selected-text-format` = "count > 3",
                             pickerOptions(liveSearch = TRUE)
                           )),
               span(tags$i(h4("Search by Peptide Sequence")), style="color:#045a8d"),
               textInputIcon(inputId = "pep_seq", label = "", placeholder = "Search the similar sequence", 
                             icon  = icon("keyboard")),
               tags$br(),
               tags$div(
                 actionGroupButtons(
                  inputIds = c("search_HLA", "search_gene", "search_seq"),
                  labels = list("HLA", "Gene", "Seq"),
                  status = "primary",fullwidth = TRUE
                 )
             )),
             mainPanel(
               fluidRow(
                 column(11,uiOutput("dt_or_plot")),
                 column(1,uiOutput("download_menu_table"))
               ),
               fluidRow(
                 uiOutput("plot_table")
               )
             )
  ))
)

ui <- bootstrapPage(
  navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
             
             HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">Neoantigen Atlas</a>'), id="nav",
             windowTitle = "Neoantigen Atlas",
        
             tabPanel("Neodb",
                      includeCSS("styles.css"),
                      neodb_switch
                      )
             )
)


server <- function(input, output, session) {
  
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
      filter(super_type == input$hla_type) %>%
      filter(hla_type != "H2")
    unique(dt$hla_type)
  })
  
  hla_alleles <- reactive({
    dt <- neodb_all %>% 
      filter(super_type == input$hla_type) %>% 
      filter(hla_type == input$hla_gene)
    unique(dt$HLA)
  })
  
  pmid <- reactive({
    dt <- neodb_all %>% 
      filter(super_type == input$hla_type) %>% 
      filter(hla_type == input$hla_gene) %>% 
      filter(HLA %in% input$hla_allels)
    c("ALL",unique(dt$PMID))
  })
  
  get_seqs <- eventReactive(input$search_seq,{
    dt <- neodb_all %>% 
      select(PMID,`Mut Epitope`,`AA change`,HLA,`WT`,Gene) %>% 
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
      filter(super_type == input$hla_type) %>% 
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
      filter(Gene %in% input$gene)
    select_neo(dt)
    }
  )
  
  observeEvent(input$hla_type,{
    updatePickerInput(session = session,inputId = "hla_gene",
                      choices = hla_genes(),
                      selected = hla_genes()[1])
  })
  
  observeEvent(c(input$hla_type,input$hla_gene),{
    updatePickerInput(session = session,inputId = "hla_allels",
                      choices = hla_alleles(),
                      selected = hla_alleles()[1])
  })
  
  observeEvent(c(input$hla_type,input$hla_gene,input$hla_allels),{
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
    seqinr::write.fasta(aa,names = names(aa),file.out = "~/MHCbindshiny/data/test_msa.fasta")
    test_msa <- readAAStringSet("~/MHCbindshiny/data/test_msa.fasta")
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
}

shinyApp(ui, server)
