
neodb_switch <- tabsetPanel(
  id = "neo_switch",
  type = "hidden",
  selected = "overview",
  tabPanel("overview",icon = icon("database"),
           br(),
           br(),
           br(),
           div(
             fluidRow(
               column(12,
                      div(slickROutput("OV",height = "400px",width = "800px"), style="text-align: center;"))
             ),
             br(),
             br(),
             br(),
             fluidRow(
               align = "center",
               br(),
               column(12,
                      actionBttn(
                        inputId = "switch_db",
                        label = "Search for validated neoantigens",
                        color = "primary",
                        style = "gradient",
                        block = FALSE,
                        size = "lg"
                      ))
             )
           )
  ),
  tabPanel("search_neo",icon = icon("database"),
           sidebarLayout(
             sidebarPanel(
               width = 3,
               span(tags$i(h4("Search by HLA")), style="color:#045a8d"),
               pickerInput("hla_suptyper", "HLA Type: ",   
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
                           options = pickerOptions(
                             `actions-box` = TRUE, 
                             size = 10,
                             `selected-text-format` = "count > 3",
                             liveSearch = TRUE
                           ),
                           multiple = TRUE),
               pickerInput("pmid", "PMID: ",   
                           choices = "a",
                           multiple = TRUE,
                           options = pickerOptions(
                             `actions-box` = TRUE, 
                             size = 10,
                             `selected-text-format` = "count > 3",
                             liveSearch = TRUE
                           )),
               span(tags$i(h4("Search by Genes")), style="color:#045a8d"),
               pickerInput("gene", "Gene: ",   
                           choices = unique(neodb_all$Gene),
			   #selected = unique(neodb_all$Gene)[1],
                           multiple = TRUE,
                           options = pickerOptions(
                             `actions-box` = TRUE, 
                             size = 10,
                             `selected-text-format` = "count > 3",
                             liveSearch = TRUE
                           )),
               span(tags$i(h4("Search by Peptide Sequence")), style="color:#045a8d"),
               textInputIcon(inputId = "pep_seq", label = "Search the similar sequence", 
                             value = "TRAATGRMV", 
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

tab_neodb <- tabPanel("Neodb",
                      neodb_switch
)
