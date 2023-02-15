tab_driver <- tabPanel("Driver search",icon = icon("search"),
                       useShinyFeedback(),
                       includeCSS("styles.css"),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      
                                      span(tags$i(h5("Select a gene and then choose one animo acid change to see mutated peptide aised")), style="color:#045a8d"),
                                      pickerInput("gene_name", "Gene name: ",   
                                                  choices = gene, 
                                                  selected = c("PTEN"),
                                                  multiple = FALSE,
						  pickerOptions(
                                                    `actions-box` = TRUE,
                                                    size = 10,
                                                    liveSearch = TRUE
                                                  )),
                                      pickerInput("aa_change", "Amino acid change: ",   
                                                  choices = "p.R14M",
                                                  multiple = FALSE),
                                      
                                      span(tags$i(h5("Or you can input amino acid change directly (Selected amino acid changes are preferred unless user input is selected).")), style="color:#045a8d"),
                                      textInput("aa_change_text",
                                                "Amino acid change:",width="100%"),
                                      
                                      #span(tags$i(h5("There could be alternate transcripts associated with this Amino acid change, please select one.")), style="color:#045a8d"),
                                      #pickerInput("enst", "Alternate transcripts: ",   
                                      #            choices = "ENST00000371953",
                                      #            multiple = FALSE),
                                      
                                      tags$div(actionBttn(
                                        inputId = "search",
                                        label = "Search",
                                        color = "primary",
                                        style = "bordered",icon = icon("search",lib = "glyphicon")
                                      ),  style="display:inline-block"),
                                      
                                      tags$div("                ",  style="display:inline-block"),
                                      
                                      # tags$div(uiOutput("download_menu"),
                                      #          style="display:inline-block"),
                         ),
                         
                         mainPanel(
                           tabsetPanel(type = "tabs",
                                       tabPanel("Heatmap",
                                                
                                                hr(style = "border-top: 0.2px solid #FFFFFF;"),
                                                fluidRow(
                                                  column(1,
                                                         dropdownButton(
                                                           span(tags$i(h5("Set other parameters")), style="color:#045a8d"),
                                                           tags$div(selectInput("hla", "Select HLA type: ",   
                                                                                choices = c("A","B","C"),
                                                                                multiple = FALSE,width = "100%"), style="display:inline-block"),
                                                           tags$div(radioGroupButtons(
                                                             inputId = "ic50",
                                                             label = "Display IC50 or Rank value: ",
                                                             choices = c("IC50", 
                                                                         "Rank"),
                                                             individual = TRUE,
                                                             checkIcon = list(
                                                               yes = icon("check-circle",style = "color: steelblue"),
                                                               no = icon("circle-o",style = "color: steelblue"))
                                                           ), style="display:inline-block"),
                                                           
                                                           #tags$i(h3("Display IC50 or Rank Threshold?")),
                                                           tags$div(prettySwitch(
                                                             inputId = "need_thre",
                                                             label = "Display IC50 or Rank Threshold?", 
                                                             status = "success",
                                                             fill = TRUE
                                                           ), style="display:inline-block"),
                                                           
                                                           uiOutput("threshold"),
                                                           actionBttn(
                                                             inputId = "ok",
                                                             label = "OK",
                                                             color = "primary",
                                                             style = "pill",icon = icon("ok",lib = "glyphicon"),size = "sm"),
                                                           circle = TRUE, status = "danger", icon = icon("cog"), width = "300px",
                                                           tooltip = tooltipOptions(title = "Click to set other parameters!")
                                                         )),
                                                  column(10,
                                                         shinycssloaders::withSpinner(plotlyOutput("heatmap",height="720px",width = "870px"))),
                                                  
                                                  column(1,
                                                         uiOutput("download_menu_image"))
                                                  
                                                )
                                       ),
                                       tabPanel("Table",
                                                hr(style = "border-top: 0.2px solid #FFFFFF;"),
                                                fluidRow(
                                                  column(11,shinycssloaders::withSpinner(DT::DTOutput('tbl2'))),
                                                  column(1,uiOutput("download_menu_table_driver"))
                                                )
                                       ),
				       tabPanel("Candidate Neoantigen",
                                                hr(style = "border-top: 0.2px solid #FFFFFF;"),
                                                fluidRow(
                                                  column(11,shinycssloaders::withSpinner(DT::DTOutput('tbl3'))),
                                                  column(1,uiOutput("download_menu_table_candidate"))
                                                )
                                       )
                           )
                         )
                       )
)
