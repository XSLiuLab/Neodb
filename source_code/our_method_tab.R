immuno_gnn <- function(input_type="single", input_file, pep, hla,
                       immuno_gnn_env, immuno_gnn_path, temp_dir){
  input_type <- match.arg(arg = input_type, choices = c("single","multiple"))
  pseudu <- MHCbinding::pseudu
  if (input_type == "multiple"){
    pep_file <- read.csv(input_file)
  }else{
    pep_file <- data.frame(pep=pep,HLA_Allele=hla)
  }
  pep_file <- pep_file %>%
    rowwise() %>%
    mutate(hla_seq=unique(pseudu$V2[pseudu$V1==HLA_Allele])) %>%
    ungroup() %>%
    mutate(type=0)
  write.csv(pep_file,file = paste0(temp_dir,"/pep_file.csv"),quote = FALSE,row.names = FALSE)
  system(paste0(". ~/.bashrc;conda run -n ",immuno_gnn_env," python ",
                normalizePath(immuno_gnn_path),"/Immuno_gnn.py ",
                "-a ",normalizePath(immuno_gnn_path),"/aaindex1_pca.csv ",
                "-m ",normalizePath(immuno_gnn_path),"/last_model.pt ",
                " -i ",paste0(temp_dir,"/pep_file.csv"),
                " -o ",temp_dir))
  res <- read.csv(paste0(temp_dir,"/pred_res.csv"))
  res <- res %>% select(-X,-label,-pred,-hla_seq,-type)
  colnames(res)[3] <- "immunogenicity"
  return(res)
}


immuno_gnn_tab_ui <- function(id){
  ns <- NS(id)
  alleles <- MHCbinding::available_alleles("Immuno","I","Immuno-GNN")
  input_tabs <- tabsetPanel(
    id = ns("upload"),
    type = "hidden",
    tabPanel("Peptides (paste)",
             pickerInput(ns("HLA_allele"), "HLA allele:",   
                         choices = alleles, 
                         options = list(
                           `actions-box` = TRUE, 
                           size = 10,
                           `selected-text-format` = "count > 3"
                         ),
                         multiple = FALSE),
             textInput(ns("pep"),"Peptides: ",value = "SLYNTVATLY")
    ),
    tabPanel("Peptide (file)", 
             fileInput(ns("upload_file_pep"), 
                       "Please upload file (csv file with pep and HLA_Allele columns): ", 
                       buttonLabel = "Upload file ...")
    )
  )
  tabPanel(id,
           sidebarLayout(
             sidebarPanel(width = 3,
                          
                          pickerInput(ns("input_type"), "Input type: ",
                                      choices = c("Peptides (paste)", "Peptide (file)"), 
                                      selected = c("Peptides (paste)"),
                                      multiple = FALSE),
                          
                          input_tabs,
                          
                          tags$div(actionBttn(
                            inputId = ns("bttn1"),
                            label = "Go!",
                            color = "primary",
                            style = "bordered"
                          ),  style="display:inline-block"),
                          
                          tags$div("                ",  style="display:inline-block"),
                          
                          tags$div(downloadBttn(
                            outputId = ns("downloadData"),
                            style = "bordered",
                            color = "primary"
                          ),  style="display:inline-block"),
                          
                          span(tags$i(h5("When select 'peptide(paste)' input, please input the peptide sequence in the input box (separated by comma, if there are more than two sequences), otherwise please upload corresponding file.")), style="color:#045a8d"),
             ),
             
             mainPanel(
               shinycssloaders::withSpinner(DT::DTOutput(ns('tbl')))
             )
           ) 
           
  )
}

immuno_gnn_tab <- immuno_gnn_tab_ui("Immuno-GNN")

immuno_gnn_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    data_upload <- reactive({
      
      if (input$input_type == "Peptide (file)"){
        req(input$upload_file_pep)
        ext <- tools::file_ext(input$upload_file_pep$name)
        switch(ext,
               csv = input$upload_file_pep$datapath,
               validate("Invalid file; Please upload a .csv file")
        )
      }
    })
    
    observeEvent(input$input_type, {
      updateTabsetPanel(inputId = "upload", selected = input$input_type)
    })
    
    ##输出
    out_dt <- eventReactive(input$bttn1, {
      
      ##参数
      if (input$input_type == "Peptides (paste)"){
        input_type <- "single"
      }else{
        input_type <- "multiple"
      }
      
      input_file <- data_upload()
      pep <- input$pep %>% strsplit(.,split = ",") %>% unlist()
      hla <- input$HLA_allele
      immuno_gnn_env <- "dl_gpu"
      immuno_gnn_path <- "/home/ubuntu/software/MHCbinding/data-raw/"
      temp_dir <- tempdir()
      
      if (input_type == "single"){
        res <- immuno_gnn(pep = pep,
                          hla = hla, immuno_gnn_env=immuno_gnn_env,
                          immuno_gnn_path=immuno_gnn_path,
                          temp_dir=temp_dir)
      }else{
        res <- immuno_gnn(input_type = "multiple",input_file=input_file,
                          immuno_gnn_env=immuno_gnn_env,
                          immuno_gnn_path=immuno_gnn_path,
                          temp_dir=temp_dir)
      }
      return(res)
    }
    )
    
    
    output$tbl <- renderDT({
      
      if (is.null(input$pep)){
        validate("Please input peptide sequence !")
      }
      if (input$input_type == "Peptides (paste)"){
        out_dt()
      }
    },
    options = list(lengthChange = TRUE,pageLength = 15), selection = 'none',rownames= FALSE
    )
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('query_data-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(out_dt(), con)
      }
    )
  })
}
