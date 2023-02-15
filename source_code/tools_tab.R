options(digits=3)
##################################################.             
##############NavBar Menu----
###############################################.
tabpanel_ui <- function(id, only_hlai=FALSE){
  ns <- NS(id)
  input_tabs <- tabsetPanel(
    id = ns("upload"),
    type = "hidden",
    tabPanel("Peptides (paste)",
             textInput(ns("pep"),"Peptides: ",value = "SLYNTVATLYY")
    ),
    tabPanel("Peptide (file)", 
             fileInput(ns("upload_file_pep"), "Please upload Peptide file (with extension .pep): ", buttonLabel = "Upload Peptide ...")
    ),
    tabPanel("VCF",
             fileInput(ns("upload_file_vcf"), "Please upload VCF file (with extension .vcf): ", buttonLabel = "Upload VCF ..."),
             pickerInput(ns("genome_version_vcf"), "Genome version for VCF: ",
                         choices = c("hg19","hg38"), 
                         selected = c("hg38"),
                         multiple = FALSE),
             textInput(ns("need_samples_vcf"),"Samples need to be processed in VCF ('ALL' for all samples):",value = "ALL")
    ),
    tabPanel("MAF",
             fileInput(ns("upload_file_maf"), "Please upload MAF file (with extension .maf): ", buttonLabel = "Upload MAF ..."),
             textInput(ns("need_samples_maf"),"Samples need to be processed in MAF ('ALL' for all samples):",value = "ALL")
    ),
    tabPanel("TXT",
             fileInput(ns("upload_file_txt"), "Please upload TXT file (with extension .txt): ", buttonLabel = "Upload TXT ..."),
             pickerInput(ns("genome_version_txt"), "Genome version for TXT: ",
                         choices = c("hg19","hg38"), 
                         selected = c("hg38"),
                         multiple = FALSE)
    )
  )
  
  if (only_hlai){
    hla_choice <- c("I")
  }else{
    hla_choice <- c("I","II")
  }
  
  tabPanel(id,
           sidebarLayout(
             sidebarPanel(width = 3,
                          pickerInput(ns("hla_type"), "HLA class: ",   
                                      choices = hla_choice, 
                                      selected = c("I"),
                                      multiple = FALSE),
                          
                          pickerInput(ns("input_type"), "Input type: ",
                                      choices = c("Peptides (paste)", "Peptide (file)", "VCF", "MAF", "TXT"), 
                                      selected = c("Peptides (paste)"),
                                      multiple = FALSE),
                          
                          pickerInput(ns("method"), "Method: ",   
                                      choices = NULL,
                                      multiple = FALSE),
                          
                          pickerInput(ns("HLA_allele"), "HLA allele:",   
                                      choices = NULL, 
                                      options = list(
                                        `actions-box` = TRUE, 
                                        size = 10,
                                        `selected-text-format` = "count > 3"
                                      ),
                                      multiple = TRUE),
                          
                          pickerInput(ns("length"), "Predicted core length: ",   
                                      choices = NULL, 
                                      selected = NULL,
                                      multiple = TRUE),
                          
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
                          
                          span(tags$i(h5("When select MHC-II, the format of allele like DPA1*01/DPB1*04:01 refers to allele(s) with α and β chains.")), style="color:#045a8d"),
                          span(tags$i(h5("When select 'peptide(paste)' input, please input the peptide sequence in the input box (separated by comma, if there are more than two sequences), otherwise please upload corresponding file.")), style="color:#045a8d"),
                          span(tags$i(h5("When upload file, you can specify the samples to be processed or use 'ALL' for processing all samples")),style="color:#045a8d")
             ),
             
             mainPanel(
               shinycssloaders::withSpinner(DT::DTOutput(ns('tbl')))
             )
           ) 
           
  )
}

tab_tools <- navbarMenu("Tools",icon = icon("toolbox"),
                        tabpanel_ui("Binding", only_hlai=FALSE),
                        tabpanel_ui("Processing", only_hlai=TRUE),
                        tabpanel_ui("Immuno", only_hlai=TRUE))

tools_server <- function(id,module_id) {
  moduleServer(id, function(input, output, session) {
    pre_method <- reactive({
      dt <- MHCbinding::available_methods(module_id)
      if (module_id == "Binding"){
        if (input$hla_type == "I"){
          dt <- dt[[1]]
        }else{
          dt <- dt[[2]]
        }
      }
      if (module_id == "Immuno"){
        dt <- dt[1:4]
      }

      return(dt)
    })
    
    alleles <- reactive({
      dt <- MHCbinding::available_alleles(pre_type = module_id, 
                                          HLA_type = input$hla_type,
                                          pre_method = input$method)
      return(dt)
    })
    
    data_upload <- reactive({
      
      if (input$input_type == "Peptide (file)"){
        req(input$upload_file_pep)
        ext <- tools::file_ext(input$upload_file_pep$name)
        switch(ext,
               pep = input$upload_file_pep$datapath,
               validate("Invalid file; Please upload a .pep file")
        )
      }else if (input$input_type == "VCF"){
        req(input$upload_file_vcf)
        ext <- tools::file_ext(input$upload_file_vcf$name)
        switch(ext,
               vcf = input$upload_file_vcf$datapath,
               validate("Invalid file; Please upload a .vcf file")
        )
      }else if (input$input_type == 'MAF'){
        req(input$upload_file_maf)
        ext <- tools::file_ext(input$upload_file_maf$name)
        switch(ext,
               maf = input$upload_file_maf$datapath,
               validate("Invalid file; Please upload a .maf file")
        )
      }else{
        req(input$upload_file_txt)
        ext <- tools::file_ext(input$upload_file_txt$name)
        switch(ext,
               txt = input$upload_file_txt$datapath,
               validate("Invalid file; Please upload a .txt file")
        )
      }
    })
    
    observeEvent(input$hla_type,{
      updatePickerInput(session = session,inputId = "method",
                        choices = pre_method(),
                        selected = pre_method()[1])
    })
    
    observeEvent(input$method,{
      if (input$method == "Netchop"){
        file_choice <- c("Peptides (paste)")
        allele_choice <- c("Netchop support any alleles")
      }else{
        file_choice <- c("Peptides (paste)", "Peptide (file)", "VCF", "MAF", "TXT")
        allele_choice <- alleles()
      }
      updatePickerInput(session = session,inputId = "input_type",
                        choices = file_choice,
                        selected = file_choice[1])
      updatePickerInput(session = session,inputId = "HLA_allele",
                        choices = allele_choice,
                        selected = allele_choice[1])
    })
    
    observeEvent(input$input_type, {
      updateTabsetPanel(inputId = "upload", selected = input$input_type)
    })
    
    observeEvent(c(input$hla_type,input$HLA_allele,input$method), {
      ll <- reactive({
        if (input$hla_type == "I"){
          lens <- available_len(pre_method = input$method,pre_allele = input$HLA_allele)
          req(input$method)
          if (input$method == "Netchop"){
            lens <- "Netchop support any length"
          }
        }else{
          lens <- c("11","12","13","14","15","16","17","18","19",
                    "20","21","22","23","24","25","26","27","28","29","30","asis")
        }
        return(lens)
      })
      
      updatePickerInput(session = session,inputId = "length", label = "Input core peptide length: ",
                        choices = ll(),
                        selected = ll()[1]
      )
    })
    
    ##输出
    out_dt <- eventReactive(input$bttn1, {
      
      ##参数
      annovar_path <- "~/software/annovar/"
      client_path <- ifelse(input$hla_type == "I" ,"~/software/mhc_i/src/","~/software/mhc_ii/")
      tmp_dir <- tempdir()
      num_thread <- 1
      mhcflurry_env <- "mhcflurry-env"
      mhcnuggets_env <- "mhcnuggets"
      netchop_path <- "~/software/netchop/"
      Immuno_IEDB_path <- "~/software/immunogenicity/"
      Immuno_Deepimmuno_path <- "~/software/DeepImmuno/"
      Deepimmuno_env <- "DeepImmuno"
      MixMHCpred_path <- "~/software/MixMHCpred/MixMHCpred"
      PRIME_path <- "~/software/PRIME-2.0/"
      seq2neo_env <- "DeepImmuno"
      seq2neo_path <- "~/software/Seq2Neo/seq2neo/function/immuno_Prediction/"
      
      req(input$length)
      if (input$method == "Netchop"){
        res <- netchop_processing(pep = input$pep,temp_dir = tempdir(), netchop_path = netchop_path)
        res <- res %>% select(-sequence_id)
        return(res)
      }
      
      if (input$input_type == "Peptides (paste)"){
        pep <- input$pep
        pep <- strsplit(pep,split = ",")[[1]]
        
        res <- MHCbinding:::general_mhcbinding(hla_type = input$hla_type, length = input$length,
                                               allele = input$HLA_allele,pre_method = input$method,
                                               method_type= module_id,
                                               peptide = pep,
                                               tmp_dir=tmp_dir, 
                                               mhcflurry_env = mhcflurry_env,
                                               mhcnuggets_env = mhcnuggets_env,
                                               netchop_path = netchop_path,
                                               Immuno_IEDB_path=Immuno_IEDB_path,
                                               Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                               Deepimmuno_env=Deepimmuno_env,
                                               MixMHCpred_path=MixMHCpred_path,
                                               PRIME_path=PRIME_path,
                                               client_path = client_path,
                                               seq2neo_env = seq2neo_env,
                                               seq2neo_path = seq2neo_path)
        
      }else if (input$input_type == "Peptide (file)"){
        
        res <- MHCbinding::batchpep_binding(pep_file=data_upload(),
                                            hla_type = input$hla_type, pep_length = input$length,
                                            allele = input$HLA_allele,pre_method = input$method,
                                            method_type= module_id,
                                            tmp_dir=tmp_dir, 
                                            mhcflurry_env = mhcflurry_env,
                                            mhcnuggets_env = mhcnuggets_env,
                                            netchop_path = netchop_path,
                                            Immuno_IEDB_path=Immuno_IEDB_path,
                                            Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                            Deepimmuno_env=Deepimmuno_env,
                                            MixMHCpred_path=MixMHCpred_path,
                                            PRIME_path=PRIME_path,
                                            client_path = client_path,
                                            seq2neo_env = seq2neo_env,
                                            seq2neo_path = seq2neo_path)
        
        
      }else if (input$input_type == "VCF"){
        if ( input$need_samples_vcf == "ALL"){
          need_allsamples <- TRUE
        }else{
          need_allsamples <- FALSE
        }
        res <- MHCbinding::vcf2binding(vcf_path = data_upload(),genome_version=input$genome_version_vcf,
                                       need_allsamples=need_allsamples,need_samples=input$need_samples_vcf,
                                       hla_type = input$hla_type, pep_length = input$length,
                                       allele = input$HLA_allele,pre_method = input$method,
                                       method_type= module_id,
                                       tmp_dir=tmp_dir, 
                                       mhcflurry_env = mhcflurry_env,
                                       mhcnuggets_env = mhcnuggets_env,
                                       netchop_path = netchop_path,
                                       Immuno_IEDB_path=Immuno_IEDB_path,
                                       Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                       Deepimmuno_env=Deepimmuno_env,
                                       MixMHCpred_path=MixMHCpred_path,
                                       PRIME_path=PRIME_path,
                                       client_path = client_path,
                                       seq2neo_env = seq2neo_env,
                                       seq2neo_path = seq2neo_path,
				       annovar_path = annovar_path,num_thread = num_thread)
      }else if (input$input_type == "MAF"){
        if ( input$need_samples_maf == "ALL"){
          need_allsamples <- TRUE
        }else{
          need_allsamples <- FALSE
        }
        res <- MHCbinding::maf2binding(maf_path = data_upload(),
                                       need_allsamples=need_allsamples,need_samples=input$need_samples_maf,
                                       hla_type = input$hla_type, pep_length = input$length,
                                       allele = input$HLA_allele,pre_method = input$method,
                                       method_type= module_id,
                                       tmp_dir=tmp_dir, 
                                       mhcflurry_env = mhcflurry_env,
                                       mhcnuggets_env = mhcnuggets_env,
                                       netchop_path = netchop_path,
                                       Immuno_IEDB_path=Immuno_IEDB_path,
                                       Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                       Deepimmuno_env=Deepimmuno_env,
                                       MixMHCpred_path=MixMHCpred_path,
                                       PRIME_path=PRIME_path,
                                       client_path = client_path,
                                       seq2neo_env = seq2neo_env,
                                       seq2neo_path = seq2neo_path,
				       annovar_path = annovar_path,num_thread = num_thread)
      }else{
        res <- MHCbinding::txt2binding(txt_path = data_upload(),
                                       genome_version = input$genome_version_txt,
                                       hla_type = input$hla_type, pep_length = input$length,
                                       allele = input$HLA_allele,pre_method = input$method,
                                       method_type= module_id,
                                       tmp_dir=tmp_dir, 
                                       mhcflurry_env = mhcflurry_env,
                                       mhcnuggets_env = mhcnuggets_env,
                                       netchop_path = netchop_path,
                                       Immuno_IEDB_path=Immuno_IEDB_path,
                                       Immuno_Deepimmuno_path=Immuno_Deepimmuno_path,
                                       Deepimmuno_env=Deepimmuno_env,
                                       MixMHCpred_path=MixMHCpred_path,
                                       PRIME_path=PRIME_path,
                                       client_path = client_path,
                                       seq2neo_env = seq2neo_env,
                                       seq2neo_path = seq2neo_path,
				       annovar_path = annovar_path,num_thread = num_thread)
        
        
      }
      return(res)
    }
    )
    
    
    output$tbl <- renderDT({

      if (is.null(input$pep)){
        validate("Please input peptide sequence !")
      }
      if (grepl("Peptide",input$input_type)){
        out_dt()
      }else{
	req(out_dt())
	req("ext_seqs_mt" %in% colnames(out_dt()))
        out_dt() %>%
          select(-ext_seqs_mt,-ext_seqs_wt,-cdna,-transcript,-pep_start,-pep_end)
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

# server <- function(input, output, session) {
#   tools_server("Binding",module_id = "Binding")
#   tools_server("Processing",module_id = "Processing")
#   tools_server("Immuno",module_id = "Immuno")
# }

# shinyApp(ui, server)
