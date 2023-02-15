library(shiny)
library(ggplot2)
library(purrr)
library(shinythemes)
library(shinyWidgets)
library(MHCbinding)
library(DT)
library(dplyr)
library(waiter)

waiter_set_theme(html = spin_rotating_plane(), color = "#333e48")
input_tabs <- tabsetPanel(
  id = "upload",
  type = "hidden",
  tabPanel("Peptides (paste)",
           textAreaInput("pep","Peptides: ")
  ),
  tabPanel("Peptide (file)", 
           fileInput("upload_file_pep", "Please upload Peptide file (with extension .pep): ", buttonLabel = "Upload Peptide ...")
  ),
  tabPanel("VCF",
           fileInput("upload_file_vcf", "Please upload VCF file (with extension .vcf): ", buttonLabel = "Upload VCF ..."),
           pickerInput("genome_version_vcf", "Genome version for VCF: ",
                       choices = c("hg19","hg38"), 
                       selected = c("hg38"),
                       multiple = FALSE),
           textAreaInput("need_samples_vcf","Samples need to be processed in VCF ('ALL' for all samples):")
  ),
  tabPanel("MAF",
           fileInput("upload_file_maf", "Please upload MAF file (with extension .maf): ", buttonLabel = "Upload MAF ..."),
           textAreaInput("need_samples_maf","Samples need to be processed in MAF ('ALL' for all samples):")
  ),
  tabPanel("TXT",
           fileInput("upload_file_txt", "Please upload TXT file (with extension .txt): ", buttonLabel = "Upload TXT ..."),
           pickerInput("genome_version_txt", "Genome version for TXT: ",
                       choices = c("hg19","hg38"), 
                       selected = c("hg38"),
                       multiple = FALSE)
  )
)

ui <- bootstrapPage(
  navbarPage(
    waiter::use_waiter(),
    tabPanel("Driver search",
             sidebarLayout(
               sidebarPanel(
                 pickerInput("hla_type", "MHC class: ",   
                             choices = c("MHC-I", "MHC-II"), 
                             selected = c("MHC-I"),
                             multiple = FALSE),
                 
                 pickerInput("input_type", "Input type: ",
                             choices = c("Peptides (paste)", "Peptide (file)", "VCF", "MAF", "TXT"), 
                             selected = c("Peptides (paste)"),
                             multiple = FALSE),
                 
                 pickerInput("method", "Method: ",   
                             choices = c("a"), 
                             #selected = c("HLA-I"),
                             multiple = FALSE),
                 
                 pickerInput("HLA_allele", "HLA allele:",   
                             choices = c("a"), 
                             options = list(
                               `actions-box` = TRUE, 
                               size = 10,
                               `selected-text-format` = "count > 3"
                             ),
                             multiple = TRUE),
                 
                 pickerInput("length", "Predicted core length: ",   
                             choices = c("11"), 
                             selected = c("11"),
                             multiple = TRUE),
                 
                input_tabs,
                
                tags$div(actionBttn(
                  inputId = "bttn1",
                  label = "Go!",
                  color = "primary",
                  style = "bordered"
                ),  style="display:inline-block"),
                
                tags$div("                ",  style="display:inline-block"),
                
                tags$div(downloadBttn(
                  outputId = "downloadData",
                  style = "bordered",
                  color = "primary"
                ),  style="display:inline-block"),
                
                span(tags$i(h5("When select MHC-II, the format of allele like DPA1*01/DPB1*04:01 refers to allele(s) with α and β chains.")), style="color:#045a8d"),
                span(tags$i(h5("When select 'peptide(paste)' input, please input the peptide sequence in the input box (separated by comma, if there are more than two sequences), otherwise please upload corresponding file.")), style="color:#045a8d"),
                span(tags$i(h5("When upload file, you can specify the samples to be processed or use 'ALL' for processing all samples")),style="color:#045a8d")
               ),
               
               mainPanel(
                 DT::DTOutput('tbl')
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  #waiter <- waiter::Waiter$new("tbl")
  
  toListen <- reactive({
    list(input$hla_type,input$HLA_allele)
  })
  
  pre_method <- reactive({
    dt <- MHCbinding::available_methods(get_methods = "api",pre_type = input$hla_type)
    dt[which(dt=="recommended")] <- "IEDB_recommended"
    return(dt)
  })
  
  alleles <- reactive({
    dt <- MHCbinding::available_alleles(pre_type = input$hla_type,pre_method = input$method)
    return(dt$alleles)
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
  
  out_dt <- eventReactive(input$bttn1, {
    waiter <- waiter::Waiter$new()
    waiter$show()
    on.exit(waiter$hide())
    
    if (input$input_type == "Peptides (paste)"){
      pep <- input$pep
      pep <- strsplit(pep,split = ",")[[1]]
      
      res <- MHCbinding:::general_mhcbinding(mhc_type = input$hla_type,peptide = pep,
                                             allele = input$HLA_allele,length = input$length,
                                             pre_method = input$method,get_method = "client")
      return(res)
    }else if (input$input_type == "Peptide (file)"){
      
      res <- MHCbinding:: batchpep_binding(pep_file=data_upload(),
                                           mhc_type=input$hla_type,pep_length= input$length,
                                           allele=input$HLA_allele,pre_method=input$method)
      
      
    }else if (input$input_type == "VCF"){
      if ( input$need_samples_vcf == "ALL"){
        need_allsamples <- TRUE
      }else{
        need_allsamples <- FALSE
      }
      res <- MHCbinding::vcf2binding(annovar_path = "~/software/annovar/",vcf_path = data_upload(),
                                     genome_version = input$genome_version_vcf,need_allsamples = need_allsamples,
                                     need_samples = input$need_samples_vcf,
                                     mhc_type = input$hla_type,pep_length = input$length,
                                     allele = input$HLA_allele,pre_method = input$method)
      return(res)
    }else if (input$input_type == "MAF"){
      if ( input$need_samples_maf == "ALL"){
        need_allsamples <- TRUE
      }else{
        need_allsamples <- FALSE
      }
      res <- MHCbinding::maf2binding(annovar_path = "~/software/annovar/",maf_path = data_upload(),
                                     need_allsamples = need_allsamples,need_samples = input$need_samples_maf,
                                     mhc_type = input$hla_type,pep_length = input$length,
                                     allele = input$HLA_allele,pre_method = input$method)
    }else{
      res <- MHCbinding::txt2binding(annovar_path = "~/software/annovar/",txt_path = data_upload(),
                                     genome_version = input$genome_version_txt,
                                     mhc_type = input$hla_type,pep_length = input$length,
                                     allele = input$HLA_allele,pre_method = input$method)
      
      
    }
  }
)
  
  observeEvent(input$hla_type,{
    updatePickerInput(session = session,inputId = "method",
                      choices = pre_method(),
                      selected = pre_method()[1])
  })
  
  observeEvent(input$method,{
    updatePickerInput(session = session,inputId = "HLA_allele",
                      choices = alleles(),
                      selected = alleles()[1])
  })
  
  observeEvent(input$input_type, {
    updateTabsetPanel(inputId = "upload", selected = input$input_type)
  })
  
  observeEvent(toListen(), {
    ll <- reactive({
      if (input$hla_type == "MHC-I"){
        c("8","9", "10", "11", "12", "13", "14", "15")
      }else{
        c("11","12","13","14","15","16","17","18","19",
          "20","21","22","23","24","25","26","27","28","29","30","asis")
      }
    })
    
    # max_selection <- reactive({
    #   if (input$hla_type == "MHC-I"){
    #     length(input$HLA_allele)
    #   }else{
    #     21
    #   }
    # })
    # 
    # lable_mhci <- reactive({
    #   if (input$hla_type == "MHC-I"){
    #     "Selectd MHC-I, need paired core lengths: "
    #   }else{
    #     "Input core peptide length: "
    #   }
    # })
    
    updatePickerInput(session = session,inputId = "length", label = "Input core peptide length: ",
                      choices = ll(),
                      selected = ll()[1],
                      # options =  list(
                      #   "max-options" = max_selection(),
                      #   "max-options-text" = "If select MHC-I, the length must paired with number of alleles !"
                      # )
                      )
  })
  
  output$tbl = renderDT({
    
    if (is.null(input$pep)){
      validate("Please input peptide sequence !")
    }
    
    out_dt() 
  }, 
    options = list(lengthChange = TRUE,pageLength = 15), selection = 'none'
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('query_data-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(out_dt(), con)
    }
  )
  
}
shinyApp(ui, server)
