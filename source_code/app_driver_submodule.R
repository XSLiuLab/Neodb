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

anchar <- readRDS("data/example_anchar.rds")
gene <- unique(anchar$genes)

plot_heatmap <- function(dt,hla_type,need_value,legend_title,need_thre=FALSE){
  dt <- dt %>% 
    select(length,allele,peptide,need_value) %>% 
    filter(grepl(paste0("HLA-",hla_type),allele))
  res <- vector("list",7)
  names(res) <- unique(dt$length)
  for (i in seq_along(res)){
    df <- dt %>% 
      filter(length == names(res)[i]) %>% 
      select(-length) %>% 
      tidyr::pivot_wider(names_from = "allele",values_from = need_value)
    df <- as.data.frame(df)
    rownames(df) <- df$peptide
    df <- df %>% select(-peptide)
    if (need_thre){
      col_fun = colorRamp2(c(0, max(dt[,need_value])), c("white", "blue"))
      h1 <- Heatmap(df,col = col_fun,cluster_rows = F,cluster_columns = F,
                    rect_gp = gpar(col = "grey", lwd = 2),row_names_side = "left",
                    show_heatmap_legend=FALSE,row_names_gp = gpar(fontsize = 8))
    }else{
      col_fun = colorRamp2(c(0, max(dt[,need_value])), c("white", "red"))
      if (i == 1){
        h1 <- Heatmap(df,col = col_fun,cluster_rows = F,cluster_columns = F,
                      rect_gp = gpar(col = "grey", lwd = 2),row_names_side = "left",
                      heatmap_legend_param=list(title = legend_title,
                                                direction = "horizontal",
                                                legend_width = unit(3, "cm"),
                                                title_position = "topcenter"),row_names_gp = gpar(fontsize = 8))
      }else{
        h1 <- Heatmap(df,col = col_fun,cluster_rows = F,cluster_columns = F,
                      rect_gp = gpar(col = "grey", lwd = 2),row_names_side = "left",
                      show_heatmap_legend=FALSE,row_names_gp = gpar(fontsize = 8))
      }
    }
    res[[i]] <- h1
  }
  
  ht_list <- res[[1]] %v%  res[[2]] %v%  res[[3]] %v%  res[[4]] %v%  res[[5]] %v%  res[[6]] %v%  res[[7]] 
  return(ht_list)
}


ui <- bootstrapPage(
  navbarPage(
    useShinyFeedback(),
    tabPanel("Driver search",
             includeCSS("styles.css"),
             sidebarLayout(
               sidebarPanel(width = 5,
                 
                 span(tags$i(h5("Select a gene and then choose one animo acid change to see mutated peptide aised")), style="color:#045a8d"),
                 pickerInput("gene_name", "Gene name: ",   
                             choices = gene, 
                             selected = c("PTEN"),
                             multiple = FALSE),
                 pickerInput("aa_change", "Amino acid change: ",   
                             choices = "p.R14M",
                             multiple = FALSE),
                 
                 span(tags$i(h5("Or you can input amino acid change directly (Selected amino acid changes are preferred unless user input is selected).")), style="color:#045a8d"),
                 textAreaInput("aa_change_text",
                               "Amino acid change:"),
                 
                 span(tags$i(h5("There could be alternate transcripts associated with this Amino acid change, please select one.")), style="color:#045a8d"),
                 pickerInput("enst", "Alternate transcripts: ",   
                             choices = "ENST00000371953",
                             multiple = FALSE),
                 
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
                                        column(2,
                                               dropdownButton(
                                                 span(tags$i(h5("Set other parameters")), style="color:#045a8d"),
                                                 tags$div(selectInput("hla", "Select HLA type: ",   
                                                                      choices = c("A","B","C"),
                                                                      multiple = FALSE,width = "100%"), style="display:inline-block"),
                                                 tags$div(radioGroupButtons(
                                                   inputId = "ic50",
                                                   label = "Display IC50 or Rank value: ",
                                                   choices = c("IC50", 
                                                               "Rank"),selected="IC50",
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
                                        column(9,
                                               shinycssloaders::withSpinner(plotOutput("heatmap",height="800px",width = "600px"))),
                                        
                                        column(1,
                                               uiOutput("download_menu_image"))
                                        
                                      )
                                      ),
                             tabPanel("Table",
                                      hr(style = "border-top: 0.2px solid #FFFFFF;"),
                                      fluidRow(
                                        column(11,shinycssloaders::withSpinner(DT::DTOutput('tbl2'))),
                                        column(1,uiOutput("download_menu_table"))
                                      )
                                    )
                             )
               )
             )
    )
  )
)


server <- function(input, output, session) {
  
  ###一个发生改变都要更新
  # toListen_enst <- reactive({
  #   list(input$gene_name,input$aa_change,input$aa_change_text)
  # })
  
  aa <- reactive({
    dt <- anchar %>% 
      filter(genes == input$gene_name)
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
    updateTextAreaInput(session, "aa_change_text",
                        value = aa()[2])
  })
  
  observeEvent(input$aa_change_text,{
    currentaa <- anchar[anchar$genes == input$gene_name,"aa"]
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
  
  observeEvent(c(input$input_thre),{
    range <- ifelse(input$ic50 == "IC50",Inf,100)
    print(range)
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
  
  data <- reactive({
    currentaa <- anchar[anchar$genes == input$gene_name,"aa"]
    req((aa_change() %in% currentaa$aa))
    dt <- readRDS(paste0("data/test/",gene_name(),"/",aa_change(),".rds"))
    return(dt)
  })

  observeEvent(c(input$gene_name,input$aa_change,input$aa_change_text),{
    
    dt <- data() %>% as.data.frame()
    enst <- unique(dt$contain_trans)
    updatePickerInput(session = session,inputId = "enst",
                      choices = enst,
                      selected = enst[1])
  })
  
  data_dt <- reactive({
    dt <- readRDS(paste0("data/test/",gene_name(),"/",aa_change(),".rds"))
    dt <- dt %>% filter(contain_trans == input$enst)
    return(dt)
  })
  
  out_plot <- eventReactive(c(input$search,input$hla,input$ok),{
    gene_name <- reactive({
      input$gene_name
    })
    aa_change <- reactive({
      ifelse(input$aa_change == "User input",
             input$aa_change_text,input$aa_change)
    })
    
    dt <- readRDS(paste0("data/test/",gene_name(),"/",aa_change(),".rds"))
    dt <- dt %>% filter(contain_trans == input$enst)
    dt <- dt %>% 
      arrange(length)
    dt1 <- dt
    dt$ic50 <- 1 -(log(dt$ic50)/log(50000))
    dt$rank <- 100 - dt$rank
    
    thre <- input$input_thre
    if(is.null(thre)){thre <- 500}
    if(input$ic50 == "IC50"){
      ht_list <- plot_heatmap(dt = dt,hla_type = input$hla,
                              need_value = "ic50",legend_title = "1-log(IC50)/log(50000)")
      if (input$need_thre){
        dt1$ic50 <- ifelse(dt1$ic50 < as.numeric(thre),1,0)
        ht_list1 <- plot_heatmap(dt = dt1,hla_type = input$hla,
                                 need_value = "ic50",need_thre = TRUE)
        p1 <- grid::grid.grabExpr(draw(ht_list,merge_legends =T,heatmap_legend_side = "top"))
        p2 <- grid::grid.grabExpr(draw(ht_list1,merge_legends =T,heatmap_legend_side = "top"))
        plot_grid(p1,p2,ncol = 1)
        
      }else{
        draw(ht_list,merge_legends =T,heatmap_legend_side = "top")
      }
    }else{
      ht_list <- plot_heatmap(dt = dt,hla_type = input$hla,
                              need_value = "rank",legend_title = "1 - %Rank")
      if (input$need_thre){
        dt1$rank <- ifelse(dt1$rank < as.numeric(thre),1,0)
        ht_list1 <- plot_heatmap(dt = dt1,hla_type = input$hla,
                                 need_value = "rank",need_thre = TRUE)
        p1 <- grid::grid.grabExpr(draw(ht_list,merge_legends =T,heatmap_legend_side = "top"))
        p2 <- grid::grid.grabExpr(draw(ht_list1,merge_legends =T,heatmap_legend_side = "top"))
        plot_grid(p1,p2,ncol = 1)
        
      }else{
        draw(ht_list,merge_legends =T,heatmap_legend_side = "top")
      }
    }
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
    
    output$download_menu_table <- renderUI({
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
    ifelse(input$hla == "B", 800,600)
  })
  height <- reactive({
    if (input$need_thre){
      1600
    }else{
      800
    }
  })
  
  observe({
    output$heatmap <- renderPlot({
      out_plot()
    },width = width(),height = height())
  }
  )
  
  output$ok_download_image <- downloadHandler(
    filename = function(){
      paste(paste0(gene_name(),"_",aa_change()), tolower(input$download_image), sep=".")
    },
    content = function(file){
      if(input$download_image == "PNG"){
        png(file,width=width(),height = height()) 
      }else{
        pdf(file,width=width()/100,height = height()/86)
      }
      print(out_plot())
      dev.off()
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
  
  output$tbl2 = renderDT({
    
    data_dt() %>% 
      select(-c("ext_seqs_mt","pep_start","pep_end","core","icore"))
  }, 
  options = list(lengthChange = TRUE,pageLength = 10), selection = 'none',width=400
  )
  
}


shinyApp(ui, server)