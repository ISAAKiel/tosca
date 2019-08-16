library(shiny)
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyr)
library(ggrepel)
library(seriation)

data_input_ui <- function() {
  fluidPage(
    shinyjs::useShinyjs(),
    
    # App title ----
    titlePanel("Uploading Files"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Select a file ----
        fileInput("file1", "Choose CSV File",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        # Horizontal line ----
        tags$hr(),
        fluidRow(
          column(6,  
        # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE)),
        
          column(6,
        # Input: Checkbox if file has header ----
        checkboxInput("rownames", "Rownames", TRUE))),
        # Input: Select separator ----
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ",", inline = T),
        # Input: Select quotes ----
        radioButtons("quote", "Quote",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = '"', inline=T),
        
        # Horizontal line ----
        tags$hr(),
        
        # Input: Select number of rows to display ----
        radioButtons("disp", "Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head", inline=T),
        actionButton(inputId = "go", label = "Go")
        
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        
        # Output: Data file ----
        div(style = 'overflow-x: auto', tableOutput("contents"))
      )
      
    )
  )
}

seriation_ui <- function(x) {
  sidebarLayout(
    sidebarPanel(
      selectInput("seriation_method",
                  "Seriation Method",
                  choices = list("Reciprocal Averaging" = 1,
                                 "Correspondence Analysis" = 2,
                                 "Original Data" = 99),
                  selected = 1),
      
      # Input: Choose dataset ----
      selectInput("download_format_ser",
                  "Download as:",
                  choices = c("svg",
                              "png",
                              "pdf",
                              "csv")),
      
      # Button
      downloadButton("downloadSerPlot",
                     "Download")
    ),
    mainPanel(
      # Output: Data file ----
      plotOutput("plot_ser")
    )
  )
}

analysis_ui <- function(x) {
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(6,
      numericInput("x_axis_dimension",
                   "Dimension x-axis",
                   value = 1,
                   min=1)),
        column(6,
      numericInput("y_axis_dimension",
                   "Dimension y-axis",
                   value = 2,
                   min=1))),
      fluidRow(
      column(4,
      checkboxGroupInput("display_what",
                         "What to display",
                         choices = list("Sites"="sites",
                                        "Types"="species"),
                         selected = c("sites",
                                      "species"))),
      column(4,
      checkboxGroupInput("display_how",
                         "Plot Elements",
                         choices = list("Points"="points",
                                        "Text"="text"),
                         selected = c("points",
                                      "text"))),
      column(4,
      checkboxInput("spead_labels",
                    "Spread labels"))),
      sliderInput("xlim", label = "X-Axis", min = 0, 
                  max = 1, value = c(0, 1)),
      sliderInput("ylim", label = "Y-Axis", min = 0, 
                  max = 1, value = c(0, 1)),
      # Input: Choose dataset ----
      selectInput("download_format",
                  "Download as:",
                  choices = c("svg",
                              "png",
                              "pdf")),
      
      # Button
      downloadButton("downloadPlot",
                     "Download")
    ),
    mainPanel(
      
      # Output: Data file ----
      plotOutput("plot_ca")
    )
  )
}

download_ui <- function() {
  mainPanel(
    
    downloadButton("downloadResult",
                   "Download Coordinates of Points")
    
  )
}

ui <- shinyUI(navbarPage(id="navbar",
                         title = "Tool for Online
                         Seriation and Correspondence Analysis (TOSCA)",
                         tabPanel(id="data_input",
                                  title = "Data Input",
                                  data_input_ui()),
                         tabPanel("CA Analysis",
                                  id="analysis",
                                  analysis_ui()),
                         tabPanel("Seriation",
                                  id="Seriation",
                                  seriation_ui()),
                         tabPanel("Export",
                                  download_ui())
))

# Define server logic to read selected file ----
server <- function(input, output, session) {
  observeEvent(input$go, {
    updateTabsetPanel(session = session,
                      inputId = "navbar",
                      selected = "CA Analysis")
  })
  observe({
    updateCheckboxInput(session,
                        "spead_labels",
                        value = F)
    shinyjs::toggleState("spead_labels",
                         'text' %in% input$display_how)
  })
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   row.names = if(input$rownames==TRUE) 1 else NULL,
                   sep = input$sep,
                   quote = input$quote)
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  },
  rownames=T)
  
  ser <- eventReactive(input$seriation_method, {
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   row.names = if(input$rownames==TRUE) 1 else NULL,
                   sep = input$sep,
                   quote = input$quote)
    
    if(input$seriation_method==1)
    {
      untidy <- TRUE
      max_it <- 1000
      counter <- 0
      rva <- as.matrix(df)
      while(untidy & counter < max_it)
      {
        counter <- counter + 1
        rva_old <- rva
        row_ind_mean <- rowSums(t(t(rva) * 1:ncol(rva)))/rowSums(rva)
        col_ind_mean <- colSums(rva * 1:nrow(rva))/colSums(rva)
        rva <- rva[order(row_ind_mean),order(col_ind_mean)]
        if(all(rva==rva_old)) tidy=TRUE
      }
    } else if(input$seriation_method==2)  {
      ca_result <- cca(as.matrix(df))
      order_sites <- order(scores(ca_result,1,"sites"))
      order_species <- order(scores(ca_result, 1, "species"))
      rva <- as.matrix(df[order_sites,order_species])
    } else {
      rva <- as.matrix(df)
    }
    return(rva)
  })
  
  ca <- eventReactive(input$go, {
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   row.names = if(input$rownames==TRUE) 1 else NULL,
                   sep = input$sep,
                   quote = input$quote)

    rva <- cca(as.matrix(df))
    
    ranges <- lapply(scores(rva),function(x){apply(x,2,range)})
    
    x_range <- range(
      unlist(
        lapply(
          ranges,
          function(x){
            x[,1]
          }
        )
      )
    )
    
    y_range <- range(
      unlist(
        lapply(
          ranges,
          function(x){
            x[,2]
          }
        )
      )
    )
    
    x_range[1] <- floor(x_range[1])
    x_range[2] <- ceiling(x_range[2])
    y_range[1] <- floor(y_range[1])
    y_range[2] <- ceiling(y_range[2])
    
    updateSliderInput(session, "xlim",
                      min=x_range[1],
                      max=x_range[2],
                      value=x_range)
    
    updateSliderInput(session, "ylim",
                      min=y_range[1],
                      max=y_range[2],
                      value=y_range)
    
    return(rva)
  })
  
  plot_ca <- function(ca) {
    ca_scores <- melt(
      scores(ca,
             choices = c(input$x_axis_dimension,
                         input$y_axis_dimension),
             display = input$display_what,
             scaling = 3)
    )

    ca_scores<-spread(ca_scores,
                            key=Var2,
                            value = value)

    test_for_ggplot <- data.frame(name = ca_scores$Var1,
                                  x = ca_scores[,ncol(ca_scores)-1],
                                  y=ca_scores[,ncol(ca_scores)])
    

    test_for_ggplot$type <- if("L1" %in% names(ca_scores))
    {
      ca_scores$L1
    }
    else
    {
      input$display_what
    }

    my_plot <- ggplot(test_for_ggplot)
    
    if("points" %in% input$display_how) {
      my_plot <- my_plot + geom_point(
        aes(x = x,
            y = y,
            col = type)
      )
    }
    if("text" %in% input$display_how) {
      if(input$spead_labels) {
        my_plot <- my_plot + geom_text_repel(
          aes(x = x,
              y = y,
              label = name),
          size=2)
      } else {
        my_plot <- my_plot + geom_text(
          aes( x = x,
               y = y,
               label = name)
          )
      }
    }
    
    x_range <- c(floor(min(test_for_ggplot$x)),
                 ceiling(max(test_for_ggplot$x)))
    
    y_range <- c(floor(min(test_for_ggplot$y)),
                 ceiling(max(test_for_ggplot$y)))
    
    updateSliderInput(session, "xlim",
                      min=x_range[1],
                      max=x_range[2])
    
    updateSliderInput(session, "ylim",
                      min=y_range[1],
                      max=y_range[2])
    
    my_plot <- my_plot + xlim(input$xlim) + ylim(input$ylim)
    
    my_plot
  }
  
  plot_ser <- function(ser) {
    my_ser <- ser()
    
    #cat(file=stderr(), head(my_ser), "\n")
    
    ser.df <- data.frame(melt(my_ser))
    colnames(ser.df) <- c("sites","types","abundance")
    ser.df$sites <- factor(ser.df$sites, levels=unique(ser.df$sites))
    ser.df$types <- factor(ser.df$types, levels=unique(ser.df$types))
    ser.df$abundance <- ifelse(ser.df$abundance==0, NA, ser.df$abundance)
    
    ggplot(data=ser.df, aes(sites,
                            types,
                            size=abundance)) +
      geom_point() +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1))
  }
  
  output$plot_ca <- renderPlot({
    
    if (is.null(ca)) return()
    plot_ca(ca())
  }, height = 800, width = 800)
  
  output$plot_ser <- renderPlot({
    
    if (is.null(ser)) return()
    plot_ser(ser())
  }, height = 800, width = 800)
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste(input$file1$name,
            '.',
            input$download_format, sep='')
      },
    content = function(file) {
      ggsave(file,
             plot = plot_ca(ca()),
             device = input$download_format)
    }
  )
  
  output$downloadSerPlot <- downloadHandler(
    filename = function() {
      paste(input$file1$name,
            '.',
            input$download_format_ser, sep='')
    },
    content = function(file) {
      if(input$download_format_ser=="csv") {
        write.csv(ser(),file,row.names = T)
      } else {
      ggsave(file,
             plot = plot_ser(ser()),
             device = input$download_format_ser,
             width = 40, height = 20, units = "cm")
      }
    }
  )
  
  output$downloadResult <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      test_for_export <- melt(
        scores(ca(),
               choices = c(input$x_axis_dimension,
                           input$y_axis_dimension),
               display = input$display_what,
               scaling = 3)
        )
      test_for_export <- spread(test_for_export,
                                key=Var2,
                                value = value)
      
      test_for_export <- data.frame(name = test_for_export$Var1,
                                    x = test_for_export[,ncol(test_for_export)-1],
                                    y=test_for_export[,ncol(test_for_export)])
      
      write.csv(test_for_export,
                file,
                row.names = FALSE)
    }
  )
}
# Run the app ----
shinyApp(ui, server)
