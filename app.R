library(shiny)
library(vegan)
library(ggplot2)
library(reshape2)
library(tidyr)
library(ggrepel)

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
        
        # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),
        
        # Input: Checkbox if file has header ----
        checkboxInput("rownames", "Rownames", TRUE),
        
        # Input: Select separator ----
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        # Input: Select quotes ----
        radioButtons("quote", "Quote",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = '"'),
        
        # Horizontal line ----
        tags$hr(),
        
        # Input: Select number of rows to display ----
        radioButtons("disp", "Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head"),
        actionButton(inputId = "go", label = "Go")
        
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        
        # Output: Data file ----
        tableOutput("contents")
      )
      
    )
  )
}

analysis_ui <- function(x) {
  sidebarLayout(
    sidebarPanel(
      numericInput("x_axis_dimension",
                   "Dimension x-axis",
                   value = 1,
                   min=1),
      numericInput("y_axis_dimension",
                   "Dimension y-axis",
                   value = 2,
                   min=1),
      checkboxGroupInput("display_what",
                         "What to display",
                         choices = list("Sites"="sites",
                                        "Types"="species"),
                         selected = c("sites",
                                      "species")),
      checkboxGroupInput("display_how",
                         "Plot Elements",
                         choices = list("Points"="points",
                                        "Text"="text"),
                         selected = c("points",
                                      "text")),
      checkboxInput("spead_labels",
                    "Spread labels")
    ),
    mainPanel(
      
      # Output: Data file ----
      plotOutput("plot"),
      
      # Input: Choose dataset ----
      selectInput("download_format",
                  "Download as:",
                  choices = c("svg",
                              "png",
                              "pdf")),
      
      # Button
      downloadButton("downloadPlot",
                     "Download")
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
                         tabPanel("Analysis",
                                  id="analysis",
                                  analysis_ui()),
                         tabPanel("Export",
                                  download_ui())
))

# Define server logic to read selected file ----
server <- function(input, output, session) {
  observeEvent(input$go, {
    updateTabsetPanel(session = session,
                      inputId = "navbar",
                      selected = "Analysis")
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
  
  ca <- eventReactive(input$go, {
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    cca(as.matrix(df))
  })
  
  plot_ca <- function(ca) {
    test_for_ggplot <- melt(
      scores(ca,
             choices = c(input$x_axis_dimension,
                         input$y_axis_dimension),
             display = input$display_what,
             scaling = 3)
    )
    
    test_for_ggplot<-spread(test_for_ggplot,
                            key=Var2,
                            value = value)
    
    test_for_ggplot <- data.frame(name = test_for_ggplot$Var1,
                                  x = test_for_ggplot[,ncol(test_for_ggplot)-1],
                                  y=test_for_ggplot[,ncol(test_for_ggplot)])
    
    test_for_ggplot$type <- if("L1" %in% names(test_for_ggplot))
    {
      test_for_ggplot$L1
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
    
    my_plot
  }
  output$plot <- renderPlot({
    
    if (is.null(ca)) return()
    plot_ca(ca())
  })
  
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
