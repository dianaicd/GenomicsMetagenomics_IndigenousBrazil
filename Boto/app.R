#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(cowplot)
library(ggplot2)
source("~/Projects/Botocudos/Scripts/translate_ids.R")
source("~/Projects/Botocudos/Scripts/length_distribution_plot.R")
source("~/Projects/Botocudos/Scripts/preseq_only_plots.R")
boto <- read.table("~/Projects/Botocudos/Files/Summaries/Botocudos_summary_2017_09_18.table", header = T)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = "lumen",
                # Application title
                titlePanel("Botocudos"),
                
                sidebarPanel(
                  #   tags$head(
                  #   tags$style(type="text/css", "select { max-width: 50px; }")
                  # ),
                  checkboxGroupInput("individual", label = "Individual:", 
                                     choices = c("All", 
                                                 as.vector(
                                                   unique(
                                                     boto$Library[boto$Library != "*"]))), 
                                     selected = "All")
                ),
                mainPanel(#width = 10,
                  plotOutput("length_plot")
                )
)  
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$length_plot <- renderPlot({
    #plots <- list()
    plots <- reactive(lapply(as.vector(unique(boto$Library[boto$Library != "*"])), 
                             function(ind) gglength(ind)))
    
    if(input$individual[1] == "All"){
      # for(ind in (unique(boto$Library[boto$Library != "*"]))){
      #   p <- gglength(ind)
      #   plots[[length(plots)+1]] <- p
      # }
      plot_grid(plotlist = plots(), ncol = 5)
    }else{
      # for(ind in input$individual){
      #   p <- gglength(ind)
      #   plots[[length(plots)+1]] <- p
      # }
      num_col <- round(sqrt(length(input$individual)))
      plots <- lapply(as.vector(input$individuals), function(ind) gglength(ind))
      plot_grid(plotlist = plots, ncol = num_col)
    }
    
  }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

