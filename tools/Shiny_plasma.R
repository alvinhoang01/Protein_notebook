rm(list=ls())
library(shiny)
library(readr)
library(DT)
library(ggplot2)
library(tidyverse)
setwd("/Users/pan/Desktop/2023postdoc/Work/Plasma/Three kits/Shiny")
data <- read.csv("Protein_filter.csv")

# Gather data into long format and add Method column
data_long <- data %>%
  pivot_longer(cols = -Protein, names_to = "key", values_to = "value") %>%
  mutate(Method = case_when(
    str_detect(key, "Preomics") ~ "Preomics",
    str_detect(key, "Thermokit") ~ "Thermo",
    str_detect(key, "Direct") ~ "Direct",
    str_detect(key, "Seer_NPA") ~ "Seer_NPA",
    str_detect(key, "Seer_NPB") ~ "Seer_NPB",
    TRUE ~ NA_character_
  )) %>%
  mutate(Replicate = case_when(
    str_detect(key, "R1") ~ "R1",
    str_detect(key, "R2") ~ "R2",
    str_detect(key, "R3") ~ "R3",
    str_detect(key, "R4") ~ "R4",
    TRUE ~ NA_character_
  ))
# CV#
data_cv <- data_long %>%
  group_by(Protein, Method) %>%
  summarize(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            CV = (sd_value / mean_value) * 100) %>%
  ungroup()
data_long <- left_join(data_long, data_cv %>% select(Protein, Method, CV), by = c("Protein", "Method"))
# UI
ui <- fluidPage(
  titlePanel("Plasma Proteomics Data Viewer"),
  sidebarLayout(
    sidebarPanel(
      textInput("protein_search", "Enter PG.ProteinAccessions", value = ""),
      selectInput("protein_select", "Select PG.ProteinAccessions", 
                  choices = unique(data_long$Protein)),
      actionButton("search", "Search")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plots",
                 fluidRow(
                   column(6, plotOutput("intensityPlot")),
                   column(6, plotOutput("cvPlot"))
                 )
        ),
        tabPanel("Table",
                 DTOutput("dataTable"))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  observeEvent(input$search, {
    protein <- ifelse(input$protein_search != "", 
                      input$protein_search, 
                      input$protein_select)
    filtered_data <- data_long %>% filter(Protein == protein)
    filtered_data_cv <- data_cv %>% filter(Protein == protein)
    
    output$intensityPlot <- renderPlot({
      ggplot(filtered_data, aes(x = Method, y = value, fill = Method)) +
        geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7, width = 0.7) +
        geom_point(aes(color = Replicate), position = position_jitter(width = 0.2, height = 0), size = 2) +
        theme_minimal() +
        labs(title = paste("MS Intensity for", protein),
             x = "Method",
             y = "MS Intensity") +
        scale_fill_manual(values = c("Preomics" = "#b5b59d", "Thermo" = "#427277", "Direct" = "#42293b", "Seer_NPA" = "#87CEFF", "Seer_NPB" = "#009ACD")) +
        scale_color_manual(values = c("R1" = "black", "R2" = "grey", "R3" = "brown", "R4" = "pink"))
    })

    output$cvPlot <- renderPlot({
      ggplot(filtered_data_cv, aes(x = Method, y = CV, fill = Method)) +
        geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7, width = 0.7) +
        theme_minimal() +
        labs(title = paste("CV for", protein),
             x = "Method",
             y = "Coefficient of Variation (CV)") +
        scale_fill_manual(values = c("Preomics" = "#b5b59d", "Thermo" = "#427277", "Direct" = "#42293b", "Seer_NPA" = "#87CEFF", "Seer_NPB" = "#009ACD"))
    })
  })
  
  output$dataTable <- renderDT({
    data_long
  })
}

# Shiny App
shinyApp(ui = ui, server = server)
###for publish#

install.packages("rsconnect")
library(rsconnect)
setwd("/Users/pan/Desktop/2023postdoc/Work/Plasma/Three kits/Shiny")
rsconnect::setAccountInfo(name='myshinyzhenyuproteomics', token='7C4FEDD141AB820C44D727F5EACFEFC9', secret='8A9vp5KP4EeZBOU/38z6TJ9QqBF/DDr98kKId8OG')
rsconnect::deployApp()
