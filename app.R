# Load packages ----
set.seed(1)
library(shiny)
library(tidyverse)

# Load data ----
gtex_hr <- read_tsv('data/GTEx_HR_combined.tsv') # GTEx HR表达数据
df <- readxl::read_excel('data/hormone_receptor_list.xlsx') # 读取激素和HR的信息
df$Hormone_name1 <- recode(df$Hormone_name1, 'orphan'='Unknown')
df$Hormone_chemical_classes <- recode(df$Hormone_chemical_classes, 'NA'='Unknown')
df$Hormone1_tissue <- recode(df$Hormone1_tissue, 'NA'='Unknown')
df <- df %>% separate_rows(Hormone1_tissue, sep = ',')
df <- df %>% filter(Hormone1_tissue != 'Unknown') # 分析组织通讯关系时候不考虑孤儿受体
vec_receptor <- df$Gene_symbol
vec_hormone <- sort(df$Hormone_name1 %>% unique())
vec_tissue <- sort(gtex_hr$SMTS %>% unique())

# Source helper functions -----
source("sankey.R")

# Some settings
options(dplyr.summarise.inform = FALSE)

# User interface ----
ui <- fluidPage(
  # titlePanel("HTHC (human tissue hormone communication)"),
  
  br(),
  br(),
  
  sidebarLayout(
    sidebarPanel(
      # helpText("Explore human tissue communication at hormone-hormone receptor level"),
      
      #设置一个重置按钮，一键重置默认值#
        
      radioButtons("sex",
                   label = "Choose sex to display",
                   choices = c("Male", "Female"),
                   selected = "Male"),
      
      selectInput("hormone_class", 
                  label = "Hormone class",
                  choices = c("All", "Amino acid derivative", "Lipid derivative", "Peptide"),
                  selected = "All"),
      
      selectInput('hormone_name',
                  label = 'Hormone name',
                  choices = c('All', vec_hormone),
                  selected='ALl'),
      
      selectInput('hormone_tissue',
                  label = 'Outgoing tissue (hormone tissue)',
                  choices = c('All', vec_tissue),
                  selected='ALl'),
      hr(style="border-color: black;"),
      
      selectInput("receptor_class", 
                  label = "Receptor class",
                  choices = c("All", "Membrane", "Nucleus"),
                  selected = "All"),
      
      selectInput('receptor_name',
                  label = 'Receptor name (HGNC gene symbol)',
                  choices = c('All', vec_receptor),
                  selected='ALl'),
      
      selectInput('receptor_tissue',
                  label = 'Incoming tissue (receptor tissue)',
                  choices = c('All', vec_tissue),
                  selected='Brain'),
      hr(style="border-color: black;"),
      
      sliderInput("number_HRtissue", 
                  label = "Number of tissues expressing a hormone receptor:",
                  min = 1, max = 10, value = 3),
      
      sliderInput("tpm_log1p", 
                  label = "Receptor expression range (log1p TPM)",
                  min = 0, max = 6.353414, value = c(0.1, 6.353414))
    ),
    
    mainPanel(img(src = "logo.png", width="20%", height="20%", style="display: block; margin-left: auto; margin-right: auto;"),
              h1("Explore human tissue communication at hormone-hormone receptor level",
                 align='center',
                 style = "font-size:15px;"),
              hr(style="border-color: black;"),
              plotOutput("plot")
              )
  )
)


# Server logic ----
server <- function(input, output) {
  output$plot <- renderPlot({
    plot_sankey(input$sex, input$hormone_class, input$hormone_name, input$hormone_tissue,
                input$receptor_class, input$receptor_name, input$receptor_tissue,
                input$number_HRtissue, input$tpm_log1p[1], input$tpm_log1p[2],
                gtex_hr, df, vec_receptor)
  })
}

# Run app ----
shinyApp(ui, server)
