library(shiny)
library(shinythemes)
library(shinydashboard)

shinyUI(navbarPage(title = "Yale COVID-19 Genomic Surveillance Initiative", theme = shinytheme("sandstone"),

         
    tabPanel("Variant Report 2022: Last updated 05-26-2022",
      fluidPage(
        titlePanel(div(windowTitle = "CovidTrackerCT", img(src = "cover.png", width = "100%", class = "bg"),)),
        
             fluidRow(
                      column(width = 4, 
                             HTML("<p style='font-size:35px'>",
                                  "SARS-CoV-2 Variant surveillance","</p>"),
                             
                             h4("We’re using viral sequencing and phylogenetics to learn more about how SARS-CoV-2 is spreading in Connecticut."),
                             h4("Our group, in collaboration with the Connecticut Department of Health, has intensified the genomic surveillance of variants of SARS-CoV-2 in our region.")
                             ),
                      column(width = 4, 
                             infoBoxOutput("total_seq", width = 5)),
                      column(width = 4,
                             plotOutput("legend"))
                      ),
             
          plotOutput("lin_freq"),

         DT::dataTableOutput("table_1"),
         h5("Author: Anne Hane, Kien Pham")
         )
      ),
    
    tabPanel("PCR Trends",
             h5("Author: Anne Hane, Kien Pham"),
             fluidRow(column(width = 6,plotOutput("sample_plot")),
                      column(width = 5, offset = 1,
                             h3("PCR Trends"),
                             h4("Every week, we process and submit samples, sequenced by Yale Center for Genome Analysis (YCGA), to the public database GISAID."),
                             h4("All samples and monitored for notable mutations in the S-gene and ORF-gene targets that at times are indicative of novel variants"))
                      ),
             plotOutput("trendplot")
             ),

    tabPanel("Subsampler",
             h5("Author: Anne Hane, Kien Pham"),
             fluidRow(column(width = 6,plotOutput("weekly_total")),
             column(width = 5, offset = 1,
                    h3("% Sequenced"),
                    h4("The graph displays the overall COVID-19 case count in Connecticut, provided by John Hopkins University, as well as the proportion of such cases sequenced for genomics by the Yale COVID-19 Genomic Surveillance Initiative")),
             plotOutput("subsampler")
                  )
    ),
    
    tabPanel("Contact",
             tags$head(
               tags$script(async = NA, src = "https://platform.twitter.com/widgets.js")
             ),
             h3("For general CovidTrackerCT website inquiries, please email Chaney Kalinich (chaney.kalinich@yale.edu) or Kien Pham (k.pham@yale.edu)."),
             h3("Follow us on Twitter @CovidCT"),
             uiOutput("tweet")
    ))
    
)