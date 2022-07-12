library(shiny)
library(shinythemes)
library(shinydashboard)
library(plotly)

shinyUI(navbarPage(title = "Yale COVID-19 Genomic Surveillance Initiative", theme = shinytheme("sandstone"),

         
    tabPanel("Variant Report",
      fluidPage(
        titlePanel(div(windowTitle = "CovidTrackerCT", img(src = "cover.png", width = "100%", class = "bg"),)),
        
             fluidRow(
               HTML("<p style='font-size:38px; text-align: center'>",
                    "Connecticut SARS-CoV-2 variant surveillance – Report 2022-06-23","</p>"), ## CHANGE TO CURRENT DATE
              
               hr(),
               hr(),
               column(width = 5, offset = 1,
                             h4("In collaboration with the Connecticut Department of Public Health, 
                                the Centers for Disease Control and Prevention, Yale University, Jackson Laboratories, and diagnostic laboratories across the state, 
                                we are conducting surveillance for SARS-CoV-2 variants using viral sequencing and phylogenetics.")
                             ),
               column(width = 3, offset = 2, style = "background-color:#fed76e;",
                      infoBoxOutput("total_seq", width = 10)),
                      ),
        hr(),

        DT::dataTableOutput("table1"),
        HTML("<p style='font-size: 12px; text-align: left'>",
             "*Cumulative sequenced cases = the total cases for each variant that are confirmed by sequencing and reported by the Connecticut Department of Public Health. ","</p>"), 
        HTML("<p style='font-size: 12px; text-align: left'>",
             "**Percent sequenced from past 3 weeks = samples collected within three weeks of the report date and sequenced by a collaborating lab for unbiased SARS-CoV-2 surveillance. ","</p>"), 
        HTML("<p style='font-size: 12px; text-align: left'>",
             "Variant of concern (VOC) = A variant for which there is evidence of an increase in transmissibility, more severe disease (increased hospitalizations or deaths), significant reduction in neutralization by antibodies generated during previous infection or vaccination, reduced effectiveness of treatments or vaccines, or diagnostic detection failures.","</p>"), 
        HTML("<p style='font-size: 12px; text-align: left'>",
             "Variant being monitored (VBM) = A variant for which there are data indicating a potential or clear impact on approved or authorized medical countermeasures or that has been associated with more severe disease or increased transmission but are no longer detected or are circulating at very low levels in the United States, and as such, do not pose a significant and imminent risk to public health in the United States.","</p>"), 
        
        hr(),
        plotlyOutput("lin_freq"),
        hr(),
        plotlyOutput("var_freq"),
        hr(),
        plotOutput("subsampler"),
        hr(),
        fluidRow(
          column(width = 8,plotOutput("CT"),
                 infoBoxOutput("Ct_caption", width = 8)
                 ),
          column(width = 4,img(src = "ctvisual.png", height = "50%",width = "50%", class = "bg"))
          ), # end row
        hr(),
        
        fluidRow(
          HTML("<p style='font-size:30px'>",
               "Additional resources","</p>"), 
          column(width = 6, 
                 
                 h4("We created a custom Nextstrain page to visualize the relatedness of sequenced COVID-19 cases in the state. Below shows a phylogenetic tree of historic sequences in Connecticut, colored by pango lineage assignments. "),
                 htmlOutput("nextstrain")
                           ),
          
          column(width = 6,
                 tags$head(
                   tags$script(async = NA, src = "https://platform.twitter.com/widgets.js")
                 ),
                 h4("For general CovidTrackerCT website inquiries, please email Kien Pham (k.pham@yale.edu)."),
                 h4("Follow us on Twitter @CovidCT"),
                 uiOutput("tweet")
          )
        ),
        hr(),
        HTML("<p style='font-size:30px'>",
             "Acknowledgement","</p>"),
        h4("The data were retrieved from GISAID. We gratefully acknowledge the Authors form Originating and Submitting laboratories of sequence data on which the analysis is based.")
        
        )) # end tab
      )
)
