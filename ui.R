#### Variant Report for SARS-CoV-2 ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#initial commit date: 06/25/22
#original author: Kien Pham
#email k.pham@yale.edu or duckien242@gmail.com
#to be updated Tuesday afternoons

title = "Connecticut SARS-CoV-2 variant surveillance – Report 2022-07-13" ## CHANGE TO CURRENT DATE OF REPORTING

########################################################
library(shiny)
library(shinythemes)
library(shinydashboard)
library(plotly)

shinyUI(navbarPage(title = "Yale COVID-19 Genomic Surveillance Initiative", theme = shinytheme("sandstone"),

         
    tabPanel("Variant Report",
      fluidPage(
        titlePanel(div(windowTitle = "CovidTrackerCT", img(src = "cover.png", width = "100%", class = "bg"),)),
        
             fluidRow(
               HTML("<p style='font-size:38px; text-align: center'>", title,
                    "</p>"), 
              
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
                 infoBoxOutput("Ct_caption", width = 12)
                 ),
          column(width = 4,img(src = "ctvisual.png", height = "50%",width = "50%", class = "bg"))
          ), # end row
        hr(),
        
        fluidRow(
          HTML("<p style='font-size:30px'>",
               "Additional resources","</p>"), 
          column(width = 6, 
                 
                  infoBoxOutput("nextstrain_caption", width = 12),
                 
                 htmlOutput("nextstrain")
                           ),
          
          column(width = 5, offset = 1,
                 
                 tags$head(tags$script('!function(d,s,id){var js,fjs=d.getElementsByTagName(s)    [0],p=/^http:/.test(d.location)?\'http\':\'https\';if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+"://platform.twitter.com/widgets.js";fjs.parentNode.insertBefore(js,fjs);}}(document,"script","twitter-wjs");'))
                 ,
                 infoBoxOutput("twitter_caption", width = 10),

                 box(style='width:600px;overflow-x: scroll;height:800px;overflow-y: scroll;',
                     a("Follow us on Twitter @CovidCT", class="twitter-timeline", href = "https://twitter.com/CovidCT")
                 )
         )
        ),
        hr(),
        HTML("<p style='font-size:30px'>",
             "Acknowledgement","</p>"),
        h4("The data were retrieved from GISAID. We gratefully acknowledge the Authors form Originating and Submitting laboratories of sequence data on which the analysis is based.")
        ,
        hr(),
        HTML("<p style='font-size:30px'>",
             "Contact","</p>"),
        h4("For media inquiries, please email Anne Hahn (anne.hahn@yale.edu) and/or Nathan Grubaugh (nathan.grubaugh@yale.edu). For general CovidTrackerCT website inquiries, please email Kien Pham (k.pham@yale.edu).")
        
        )), # end tab
  
  tabPanel("Why genomics",
           fluidPage(
             titlePanel(div(windowTitle = "CovidTrackerCT", img(src = "cover.png", width = "100%", class = "bg"),)),
             HTML("<p style='font-size:38px; text-align: center'>",
                  "Why genomics?","</p>"),
             hr(),
             hr(),
     fluidRow(
               
               column(width = 5, offset = 1,
                      h4("By sequencing many cases of SARS-CoV-2, the virus that causes COVID-19, we can learn about how, where, and when it is transmitted. Sequencing is a process used to find out the order of DNA bases in a genome. 
                      All viruses (and every living thing, for that matter) makes small errors when replicating their genomes. While these errors (called “mutations”) are generally harmless, future infections will all end up with the same error. 
                         As the virus accumulates these errors, we can sequence them and use common errors to link cases together.")
               ),
               column(width = 4, offset = 1, style = "background-color:#fed76e;",
                      HTML('<iframe width="620" height="400" src="https://www.youtube.com/embed/UTzUtW3qs7M" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>')
                      ),
             
             ),
     hr(),
     navlistPanel(
              tabPanel("What we can (and cannot) learn",
                       HTML("<p style='font-size:18px'>","Genomics can be used to answer a wide variety of questions. In fact, the novel coronavirus was first discovered when a large group of people suddenly became sick with a strange pneumonia, and researchers at the Shanghai Public Health Clinical Center & School of Public Health sequenced samples from the patients to identify what was making them sick.  
                          That genome is still used by many others as a reference. We’re using genomics to answer a few different questions related to transmission in Connecticut (and surrounding states), such as:"),
                       HTML("<p style='font-size:22px'>",
                            "<ul><li>Is SARS-CoV-2 frequently imported from other outbreaks, or are most cases related to SARS-CoV-2 already circulating in Connecticut?
                            </li><li>When did the epidemic begin in Connecticut?
                            </li><li>Where was the virus introduced from?
                            </li><li>What factors are contributing to spread within Connecticut?</li></ul>"),
                       HTML("<p style='font-size:18px'>","<b>While genomics are a powerful tool for understanding the SARS-CoV-2 pandemic, there are some important caveats to consider.</b>"),
                       HTML("<p style='font-size:18px'>","First, and most important, the sequences generated are not a random sample of cases and are not representative of all of the outbreaks. The virus will only be sequenced from a small fraction of COVID-19 cases. These cases are typically closer to research or public health institutions with sequencing capabilities. This means that, while we may make general inferences about geographic spreading patterns, we can’t draw exact conclusions about where a virus came from without data from many different places. For example, that if the closest genetic relative of a SARS-CoV-2 virus sequenced in location A is an earlier sequence from in location B, we know there were lots of cases in between. 
                            The virus may have traveled from B to C, C to D, and finally D to A."),
                       img(src = "learn_genomics.png"),
                       HTML("<p style='font-size:18px'>","In addition to this, sequencing isn’t perfect. It is standard to go back and sequence samples again if we didn’t get the entire genome in one go, so sequences will be updated. 
                       <strong> Data on this site should be considered preliminary.</strong> 
                       Because our focus is on Connecticut, we don’t include all of the sequences that have been generated worldwide. 
                            We choose representative samples of outbreaks in other countries and regions in the US in order to figure out how often SARS-CoV-2 is being introduced, but we can’t draw any conclusions regarding transmission outside of Connecticut with this data."),
                       
                       ),
              tabPanel("How does it work?",
                       HTML("<p style='font-size:18px'>","The exact steps we will take for SARS-CoV-2 genomics might vary a little based on the question we want to answer, but follows the same general steps. 
                            The approach we are using to sequence is based on the <a href='https://artic.network/protocols.html'>ARTIC Network protocol</a> 
                            and is similar to what was used for real-time sequencing during the 2013-2016 Ebola virus epidemic. 
                            If you want more details, check out <a href = 'https://covidtrackerct.com/protocols/'>our protocols </a>.
                            We were able to generate a sequence of the first detected case in Connecticut within 14 hours, and plan to release more sequences and a discussion of what they mean on a weekly basis moving forward."),
                       HTML("<p style='font-size:18px'>","Once we sequence samples, we build an evolutionary tree.  An evolutionary tree, or phylogenetic tree, is a branching diagram that shows the evolutionary relationship between sequences. Because the virus mutates at a fairly standard rate, these mutations can serve as a “molecular clock” where the number of differences between sequence samples correspond to the amount of time since they had a common ancestor. This is often referred to as the most recent common ancestor or MRCA. The most recent common ancestor, which we estimate with our data, can tell us about when and where a virus was introduced."),
                       HTML("<p style='font-size:18px'>","We can also make a map showing how the virus spreads geographically. Because sequencing across the globe has not been done at the same amount, we can only make general conclusions about where the virus is coming from outside of Connecticut, like we did in our paper on the 
                            <a href = 'https://www.sciencedirect.com/science/article/pii/S0092867420304840'>Coast-to-coast spread of SARS-CoV-2 during the early epidemic in the United States </a>.
                            However, because we are sampling many viruses throughout the state, we will be able to learn about how it’s spreading in Connecticut and factors that are increasing or preventing transmission."),
              ),
              tabPanel("Why does this matter?",
                       HTML("<p style='font-size:18px'>","Tracking the spread of SARS-CoV-2 in Connecticut can allow for more educated decisions on how to control its spread. 
                       The geographic disease transmission patterns can tell us whether interventions like border closures and decreasing daily work commutes are working, and provide insight into how the virus is still spreading. In addition, as the pandemic declines, the routine evaluation of genome sequences can provide insight into whether new introductions are still occurring or if new cases are still the result of local transmission in the community, which could tell us whether it is safe to start lifting restrictions."),
              )
      ) 
    )       
  ), # end tab
  
  tabPanel("About",
           
           fluidPage(
             titlePanel(div(windowTitle = "CovidTrackerCT", img(src = "cover.png", width = "100%", class = "bg"),)),
             
             fluidRow(
               HTML("<p style='font-size:38px; text-align: center'>",
                    "About","</p>"),
               column(width = 4, 
                      img(src = "covidtracker.png", height="40%", width="40%", align="center"),
                      infoBoxOutput("about_tracker", width = 12)
               ),
               column(width = 4, style = "background-color:#fed76e;",
                      img(src = "glab_logo.png", height="50%", width="50%", align="center"),
                      infoBoxOutput("about_grubaughlab", width = 12)
               ),
               
               column(width = 4, 
                      img(src = "impact.png", height="50%", width="50%", align="center"),
                      infoBoxOutput("about_impact", width = 12)
             ),
           ),
           hr(),
           HTML("<p style='font-size:30px'>",
                "Contact","</p>"),
           h4("For media inquiries, please email Anne Hahn (anne.hahn@yale.edu) and/or Nathan Grubaugh (nathan.grubaugh@yale.edu). For general CovidTrackerCT website inquiries, please email Kien Pham (k.pham@yale.edu).")
           
      
  )
  ) # end tab

  
      )
)



