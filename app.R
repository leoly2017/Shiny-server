#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Supplementary Animated Figures"),
   h3("Global patterns in monthly activity of influenza virus, respiratory syncytial virus, parainfluenza virus and metapneumovirus: a systematic analysis"),
   hr(),
   selectInput("gifChoosing", "Select a virus: ", choices = c("Influenza Virus (IFV)",
                                                              "Respiratory Syncytial Virus (RSV)",
                                                              "Parainfluenza Virus (PIV)",
                                                              "Metapnuemovirus (MPV)"),
               selected = NULL),
   hr(),
   conditionalPanel(condition = "input.gifChoosing=='Influenza Virus (IFV)'",
                    img(src = "IFV.gif")),
   conditionalPanel(condition = "input.gifChoosing=='Respiratory Syncytial Virus (RSV)'",
                    img(src = "RSV.gif")),
   conditionalPanel(condition = "input.gifChoosing=='Parainfluenza Virus (PIV)'",
                    img(src = "PIV.gif")),
   conditionalPanel(condition = "input.gifChoosing=='Metapnuemovirus (MPV)'",
                    img(src = "MPV.gif")),
   HTML("<br>AAP = annual average percentage")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
}

# Run the application 
shinyApp(ui = ui, server = server)

