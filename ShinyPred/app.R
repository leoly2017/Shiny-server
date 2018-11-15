#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Prediction of Virus Epidemic Months"),
   h5("by You Li, You.Li2@ed.ac.uk; last update: 15-Nov-2018"),
   hr(),
   h4("Please input monthly temperature in °C"),
   flowLayout(numericInput(inputId = "temp1", label = "January", value = 4),
              numericInput(inputId = "temp2", label = "February", value = 4),
              numericInput(inputId = "temp3", label = "March", value = 6),
              numericInput(inputId = "temp4", label = "April", value = 9),
              numericInput(inputId = "temp5", label = "May", value = 12),
              numericInput(inputId = "temp6", label = "June", value = 14),
              numericInput(inputId = "temp7", label = "July", value = 16),
              numericInput(inputId = "temp8", label = "August", value = 16),
              numericInput(inputId = "temp9", label = "September", value = 14),
              numericInput(inputId = "temp10", label = "October", value = 10),
              numericInput(inputId = "temp11", label = "November", value = 7),
              numericInput(inputId = "temp12", label = "December", value = 3)
   ),
   h4("Please input monthly relative humidity in %"),
   flowLayout(numericInput(inputId = "rh1", label = "January", value = 87),
              numericInput(inputId = "rh2", label = "February", value = 88),
              numericInput(inputId = "rh3", label = "March", value = 79),
              numericInput(inputId = "rh4", label = "April", value = 77),
              numericInput(inputId = "rh5", label = "May", value = 75),
              numericInput(inputId = "rh6", label = "June", value = 77),
              numericInput(inputId = "rh7", label = "July", value = 77),
              numericInput(inputId = "rh8", label = "August", value = 77),
              numericInput(inputId = "rh9", label = "September", value = 79),
              numericInput(inputId = "rh10", label = "October", value = 84),
              numericInput(inputId = "rh11", label = "November", value = 89),
              numericInput(inputId = "rh12", label = "December", value = 89)
   ),
   actionButton("action", label = "Generate"),
   h4("Results"),
   navlistPanel(
     tabPanel("RSV", 
              tabsetPanel(
                tabPanel("Figure", plotOutput(outputId = "plot.RSV", width = 800)),
                tabPanel("Table", DTOutput('table.RSV'))
              )),
     tabPanel("IFV", 
              tabsetPanel(
                tabPanel("Figure", plotOutput(outputId = "plot.IFV", width = 800)),
                tabPanel("Table", DTOutput('table.IFV'))
              )),
     tabPanel("IFV A", 
              tabsetPanel(
                tabPanel("Figure", plotOutput(outputId = "plot.IFVA", width = 800)),
                tabPanel("Table", DTOutput('table.IFVA'))
              )),
     tabPanel("IFV B", 
              tabsetPanel(
                tabPanel("Figure", plotOutput(outputId = "plot.IFVB", width = 800)),
                tabPanel("Table", DTOutput('table.IFVB'))
              )),
     tabPanel("IFV A(H1N1)pdm", 
              tabsetPanel(
                tabPanel("Figure", plotOutput(outputId = "plot.IFVH1N1pdm", width = 800)),
                tabPanel("Table", DTOutput('table.IFVH1N1pdm'))
              )),
     tabPanel("IFV A(H3N2)", 
              tabsetPanel(
                tabPanel("Figure", plotOutput(outputId = "plot.IFVH3N2", width = 800)),
                tabPanel("Table", DTOutput('table.IFVH3N2'))
              )),
     widths = c(2,5)
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(DT.options = list(pageLength = 12))
  load("IFV.RData")
  load("IFVA.RData")
  load("IFVB.RData")
  load("IFVH1N1pdm.RData")
  load("IFVH3N2.RData")
  load("RSV.RData")
  genMC <- function(inputpred, cut.off = 75) {
    inputpred[25] <-1
    MX <- matrix(
      c(
        rnorm(2000, mean = inputpred[1], sd = inputpred[13] * inputpred[25]),
        rnorm(2000, mean = inputpred[2], sd = inputpred[14]* inputpred[25]),
        rnorm(2000, mean = inputpred[3], sd = inputpred[15]* inputpred[25]),
        rnorm(2000, mean = inputpred[4], sd = inputpred[16]* inputpred[25]),
        rnorm(2000, mean = inputpred[5], sd = inputpred[17]* inputpred[25]),
        rnorm(2000, mean = inputpred[6], sd = inputpred[18]* inputpred[25]),
        rnorm(2000, mean = inputpred[7], sd = inputpred[19]* inputpred[25]),
        rnorm(2000, mean = inputpred[8], sd = inputpred[20]* inputpred[25]),
        rnorm(2000, mean = inputpred[9], sd = inputpred[21]* inputpred[25]),
        rnorm(2000, mean = inputpred[10], sd = inputpred[22]* inputpred[25]),
        rnorm(2000, mean = inputpred[11], sd = inputpred[23]* inputpred[25]),
        rnorm(2000, mean = inputpred[12], sd = inputpred[24]* inputpred[25])
      ),
      byrow = TRUE, nrow = 12
    )
    MXs <- t(t(MX)/apply(MX, 2, sum))*100
    res <- as.data.frame(t(apply(MXs, 1, FUN = quantile, probs = c(0.5, 0.025, 0.975))))
    names(res) <- c("est", "lower", "upper")
    res[res<0] <- 0
    res <- round(res, 2)
    res$Month <- factor(month.abb, levels = month.abb)
    res$epi <- FALSE
    aapOrder <- order(res$est, decreasing = TRUE)
    for (j in 1:12) {
      if(sum(res$est[aapOrder[1:j]]) >= cut.off) {
        res$epi[aapOrder[1:j]] <- TRUE
        break
      }
    }
    return(res)
  }
  getOnsetMonth <- function(epiVector, Offset = FALSE) {
    epiVectorB <- epiVector[c(12, 1:11)]
    onsetMonth <- c(1:12)[epiVector==TRUE & epiVectorB==FALSE]
    epiVectorD <- epiVector[c(1:12, 1:12)]
    epiVectorA <- epiVectorD[(onsetMonth[1]+1):24 ]
    offsetMonth <- (onsetMonth[1]:24)[
      epiVectorD[(onsetMonth[1]):23 ]==TRUE & epiVectorA==FALSE
      ]
    offsetMonth <- offsetMonth[1:length(onsetMonth)]
    seasonLength <- offsetMonth - onsetMonth + 1
    if(Offset){
      offsetMonth[offsetMonth >12] <- offsetMonth[offsetMonth >12] - 12
      return(offsetMonth[order(seasonLength, decreasing = TRUE)[1]])
    }else{
      return(onsetMonth[order(seasonLength, decreasing = TRUE)[1]])
    }
  }
  plotSimu <- function(MCres) {
    if(sum(is.na(MCres))==0){
      colToFill <- c("#FFEDA0", "#FEB24C", "#F03B20")
      OnsetMonth <- getOnsetMonth(MCres$epi)
      MCres$epi <- as.character(MCres$epi)
      MCres$epi[OnsetMonth] <- "Onset"
      MCres$epi[MCres$epi=="TRUE"] <- "Epidemic"
      MCres$epi[MCres$epi=="FALSE"] <- "Non-epidemic"
      MCres$epi <- factor(MCres$epi, levels = c("Non-epidemic",
                                                "Onset",
                                                "Epidemic"))
      resPlot <- ggplot(data = MCres, aes(x = Month)) +
        geom_col(aes(y = MCres$est, fill = MCres$epi),
                 show.legend = TRUE) +
        geom_errorbar(data = MCres, mapping = aes(ymin = lower, ymax = upper),
                      size = 0.5, width = 0.3) +
        labs(y = "Annual Average Percentage Predicted") +
        theme(text = element_text(size = 15)) + 
        scale_fill_manual(name = NULL,
                          values=colToFill)
      return(resPlot) 
    }
  }
  observeEvent(input$action, {
    inputDF <- data.frame(
        TEMPc = c(input$temp1, input$temp2, input$temp3,
                  input$temp4, input$temp5, input$temp6,
                  input$temp7, input$temp8, input$temp9,
                  input$temp10, input$temp11, input$temp12) - mean(
                    c(
                      input$temp1, input$temp2, input$temp3,
                      input$temp4, input$temp5, input$temp6,
                      input$temp7, input$temp8, input$temp9,
                      input$temp10, input$temp11, input$temp12)
                  ),
        RHc = c(input$rh1, input$rh2, input$rh3,
                input$rh4, input$rh5, input$rh6,
                input$rh7, input$rh8, input$rh9,
                input$rh10, input$rh11, input$rh12)- mean(
                  c(
                    input$rh1, input$rh2, input$rh3,
                    input$rh4, input$rh5, input$rh6,
                    input$rh7, input$rh8, input$rh9,
                    input$rh10, input$rh11, input$rh12)
                )
      )
    predEst.RSV <- unlist(predict(fit.RSV, inputDF, se = TRUE))
    predEst.IFV <- unlist(predict(fit.IFV, inputDF, se = TRUE))
    predEst.IFVA <- unlist(predict(fit.IFVA, inputDF, se = TRUE))
    predEst.IFVB <- unlist(predict(fit.IFVB, inputDF, se = TRUE))
    predEst.IFVH1N1pdm <- unlist(predict(fit.IFVAH1N1pdm, inputDF, se = TRUE))
    predEst.IFVH3N2 <- unlist(predict(fit.IFVAH3N2, inputDF, se = TRUE))
    dfRSV <- genMC(predEst.RSV, cut.off = 55)
    dfIFV <- genMC(predEst.IFV, cut.off = 55)
    dfIFVA <- genMC(predEst.IFVA, cut.off = 55)
    dfIFVB <- genMC(predEst.IFVB, cut.off = 60)
    dfIFVH1N1pdm <- genMC(predEst.IFVH1N1pdm, cut.off = 60)
    dfIFVH3N2 <- genMC(predEst.IFVH3N2, cut.off = 60)
    output$table.RSV = renderDT(
      dfRSV[c(4,5, 1:3)],
      colnames = c("Month","Epidemic Month","Median", "Lower", "Upper"),
      extensions = 'Buttons', options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
    )
    output$plot.RSV = renderPlot({
      plotSimu(dfRSV)
    })
    output$table.IFV = renderDT(
      dfIFV[c(4,5, 1:3)],
      colnames = c("Month","Epidemic Month","Median", "Lower", "Upper"),
      extensions = 'Buttons', options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
    )
    output$plot.IFV = renderPlot({
      plotSimu(dfIFV)
    })
    output$table.IFVA = renderDT(
      dfIFVA[c(4,5, 1:3)],
      colnames = c("Month","Epidemic Month","Median", "Lower", "Upper"),
      extensions = 'Buttons', options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
    )
    output$plot.IFVA = renderPlot({
      plotSimu(dfIFVA)
    })
    output$table.IFVB = renderDT(
      dfIFVB[c(4,5, 1:3)],
      colnames = c("Month","Epidemic Month","Median", "Lower", "Upper"),
      extensions = 'Buttons', options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
    )
    output$plot.IFVB = renderPlot({
      plotSimu(dfIFVB)
    })
    output$table.IFVH1N1pdm = renderDT(
      dfIFVH1N1pdm[c(4,5, 1:3)],
      colnames = c("Month","Epidemic Month","Median", "Lower", "Upper"),
      extensions = 'Buttons', options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
    )
    output$plot.IFVH1N1pdm = renderPlot({
      plotSimu(dfIFVH1N1pdm)
    })
    output$table.IFVH3N2 = renderDT(
      dfIFVH3N2[c(4,5, 1:3)],
      colnames = c("Month","Epidemic Month","Median", "Lower", "Upper"),
      extensions = 'Buttons', options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
    )
    output$plot.IFVH3N2 = renderPlot({
      plotSimu(dfIFVH3N2)
    })
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

