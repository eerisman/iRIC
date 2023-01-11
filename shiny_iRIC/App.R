ui <- fluidPage(
  
  titlePanel("Idealized Retention Index Calculations by mrED"),
  
  fluidRow(
    column(12,
           wellPanel(paste("Warning:  This is intended to give an estimate if 
          two compouds are likely to be seperable based on an alkane ladder, 
          retention time of one compound and the retention idex of another.  
          Marjor assumptions include Gaussian peak shapes and accurate RI 
          measurements.  RI will vary (hopefully slightly) with column 
          conditions (age and current state of maintenance of the GC).  RI will 
          vary (hopefully slightly) with column brand (eg HP5, DB5, Rxi5).  
          RI can vary greatly with temperature ramp rate.  For best results use 
          with an alkane ladder containing every alkane (eg not just even or odd
          carbon number) and use cautiously when using RI generated under 
          different conditions.  No accuracy gaurantee is given or implied."))
    )
    
  ),
  fluidRow(
    
    column(4,  
           wellPanel(
             textInput(inputId = "RI1",
                       label = "Alkane 1 Retention Time (min)",
                       value = 9.9824),
             numericInput(inputId = "C1",
                          label = "Carbon count (RI/100)",
                          value = 20,
                          min = 1,
                          step = 1),
             textInput(inputId = "RI2",
                       label = "Alkane 2 Retention Time (min)",
                       value = 10.4516),
             numericInput(inputId = "C2",
                          label = "Carbon count (RI/100)",
                          value = 21,
                          min = 1,
                          step = 1)
           ),
           wellPanel(
             textInput(inputId = "cmpd1",
                       label = "Compound 1 Retention Time (min)",
                       value = 10.165),
             textInput(inputId = "FWHM1",
                       label = "Compound 1 FWHM-Full Width at Half Maximum (min)",
                       value = 0.016),
             textInput(inputId = "cmpd2RI",
                       label = "Compound 2 Retention Index",
                       value = 2116.5)
           )        
    ),
    
    column(8,
           # Output: plot ----
           plotOutput(outputId = "RIPlot"),
           plotOutput(outputId = "CMPDPlot")
           
    )
  ),
  fluidRow(
    column(4),
    column(8,
           wellPanel(
             textOutput("text"),
             textOutput("text1"),
             textOutput("text2")
           )
    )
  )
)

server <- function(input, output, session) {
  
  # RI = lower alkane carbon#*100 + (compound retention time - lower alkane retention time)*(higher alkane carbon#*100-lower alkane carbon#*100)/(higher alkane retention time - lower alkane retention time)
  
  vals <- reactiveValues()
  observe({
    vals$RT1 <- as.numeric(input$RI1)           #alkane 1 retention time
    vals$c1 <- input$C1*100                     #alkane 1 carbon #
    vals$RT2 <- as.numeric(input$RI2)           #alkane 2 retention time
    vals$c2 <- input$C2*100                     #alkane 2 carbon #
    vals$cmpd1 <- as.numeric(input$cmpd1)       #compound 1 retention time
    vals$FWHM <- as.numeric(input$FWHM1)        #compound 1 peak width at 1/2 max
    vals$cmpd2RI <- as.numeric(input$cmpd2RI)   #compound 2 retention index
    
    vals$rt.calc <- (vals$cmpd2RI - vals$c1)/(vals$c2-vals$c1)*(vals$RT2-vals$RT1)+vals$RT1 #compound 2 calculated retention time
  })
  
  output$RIPlot <- renderPlot({
    
    x <- seq(vals$RT1 - 2, vals$RT2 + 2, length=1000)
    y.1 <- dnorm(x, mean = vals$RT1, sd = 0.01) #rounded up sd based on 0.02 min FWHM
    y.1 <- y.1/max(y.1)
    y.2 <- dnorm(x, mean = vals$RT2, sd = 0.01) #rounded up sd based on 0.02 min FWHM
    y.2 <- y.2/max(y.2)
    y.total <- y.1 + y.2
    
    plot(x,y.total, type = "l", xlim = c(vals$RT1 - 1,vals$RT2 +1), main = "Alkanes", xlab = "time (min)", ylab = "abundance")
    text(x=vals$RT1 , y=max(y.total)/2, paste("RI = ",vals$c1), pos = 2)
    text(x=vals$RT1 , y=max(y.total)/2*.9, vals$RT1, pos = 2)
    text(x=vals$RT2 , y=max(y.total)/2, paste("RI = ",vals$c2), pos = 4)
    text(x=vals$RT2 , y=max(y.total)/2*.9, vals$RT2, pos = 4)
    
  })
  output$CMPDPlot <- renderPlot({
    
    x <- seq(min(vals$cmpd1,vals$rt.calc) - 2, max(vals$rt.calc,vals$cmpd1) + 2, length=1000)
    y.1 <- dnorm(x, mean = vals$cmpd1, sd = vals$FWHM/2.355) #FWHM ~2.355 * sigma for Gaussian
    y.1 <- y.1/max(y.1)
    
    y.2 <- dnorm(x, mean = vals$rt.calc, sd = vals$FWHM/2.355) #estimate based on compound
    y.2 <- y.2/max(y.2)
    y.total <- y.1 + y.2
    
    text.pos1 <- if(vals$rt.calc > vals$cmpd1){
      4
    } else{
      2}
    
    text.pos2 <- if(vals$ rt.calc > vals$cmpd1){
      2
    } else{
      4}
    
    plot(x,y.total, type = "l", xlim = c(min(vals$cmpd1,vals$rt.calc) - 1, max(vals$cmpd1,vals$rt.calc) +1), main = "Compounds", xlab = "time (min)", ylab = "abundance")
    text(x=vals$cmpd1 , y=max(y.total)/2, "Compound 1", pos = text.pos2)
    text(x=vals$cmpd1 , y=max(y.total)/2*.9, vals$cmpd1, pos = text.pos2)
    text(x=vals$rt.calc , y=max(y.total)/2, "Compound 2", pos = text.pos1, col = "red")
    text(x=vals$rt.calc , y=max(y.total)/2*.9, "Calculated retention time", pos = text.pos1, col = "red")
    text(x=vals$rt.calc , y=max(y.total)/2*.8, vals$rt.calc, pos = text.pos1, col = "red")
  })
  output$text <- renderText({
    
    cmpd.res <- round(abs(vals$cmpd1- vals$rt.calc)*1.18/(2*vals$FWHM),2)
    noquote(c("Resolution (compound 1, compound 2) = ", cmpd.res))
    
  })
  output$text1 <- renderText({
    
    min.rt <- round(1.5*2*vals$FWHM/1.18,2)
    noquote(c("Approximate minimum retention difference for approximate baseline resolution (R = 1.5):", min.rt, "min"))
  })
  output$text2 <- renderText({
    
    min.rt <- (1.5*2*vals$FWHM/1.18)
    
    min.ri <- round(min.rt * (vals$c2 - vals$c1) /(vals$RT2-vals$RT1),2)
    noquote(c("Approximate minimun RI difference for approximate baseline resolution (R= 1.5) = ", min.ri))
    
  })
  
  session$onSessionEnded(function() {
       stopApp()
  }) 
  
}

shinyApp(ui = ui, server = server)
