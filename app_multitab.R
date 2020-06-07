library(shiny)
library(stringr)
library(MASS)
library(splines)
library(survival)
library(ggplot2)
source("video_code.R")
#########################
ui <- fluidPage(
  
  titlePanel("Enhanced Kaplan-Meier Curves"),
  tabsetPanel(
    tabPanel("Arm 1", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(textAreaInput("arm_1_s", "Paste from WebPlotDigitizer"),
                            textAreaInput("arm_1_interval", "Interval length"),
                            textAreaInput("arm_1_n_at_risk", "Number at risk (comma-separated)"),
                            actionButton("do", "Run")),
               mainPanel(plotOutput("plot"),
                         numericInput("quantile", "Quantile", value = 0.5, min = 0, max = 1, step = 0.1),
                         textOutput("info_q"),
                         numericInput("prob", "Milestone", value = 1, min = 0),
                         textOutput("info_p"),
                         downloadButton('download_1', 'Download Arm 1 IPD'))
             )
             
    ),
    tabPanel("Arm 2", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(textAreaInput("arm_2_s", "Paste from WebPlotDigitizer"),
                            textAreaInput("arm_2_interval", "Interval length"),
                            textAreaInput("arm_2_n_at_risk", "Number at risk (comma-separated)"),
                            actionButton("do_2", "Run")),
               mainPanel(plotOutput("plot_2"),
                         numericInput("quantile_2", "Quantile", value = 0.5, min = 0, max = 1, step = 0.1),
                         textOutput("info_q_2"),
                         numericInput("prob_2", "Milestone", value = 1, min = 0),
                         textOutput("info_p_2"),
                         downloadButton('download_2', 'Download Arm 2 IPD'))
             )
             
    ),
    tabPanel("Both arms", fluid = TRUE,
             plotOutput("plot_both"),
             verbatimTextOutput("cox_text"),
             downloadButton('download_both', 'Download Both Arms IPD'))
  )
  
)
server <- function(input, output) {
  
  ########## Arm 1 ###################
  observeEvent(input$do, {
    wpd <- reactive({str_split(str_split(input$arm_1_s, "\n", simplify = TRUE),
                               ",", simplify = TRUE)})
    wpd_t <- reactive({as.numeric(wpd()[,1])})
    wpd_s <- reactive({as.numeric(wpd()[,2])})
    
    nrn <- reactive({as.numeric(str_split(input$arm_1_n_at_risk,",")[[1]])})
    trn <- reactive({as.numeric(input$arm_1_interval) * (seq_along(nrn()) - 1)})
    
    IPD <- reactive({getIPD(wpd_t(),wpd_s(),trn(),nrn())})
    
    IPD_survfit <- reactive({survfit(Surv(IPD()[,1], IPD()[,2])~1)})
    
    my_quantile <- reactive({round(unlist(quantile(IPD_survfit(), probs = input$quantile)),2)})
    sum_IPD <- reactive({summary(IPD_survfit(), time = input$prob)})
    
    #output$value <- renderText({ input$arm_1_s })
    #output$value2 <- renderText({ input$arm_1_interval })
    #output$value3 <- renderText({ input$arm_1_n_at_risk })
    output$plot <- renderPlot({ 
      
      #plot(survfit(Surv(IPD()[,1], IPD()[,2])~1)) })
      
      p1 <- survminer::ggsurvplot(IPD_survfit(), 
                                  data = as.data.frame(IPD()), 
                                  risk.table = TRUE, 
                                  legend.labs = "Arm 1",
                                  palette = "red",
                                  break.x.by = as.numeric(input$arm_1_interval)) 
      p1$plot <- p1$plot + 
        geom_hline(yintercept = 1 - input$quantile) +
        geom_vline(xintercept = input$prob, linetype = "dashed") 
      
      p1
    })
    
    
    output$info_q <- renderText({paste0(my_quantile()[1],
                                        " with 95% conf int (",
                                        my_quantile()[2], 
                                        ", ", 
                                        my_quantile()[3], 
                                        ")")})
    
    output$info_p <- renderText({paste0(round(sum_IPD()$surv, 2),
                                        " with 95% conf int (",
                                        round(sum_IPD()$lower, 2), 
                                        ", ", 
                                        round(sum_IPD()$upper, 2), 
                                        ")")})
    
    output$download_1 <- downloadHandler(
         filename = function() {
           "IPD_1.csv"
         },
         content = function(con) {
           write.csv(IPD(), con, row.names = FALSE)
         }
       )
    
    ########## Arm 2 ###################
    observeEvent(input$do_2, {
      wpd_2 <- reactive({str_split(str_split(input$arm_2_s, "\n", simplify = TRUE),
                                   ",", simplify = TRUE)})
      wpd_t_2 <- reactive({as.numeric(wpd_2()[,1])})
      wpd_s_2 <- reactive({as.numeric(wpd_2()[,2])})
      
      nrn_2 <- reactive({as.numeric(str_split(input$arm_2_n_at_risk,",")[[1]])})
      trn_2 <- reactive({as.numeric(input$arm_2_interval) * (seq_along(nrn_2()) - 1)})
      
      IPD_2 <- reactive({getIPD(wpd_t_2(),wpd_s_2(),trn_2(),nrn_2(), arm.id = 2)})
      
      IPD_survfit_2 <- reactive({survfit(Surv(IPD_2()[,1], IPD_2()[,2])~1)})
      
      my_quantile_2 <- reactive({round(unlist(quantile(IPD_survfit_2(), probs = input$quantile_2)),2)})
      sum_IPD_2 <- reactive({summary(IPD_survfit_2(), time = input$prob_2)})
      
      
      output$plot_2 <- renderPlot({ 
        
        p2 <- survminer::ggsurvplot(IPD_survfit_2(), 
                                    data = as.data.frame(IPD_2()), 
                                    risk.table = TRUE, 
                                    legend.labs = "Arm 2",
                                    break.x.by = as.numeric(input$arm_2_interval),
                                    palette = c("blue", "red")) 
        p2$plot <- p2$plot + 
          geom_hline(yintercept = 1 - input$quantile_2) +
          geom_vline(xintercept = input$prob_2, linetype = "dashed") 
        
        p2
      })
      
      
      output$info_q_2 <- renderText({paste0(my_quantile_2()[1],
                                            " with 95% conf int (",
                                            my_quantile_2()[2], 
                                            ", ", 
                                            my_quantile_2()[3], 
                                            ")")})
      
      output$info_p_2 <- renderText({paste0(round(sum_IPD_2()$surv, 2),
                                            " with 95% conf int (",
                                            round(sum_IPD_2()$lower, 2), 
                                            ", ", 
                                            round(sum_IPD_2()$upper, 2), 
                                            ")")})
      
      
      output$download_2 <- downloadHandler(
        filename = function() {
          "IPD_2.csv"
        },
        content = function(con) {
          write.csv(IPD_2(), con, row.names = FALSE)
        }
      )
      
      IPD_both <- reactive({rbind(IPD(), IPD_2())})
      IPD_survfit_both <- reactive({survfit(Surv(IPD_both()[,1], IPD_both()[,2]) ~ IPD_both()[,3])})
      IPD_both_df <- reactive({
        df <- as.data.frame(IPD_both())
        names(df) <- c("time", "event", "arm")
        df$arm <- as.factor(df$arm)
        df
      })
      
      output$plot_both <- renderPlot({ 
        
        pboth <- survminer::ggsurvplot(IPD_survfit_both(), 
                                       data = as.data.frame(IPD_both()), 
                                       risk.table = TRUE, 
                                       legend.labs = c("Arm 1", "Arm_2"),
                                       break.x.by = as.numeric(input$arm_2_interval),
                                       palette = c("red", "blue"),
                                       conf.int = TRUE) 
        
        
        pboth
      })
      
      output$cox_text <- renderPrint({summary(coxph(Surv(time, event) ~ arm,
                                                    data = IPD_both_df()))})
      
      
      output$download_both <- downloadHandler(
        filename = function() {
          "IPD_both.csv"
        },
        content = function(con) {
          write.csv(IPD_both_df(), con, row.names = FALSE)
        }
      )
      
    })
  })
 
  
}
shinyApp(ui, server)

