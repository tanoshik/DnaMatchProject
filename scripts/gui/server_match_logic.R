register_match_logic <- function(input, output, session, db, visible_loci,
                                 query_profile_reactive, match_result_reactive, db_count) {
  
  output$db_count <- renderText({
    paste("Database Samples : N =", db_count)
  })
  
  observeEvent(input$run_match, {
    profile_df <- query_profile_reactive()
    profile <- split(profile_df[, c("Allele1", "Allele2")], profile_df$Locus)
    
    result <- run_match(profile, db)
    
    score_df <- result$score_df
    score_df$Score <- as.integer(score_df$Score)
    
    filtered_df <- switch(input$filter_type,
                          "top_n" = {
                            n <- input$top_n
                            head(score_df[order(-score_df$Score), ], n)
                          },
                          "score_min" = {
                            min_score <- input$min_score
                            score_df[score_df$Score >= min_score, ]
                          },
                          "all" = {
                            score_df
                          }
    )
    
    match_result_reactive(filtered_df)
    updateTabsetPanel(session, "main_tabs", selected = "Result")
  })
  
  output$confirm_table <- renderTable({
    df <- query_profile_reactive()
    if (is.null(df)) return(NULL)
    
    freq_df <- freq_table_df_reactive()
    df$Freq <- sprintf("%.6f", freq_df$Freq)
    df
  })
  
  output$result_table <- renderTable({
    result <- match_result_reactive()
    if (is.null(result)) return(NULL)
    result[order(-result$Score), ]
  })
  
  output$total_freq <- renderText({
    total <- total_freq_reactive()
    if (is.null(total)) return("")
    paste("Total Frequency:", format(total, scientific = TRUE, digits = 2))
  })
  
  output$download_result <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0("match_results_", timestamp, ".csv")
    },
    content = function(file) {
      result <- match_result_reactive()
      if (is.null(result)) return(NULL)
      write.csv(result, file, row.names = FALSE)
    }
  )
}
