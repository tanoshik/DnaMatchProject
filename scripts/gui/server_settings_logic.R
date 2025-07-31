register_settings_logic <- function(input, output, session, db_rv, locus_order) {
  observeEvent(input$db_upload, {
    req(input$db_upload)
    tryCatch({
      uploaded_df <- read.csv(input$db_upload$datapath, stringsAsFactors = FALSE)
      
      new_db <- prepare_database(db_file = input$db_upload$datapath,
                                 locus_file = "data/locus_order.rds")
      db_rv(new_db$prepared)  
      
      output$db_upload_status <- renderText("Database successfully loaded.")
    }, error = function(e) {
      output$db_upload_status <- renderText(paste("Error:", e$message))
    })
  })
}
