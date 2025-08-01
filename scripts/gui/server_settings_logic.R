register_settings_logic <- function(input, output, session, db_rv, locus_order,
                                    db_file_name_rv, db_count_rv) {
  
  observeEvent(input$db_upload, {
    req(input$db_upload)
    tryCatch({
      uploaded_df <- read.csv(input$db_upload$datapath, stringsAsFactors = FALSE)
      new_db <- prepare_database_df(uploaded_df, locus_order)
      db_rv(new_db$prepared)
      
      # ← 渡された reactiveVal に反映
      db_file_name_rv(input$db_upload$name)
      db_count_rv(length(new_db$prepared))
      
      output$db_upload_status <- renderText("Database successfully loaded.")
    }, error = function(e) {
      output$db_upload_status <- renderText(paste("Error:", e$message))
    })
  })
  
  output$db_file_display <- renderText({
    fname <- db_file_name_rv()
    if (!is.null(fname)) paste("Current DB File:", fname) else ""
  })
  
  output$db_count_display <- renderText({
    count <- db_count_rv()
    if (!is.null(count)) paste("Sample Count:", count) else ""
  })
}
