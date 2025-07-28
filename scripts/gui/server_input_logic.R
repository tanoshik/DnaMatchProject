register_input_logic <- function(input, output, session, visible_loci, freq_table,
                                 query_profile_reactive, freq_table_df_reactive, total_freq_reactive) {
  observeEvent(input$query_file, {
    req(input$query_file)
    
    for (locus in visible_loci) {
      updateSelectInput(session, paste0("input_", locus, "_1"), selected = "any")
      updateSelectInput(session, paste0("input_", locus, "_2"), selected = "any")
    }
    
    df <- read.csv(input$query_file$datapath, stringsAsFactors = FALSE)
    names(df) <- tolower(names(df))
    
    required_cols <- c("locus", "allele1", "allele2")
    if (!all(required_cols %in% names(df))) {
      showModal(modalDialog(
        title = "Invalid CSV format",
        "CSV file must contain columns: Locus, Allele1, Allele2",
        easyClose = TRUE
      ))
      return()
    }
    
    if ("sampleid" %in% names(df)) {
      updateTextInput(session, "sample_name", value = df$sampleid[1])
      df <- df[df$sampleid == df$sampleid[1], ]
    }
    
    df <- df[, c("locus", "allele1", "allele2")]
    colnames(df) <- c("Locus", "Allele1", "Allele2")
    
    df$Allele1[is.na(df$Allele1) | df$Allele1 == ""] <- "any"
    df$Allele2[is.na(df$Allele2) | df$Allele2 == ""] <- "any"
    
    df <- df[df$Locus %in% visible_loci, ]
    missing_loci <- setdiff(visible_loci, df$Locus)
    if (length(missing_loci) > 0) {
      df <- rbind(df, data.frame(
        Locus = missing_loci,
        Allele1 = "any",
        Allele2 = "any"
      ))
    }
    
    df$Locus <- factor(df$Locus, levels = visible_loci)
    df <- df[order(df$Locus), ]
    
    for (locus in df$Locus) {
      a1 <- df[df$Locus == locus, "Allele1"]
      a2 <- df[df$Locus == locus, "Allele2"]
      valid_choices <- get_allele_choices(locus)
      
      if (!(a1 %in% valid_choices)) {
        showNotification(paste("Allele", a1, "is not valid for", locus, "- set to blank."))
        a1 <- ""
      }
      if (!(a2 %in% valid_choices)) {
        showNotification(paste("Allele", a2, "is not valid for", locus, "- set to blank."))
        a2 <- ""
      }
      
      updateSelectInput(session, paste0("input_", locus, "_1"), selected = a1)
      updateSelectInput(session, paste0("input_", locus, "_2"), selected = a2)
    }
  })
  
  observeEvent(input$goto_confirm, {
    updateTabsetPanel(session, "main_tabs", selected = "Confirm")
    
    query_df <- data.frame(
      Locus = visible_loci,
      Allele1 = sapply(visible_loci, function(locus) {
        val <- input[[paste0("input_", locus, "_1")]]
        if (is.null(val) || val == "") "any" else val
      }),
      Allele2 = sapply(visible_loci, function(locus) {
        val <- input[[paste0("input_", locus, "_2")]]
        if (is.null(val) || val == "") "any" else val
      }),
      stringsAsFactors = FALSE
    )
    
    prepared <- prepare_profile_df(query_df, homo_to_any = input$homo_to_any)
    query_profile_reactive(prepared)
    
    freq_df <- calc_freq_loci_df(prepared, freq_table)
    total_freq <- calc_total_freq(prepared, freq_table)
    freq_table_df_reactive(freq_df)
    total_freq_reactive(total_freq)
  })
}
