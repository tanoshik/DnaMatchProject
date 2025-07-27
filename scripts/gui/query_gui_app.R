library(shiny)

# RDSファイルを読み込む
freq_table <- readRDS("data/freq_table.rds")
locus_order <- readRDS("data/locus_order.rds")

# Amelogeninを除外し、順序を維持したまま21ローカスを取得
visible_loci <- setdiff(locus_order, "Amelogenin")

# GF仕様に基づき、左12ローカス（D8〜FGA）、右9ローカス（D22〜D2S1338）に分割
left_loci <- visible_loci[1:12]
right_loci <- visible_loci[13:21]

# アレル選択肢生成
get_allele_choices <- function(locus) {
  alleles <- freq_table$Allele[freq_table$Locus == locus]
  alleles <- alleles[!is.na(alleles) & grepl("^[0-9.]+$", alleles)]
  numeric_alleles <- as.numeric(alleles)
  sorted_alleles <- sort(unique(numeric_alleles))
  choices <- c("any", as.character(sorted_alleles))
  return(choices)
}

# 入力UI構成
renderLocusInput <- function(locus) {
  choices <- get_allele_choices(locus)
  fluidRow(
    column(2, tags$strong(locus)),
    column(5, selectInput(paste0("input_", locus, "_1"), label = NULL, choices = choices, width = "100%")),
    column(5, selectInput(paste0("input_", locus, "_2"), label = NULL, choices = choices, width = "100%"))
  )
}

# UI定義
ui <- fluidPage(
  titlePanel("Query Profile Input"),
  tabsetPanel(id = "main_tabs",
              tabPanel("Input",
                       div(style = "max-height: 100vh; overflow-y: auto;",
                           fluidRow(
                             column(6, style = "margin-top: 20px;",
                                    lapply(left_loci, renderLocusInput)
                             ),
                             column(6, style = "margin-top: 20px;",
                                    lapply(right_loci, renderLocusInput),
                                    textInput("sample_name", "Sample Name", value = "QuerySample1")
                             )
                           ),
                           hr(),
                           fluidRow(
                             column(6,
                                    fileInput("query_file", "Select Query Profile CSV")
                             ),
                             column(6,
                                    checkboxInput("homo_to_any", "Convert homozygous alleles to 'any'", value = FALSE),
                                    actionButton("goto_confirm", "Go to Confirm")
                             )
                           )
                       )
              ),
              tabPanel("Confirm",
                       fluidRow(
                         column(8,
                                h3("Prepared Query Profile"),
                                tableOutput("confirm_table")
                         ),
                         column(4,
                                div(style = "display: flex; flex-direction: column; align-items: flex-start; margin-right: 100px; margin-top: 40px;",
                                    div(
                                      actionButton("run_match", "Run Match")
                                    ),
                                    div(style = "margin-Top: 20px;",
                                          textOutput("total_freq")
                                    )
                                    
                                )
                         )
                       )
              ),
              tabPanel("Result",
                       h3("Match Results"),
                       tableOutput("result_table")
              )
  )
)

# サーバー定義
server <- function(input, output, session) {
  source("scripts/utils_profile.R")
  source("scripts/io_profiles.R")
  source("scripts/scoring.R")
  source("scripts/matcher.R")
  source("scripts/utils_freq.R")
  
  query_profile_reactive <- reactiveVal(NULL)
  match_result_reactive <- reactiveVal(NULL)
  freq_table_df_reactive <- reactiveVal(NULL)
  total_freq_reactive <- reactiveVal(NULL)
  
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
  
  observeEvent(input$run_match, {
    profile_df <- query_profile_reactive()
    profile <- split(profile_df[, c("Allele1", "Allele2")], profile_df$Locus)
    
    db <- read_db_profiles("data/database_profile.csv", locus_order, homo_to_any = FALSE)
    result <- run_match(profile, db, top_n = 10)
    result$score_df$Score <- as.integer(result$score_df$Score)
    match_result_reactive(result$score_df)
    
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
    result
  })
  
  output$total_freq <- renderText({
    total <- total_freq_reactive()
    if (is.null(total)) return("")
    paste("Total Frequency:", format(total, scientific = TRUE, digits = 2))
  })
}

app <- shinyApp(ui = ui, server = server)
if (interactive()) runApp(app)
app
