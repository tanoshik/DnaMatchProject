library(shiny)

# freq_table.rds を読み込む（前提：data/ に存在）
freq_table <- readRDS("data/freq_table.rds")

get_allele_choices <- function(locus) {
  alleles <- freq_table$Allele[freq_table$Locus == locus]
  alleles <- alleles[!is.na(alleles) & grepl("^[0-9.]+$", alleles)]
  numeric_alleles <- as.numeric(alleles)
  sorted_alleles <- sort(unique(numeric_alleles))
  choices <- c("any", as.character(sorted_alleles))
  return(choices)
}

# ローカス定義（Amelogenin除外）
left_loci <- c("D3S1358", "TH01", "D21S11", "D18S51", "Penta_E", "D5S818", "D13S317", "D7S820", "D16S539", "CSF1PO")
right_loci <- c("vWA", "D8S1179", "TPOX", "Penta_D", "D19S433", "FGA", "D2S1338", "D6S1043", "D12S391", "SE33")

# 入力UI構成（labelとselectInputの縦位置を揃える）
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
                       div(style = "max-height: 100vh; overflow-y: auto;",  # スクロールバー防止・高さ調整
                           fluidRow(
                             column(6,
                                    h4("Left Loci"),
                                    lapply(left_loci, renderLocusInput)
                             ),
                             column(6,
                                    h4("Right Loci"),
                                    lapply(right_loci, renderLocusInput)
                             )
                           ),
                           hr(),
                           fileInput("query_file", "Select Query Profile CSV"),
                           actionButton("goto_confirm", "Go to Confirm")
                       )
              ),
              tabPanel("Confirm",
                       h3("Prepared Query Profile"),
                       tableOutput("confirm_table")
              )
  )
)

# サーバー定義
server <- function(input, output, session) {
  
  # CSV選択後の読み込み＆UI更新＆Confirm遷移
  observeEvent(input$query_file, {
    req(input$query_file)
    query_df <- read.csv(input$query_file$datapath, stringsAsFactors = FALSE)
    
    for (locus in c(left_loci, right_loci)) {
      a1 <- query_df[query_df$Locus == locus, "Allele1"]
      a2 <- query_df[query_df$Locus == locus, "Allele2"]
      updateSelectInput(session, paste0("input_", locus, "_1"), selected = a1)
      updateSelectInput(session, paste0("input_", locus, "_2"), selected = a2)
    }
    
    updateTabsetPanel(session, "main_tabs", selected = "Confirm")
  })
  
  # Go to Confirmボタン（手動遷移）
  observeEvent(input$goto_confirm, {
    updateTabsetPanel(session, "main_tabs", selected = "Confirm")
  })
  
  # 確認用表示（Confirmタブ）
  output$confirm_table <- renderTable({
    data.frame(
      Locus = c(left_loci, right_loci),
      Allele1 = sapply(c(left_loci, right_loci), function(locus) input[[paste0("input_", locus, "_1")]]),
      Allele2 = sapply(c(left_loci, right_loci), function(locus) input[[paste0("input_", locus, "_2")]])
    )
  })
  
}

shinyApp(ui, server)
