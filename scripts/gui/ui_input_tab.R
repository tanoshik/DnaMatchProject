# --- get_allele_choices() の定義 ---
# freq_tableから指定ローカスのアレル一覧を取得し、"any"を先頭に配置
get_allele_choices <- function(locus, freq_table) {
  alleles <- freq_table[freq_table$Locus == locus, "Allele"]
  alleles <- unique(as.character(alleles))
  alleles <- sort(alleles[alleles != "any"])  # "any" を除いて昇順
  return(c("any", alleles))  # "any" を先頭に追加
}

# --- renderLocusInput(): ローカス単位の入力UI生成 ---
renderLocusInput <- function(locus, freq_table) {
  choices <- get_allele_choices(locus, freq_table)
  fluidRow(
    column(2, tags$div(style = "margin-top: 8px;", locus)),
    column(5, selectInput(paste0("input_", locus, "_1"), label = NULL, choices = choices, selected = "any")),
    column(5, selectInput(paste0("input_", locus, "_2"), label = NULL, choices = choices, selected = "any"))
  )
}

# --- render_input_tab(): Inputタブ全体のUI構成 ---
render_input_tab <- function(left_loci, right_loci, freq_table) {
  tabPanel("Input",
           div(style = "max-height: 100vh; overflow-y: auto;",
               fluidRow(
                 column(6, style = "margin-top: 20px;",
                        lapply(left_loci, function(locus) renderLocusInput(locus, freq_table))
                 ),
                 column(6, style = "margin-top: 20px;",
                        lapply(right_loci, function(locus) renderLocusInput(locus, freq_table)),
                        textInput("sample_name", "Sample Name", value = "QuerySample1")
                 )
               ),
               hr(),
               fluidRow(
                 column(6, fileInput("query_file", "Select Query Profile CSV")),
                 column(6,
                        checkboxInput("homo_to_any", "Convert homozygous alleles to 'any'", value = FALSE),
                        actionButton("goto_confirm", "Go to Confirm")
                 )
               )
           )
  )
}
