# アレル選択肢を取得：freq_table から該当ローカスのアレルを昇順で抽出し、"any" を末尾に追加
get_allele_choices <- function(locus, freq_table) {
  alleles <- freq_table$Allele[freq_table$Locus == locus]
  alleles <- unique(as.character(alleles))
  alleles <- alleles[alleles != "any"]
  if (length(alleles) == 0) return("any")
  alleles <- sort(as.numeric(alleles), na.last = NA)
  choices <- c(as.character(alleles), "any")
  return(choices)
}

# ローカスごとのアレル入力UI生成（1つのselectInputのみ表示）
renderLocusInput <- function(locus, freq_table) {
  tagList(
    selectInput(
      inputId = paste0("input_", locus, "_1"),
      label = paste(locus, "Allele 1"),
      choices = get_allele_choices(locus, freq_table),
      selected = "any"
    ),
    selectInput(
      inputId = paste0("input_", locus, "_2"),
      label = paste(locus, "Allele 2"),
      choices = get_allele_choices(locus, freq_table),
      selected = "any"
    )
  )
}

# Inputタブの全体UI構築
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
