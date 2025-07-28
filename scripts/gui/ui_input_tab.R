render_input_tab <- function(left_loci, right_loci) {
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
                 column(6, fileInput("query_file", "Select Query Profile CSV")),
                 column(6,
                        checkboxInput("homo_to_any", "Convert homozygous alleles to 'any'", value = FALSE),
                        actionButton("goto_confirm", "Go to Confirm")
                 )
               )
           )
  )
}
