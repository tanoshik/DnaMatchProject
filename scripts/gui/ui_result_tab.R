render_result_tab <- function() {
  tabPanel("Result",
           h3("Match Results"),
           downloadButton("download_result", "Download CSV"),
           tableOutput("result_table")
  )
}
