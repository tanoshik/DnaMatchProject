render_result_tab <- function() {
  tabPanel("Result",
           fluidRow(
             column(8,
                    h3("Match Results"),
                    div(style = "display: flex; justify-content: flex-end; gap: 10px;",
                        downloadButton("download_result", label = "Summary"),
                        downloadButton("download_detail", label = "Detail")
                    )
             ),
           ),
           tableOutput("result_table")
  )
}
