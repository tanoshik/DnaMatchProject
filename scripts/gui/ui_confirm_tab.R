render_confirm_tab <- function(visible_loci) {
  tabPanel("Confirm",
           fluidRow(
             column(8,
                    h3("Prepared Query Profile"),
                    tableOutput("confirm_table")
             ),
             column(4,
                    div(style = "display: flex; flex-direction: column; align-items: flex-start; margin-left: 40px; min-width: 200px; max-width: 260px; margin-top: 40px;",
                        div(actionButton("run_match", "Run Match")),
                        div(style = "margin-top: 30px;",
                            radioButtons("filter_type", "Display Limit",
                                         choices = c("Top N" = "top_n", "Score ≥ n" = "score_min", "All" = "all"),
                                         selected = "top_n"
                            ),
                            conditionalPanel("input.filter_type == 'top_n'",
                                             numericInput("top_n", "Top N", value = 10, min = 1, width = "150px")
                            ),
                            conditionalPanel("input.filter_type == 'score_min'",
                                             numericInput("min_score", "Minimum Score", value = length(visible_loci) * 2, min = 1, width = "150px")
                            )
                        ),
                        div(style = "margin-Top: 20px;",
                            textOutput("total_freq"),
                            textOutput("db_count")
                        )
                    )
             )
           )
  )
}
