render_settings_tab <- function() {
  tabPanel("Settings",
           fluidRow(
             column(12,
                    h4("Database Profile Upload"),
                    fileInput("db_upload", "Upload database_profile.csv",
                              accept = c(".csv")),
                    verbatimTextOutput("db_upload_status")
             )
           )
  )
}
