render_settings_tab <- function() {
  tabPanel("Settings",
           fluidRow(
             column(12,
                    h4("Database Profile Upload"),
                    fileInput("db_upload", "Upload database_profile.csv", accept = ".csv"),
                    verbatimTextOutput("db_upload_status"),
                    br(),
                    textOutput("db_file_display"),    # ← 現在のDBファイル名
                    textOutput("db_count_display")    # ← サンプル数
             )
           )
  )
}
