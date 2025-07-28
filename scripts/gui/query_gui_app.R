source("scripts/utils_profile.R")
source("scripts/io_profiles.R")
source("scripts/scoring.R")
source("scripts/matcher.R")
source("scripts/utils_freq.R")

source("scripts/gui/ui_input_tab.R")
source("scripts/gui/ui_confirm_tab.R")
source("scripts/gui/ui_result_tab.R")
source("scripts/gui/server_input_logic.R")
source("scripts/gui/server_match_logic.R")

freq_table <- readRDS("data/freq_table.rds")
locus_order <- readRDS("data/locus_order.rds")
visible_loci <- setdiff(locus_order, "Amelogenin")
left_loci <- visible_loci[1:12]
right_loci <- visible_loci[13:21]

db <- read_db_profiles("data/database_profile.csv", locus_order)
db_count <- length(db)

query_profile_reactive <- reactiveVal(NULL)
match_result_reactive <- reactiveVal(NULL)
freq_table_df_reactive <- reactiveVal(NULL)
total_freq_reactive <- reactiveVal(NULL)

ui <- fluidPage(
  titlePanel("Query Profile Input"),
  tabsetPanel(id = "main_tabs",
              render_input_tab(left_loci, right_loci),
              render_confirm_tab(visible_loci),
              render_result_tab()
  )
)

server <- function(input, output, session) {
  register_input_logic(input, output, session, visible_loci, freq_table,
                       query_profile_reactive, freq_table_df_reactive, total_freq_reactive)
  register_match_logic(input, output, session, db, visible_loci,
                       query_profile_reactive, match_result_reactive, db_count)
}

app <- shinyApp(ui = ui, server = server)
if (interactive()) runApp(app)
app
