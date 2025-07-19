# scripts/gui/query_gui_app.R

library(shiny)
library(rprojroot)

# Locate project root
root <- rprojroot::find_root(rprojroot::has_file("DnaMatchProject.Rproj"))

# Load data
locus_order <- readRDS(file.path(root, "data", "locus_order.rds"))
allele_count <- readRDS(file.path(root, "data", "freq_table.rds"))

# Prepare allele choices per locus
allele_choices <- lapply(locus_order, function(locus) {
  alleles <- allele_count[[locus]]
  sort(unique(c("any", names(alleles))))
})
names(allele_choices) <- locus_order

# Split loci left/right
locus_left <- locus_order[1:11]
locus_right <- locus_order[12:21]

ui <- fluidPage(
  titlePanel("Query Profile Input (GlobalFiler, Japanese)"),
  fluidRow(
    column(6, lapply(locus_left, function(locus) {
      fluidRow(
        column(4, strong(locus)),
        column(4, selectInput(paste0("a1_", locus), NULL, choices = allele_choices[[locus]], selected = "any")),
        column(4, selectInput(paste0("a2_", locus), NULL, choices = allele_choices[[locus]], selected = "any"))
      )
    })),
    column(6, lapply(locus_right, function(locus) {
      fluidRow(
        column(4, strong(locus)),
        column(4, selectInput(paste0("a1_", locus), NULL, choices = allele_choices[[locus]], selected = "any")),
        column(4, selectInput(paste0("a2_", locus), NULL, choices = allele_choices[[locus]], selected = "any"))
      )
    }))
  ),
  br(),
  fluidRow(
    column(4, fileInput("file_input", "Select Query Profile CSV", accept = ".csv")),
    column(2, radioButtons("q_homo_to_any", "Homozygous to 'any'", choices = c("No" = FALSE, "Yes" = TRUE), inline = TRUE)),
    column(2, actionButton("submit", "OK")),
    column(2, actionButton("cancel", "Cancel"))
  ),
  hr(),
  tableOutput("prepared")
)

server <- function(input, output, session) {
  query_data <- reactiveVal(NULL)
  
  observeEvent(input$file_input, {
    req(input$file_input)
    df <- read.csv(input$file_input$datapath, stringsAsFactors = FALSE)
    query_data(df)
  })
  
  observeEvent(input$cancel, {
    query_data(NULL)
    shinyjs::reset("file_input")
  })
  
  observeEvent(input$submit, {
    if (is.null(query_data())) {
      df <- data.frame(
        Locus = locus_order,
        allele1 = sapply(locus_order, function(l) input[[paste0("a1_", l)]]),
        allele2 = sapply(locus_order, function(l) input[[paste0("a2_", l)]]),
        stringsAsFactors = FALSE
      )
      query_data(df)
    }
    
    df <- query_data()
    prep <- lapply(locus_order, function(locus) {
      row <- df[df$Locus == locus, ]
      a1 <- ifelse(nrow(row) == 0 || is.na(row$allele1) || row$allele1 == "", "any", row$allele1)
      a2 <- ifelse(nrow(row) == 0 || is.na(row$allele2) || row$allele2 == "", "any", row$allele2)
      if (isTRUE(as.logical(input$q_homo_to_any)) && a1 == a2 && a1 != "any") a2 <- "any"
      c(a1, a2)
    })
    names(prep) <- locus_order
    
    output$prepared <- renderTable({
      data.frame(
        Locus = locus_order,
        allele1 = sapply(prep, `[[`, 1),
        allele2 = sapply(prep, `[[`, 2),
        stringsAsFactors = FALSE
      )
    })
  })
}

shinyApp(ui, server)
