# Shiny Application
library(shiny)

ui <- fluidPage(
  h1("Shiny Application")
)

server <- function(input, output) {
  # Server logic
}

shinyApp(ui = ui, server = server)

