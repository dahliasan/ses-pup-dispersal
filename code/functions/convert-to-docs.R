library(here)

# Create code directory if it doesn't exist
dir.create(here("code"), recursive = TRUE, showWarnings = FALSE)

# Convert functions.R to qmd
functions_content <- readLines(here("code", "functions", "functions.R"))
functions_qmd <- c(
    "---",
    "title: Utility Functions",
    "format:",
    "  html:",
    "    code-fold: show",
    "    code-tools: true",
    "---",
    "",
    "## R Functions",
    "",
    "```{r}",
    "#| echo: true",
    "#| output: false",
    functions_content,
    "```"
)
writeLines(functions_qmd, here("code", "functions.qmd"))

message("Created functions.qmd in code directory")
