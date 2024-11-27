library(knitr)
library(here)

# Create directories if they don't exist
dir.create(here("publication", "code"), recursive = TRUE, showWarnings = FALSE)

# Function to add frontmatter to markdown files
add_frontmatter <- function(file_path, title) {
    content <- readLines(file_path)
    frontmatter <- c(
        "---",
        paste0("title: ", title),
        "layout: default",
        "---",
        "",
        content
    )
    writeLines(frontmatter, file_path)
}

# Convert dispersal-analysis.qmd
knitr::knit(
    input = here("code", "dispersal-analysis.qmd"),
    output = here("publication", "code", "dispersal-analysis.md")
)

# Convert export-data.qmd
knitr::knit(
    input = here("code", "export-data.qmd"),
    output = here("publication", "code", "export-data.md")
)

# Convert functions.R to markdown
functions_content <- readLines(here("code", "functions", "functions.R"))
functions_md <- c(
    "# Utility Functions",
    "",
    "```r",
    functions_content,
    "```"
)
writeLines(functions_md, here("publication", "code", "functions.md"))

# Add frontmatter to each file
add_frontmatter(
    here("publication", "code", "dispersal-analysis.md"),
    "Dispersal Analysis"
)
add_frontmatter(
    here("publication", "code", "export-data.md"),
    "Data Export"
)
add_frontmatter(
    here("publication", "code", "functions.md"),
    "Utility Functions"
)