library(knitr)
library(here)
library(stringr)

# Create directories if they don't exist
dir.create(here("docs", "code"), recursive = TRUE, showWarnings = FALSE)

# Function to process qmd content and respect chunk options
process_qmd <- function(input_file) {
    content <- readLines(input_file)
    in_chunk <- FALSE
    output_lines <- character()
    skip_chunk <- FALSE

    for (line in content) {
        if (str_detect(line, "^```\\{.*\\}")) {
            # Start of code chunk
            in_chunk <- TRUE
            chunk_options <- str_extract(line, "\\{.*\\}")
            skip_chunk <- str_detect(chunk_options, "output:\\s*false")
            output_lines <- c(output_lines, line)
        } else if (str_detect(line, "^```$") && in_chunk) {
            # End of code chunk
            in_chunk <- FALSE
            skip_chunk <- FALSE
            output_lines <- c(output_lines, line)
        } else if (!in_chunk || !skip_chunk) {
            # Regular line or code chunk that should be included
            output_lines <- c(output_lines, line)
        }
    }

    return(output_lines)
}

# Function to add frontmatter
add_frontmatter <- function(file_path, title) {
    content <- readLines(file_path)
    frontmatter <- c(
        "---",
        paste0("title: ", title),
        "layout: default",
        "nav_order: 2",
        "---",
        "",
        content
    )
    writeLines(frontmatter, file_path)
}

# Convert dispersal-analysis.qmd
dispersal_content <- process_qmd(here("code", "dispersal-analysis.qmd"))
temp_file <- tempfile(fileext = ".qmd")
writeLines(dispersal_content, temp_file)
knitr::knit(
    input = temp_file,
    output = here("docs", "code", "dispersal-analysis.md")
)

# Convert export-data.qmd
export_content <- process_qmd(here("code", "export-data.qmd"))
temp_file <- tempfile(fileext = ".qmd")
writeLines(export_content, temp_file)
knitr::knit(
    input = temp_file,
    output = here("docs", "code", "export-data.md")
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
writeLines(functions_md, here("docs", "code", "functions.md"))

# Add frontmatter to each file
add_frontmatter(
    here("docs", "code", "dispersal-analysis.md"),
    "Dispersal Analysis"
)
add_frontmatter(
    here("docs", "code", "export-data.md"),
    "Data Export"
)
add_frontmatter(
    here("docs", "code", "functions.md"),
    "Utility Functions"
)
