# Calculate correlation matrix
cor_matrix <- cor(modal_data[, predictor_vars], use = "pairwise.complete.obs")

# Visualize correlation matrix
png(paste0(output_path, "/correlation_matrix.png"), width = 800, height = 800)
corrplot(cor_matrix,
    method = "color", type = "upper", order = "hclust",
    tl.col = "black", tl.srt = 45
)
dev.off()

# Function to find highly correlated variable pairs
find_high_cor <- function(cor_matrix, threshold = 0.7) {
    cor_matrix %>%
        as.data.frame() %>%
        rownames_to_column("var1") %>%
        pivot_longer(-var1, names_to = "var2", values_to = "correlation") %>%
        filter(abs(correlation) > threshold, correlation < 1, var1 < var2) %>%
        mutate(correlation = round(correlation, 2)) %>%
        arrange(desc(abs(correlation)))
}

# Find highly correlated pairs
high_cor_pairs <- find_high_cor(cor_matrix)

cat("\n\nHighly correlated predictor pairs:\n")
print(high_cor_pairs)

# Function to remove one variable from each highly correlated pair
get_variables_to_remove <- function(cor_pairs, cor_matrix) {
    vars_to_remove <- cor_pairs %>%
        rowwise() %>%
        mutate(
            cor_var1 = mean(abs(cor_matrix[var1, ])),
            cor_var2 = mean(abs(cor_matrix[var2, ])),
            to_remove = ifelse(cor_var1 > cor_var2, var1, var2)
        ) %>%
        pull(to_remove) %>%
        unique()

    return(vars_to_remove)
}

cat("\n\n--- Variables to remove ---\n")
print(get_variables_to_remove(high_cor_pairs, cor_matrix))
