# New function to print and save results
print_and_save_results <- function(results, output_path) {
    # Start capturing console output
    sink(file.path(output_path, "console_output.txt"))

    # Print summary statistics and test results with added line breaks
    cat("\n\n--- Summary of critical period ---\n")
    print(psych::describe(results$critical_period))
    print(Hmisc::describe(results$critical_period))

    cat("\n\n--- T-test results ---\n")
    cat("(if p < 0.05, then the mean maximum correlation is significantly different from 0)\n")
    print(results$t_test_result)

    cat("\n\n--- Correlation test results ---\n")
    cat("(if p < 0.05, then the correlation is significantly different from 0)\n")
    print(results$correlation_test)

    cat("\n\n--- GAMM model summary ---\n")
    print(summary(results$gamm_model$gam))

    cat("\n\n--- Individual seal bearing compared to particle mean bearing summary ---\n")
    print(Hmisc::describe(results$seal_following_particle %>% select(-data)))

    cat("\n\n--- Following x Survived contingency table ---\n")
    print(results$contingency_table)

    # Stop capturing console output
    sink()

    ## Save model plots
    png(file.path(output_path, "gamm_model_plots.png"), width = 6, height = 3.5, units = "in", res = 300, pointsize = 6)
    par(mfrow = c(2, 2))
    plot(results$gamm_model$gam, all.terms = TRUE, pages = 1)
    dev.off()

    # Save plots as individual files
    for (plot_name in names(results$plots)) {
        ggsave(file.path(output_path, paste0(plot_name, ".png")), results$plots[[plot_name]], width = 10, height = 8, bg = "white")
    }

    cat("Results printed and saved to", output_path, "\n")
}
