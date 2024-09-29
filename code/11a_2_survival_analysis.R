# Function to generate plots and perform tests based on survival
plot_survival_analysis <- function(indiv_summ, masterData, loc_early, locf_early, pt_early) {
    # Purpose: This function analyzes and visualizes the relationship between dispersal direction and survival rates for seal pups

    # Step 1: Merge survival data with individual summaries
    d <- readRDS("./Output/all_data_combined.rds")
    ds <- d %>%
        group_by(id) %>%
        filter(trip == 1, sim == 0) %>%
        summarise(
            birthyear = first(birthyear),
            surviveTrip1 = ifelse(first(is_trip_complete) == TRUE | first(seen_6m) == TRUE, TRUE, FALSE),
            surviveYear1 = ifelse(first(seen_1y) == TRUE, TRUE, FALSE),
            weanmass = first(weanmass)
        )



    b <- indiv_summ %>%
        ungroup() %>%
        left_join(ds)

    # Step 2: Categorize seals by weight
    b <- b %>% mutate(
        sizeCat = case_when(
            weanmass > 135 ~ "heavy",
            weanmass < 96 ~ "light",
            TRUE ~ "avg"
        )
    )

    # Step 3: Analyze survival for the first trip
    b <- b %>% mutate(survive = surviveTrip1)

    # Create circular data for seals that survived and didn't survive the first trip
    circ_live <- circular(b %>% filter(survive == TRUE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
    circ_dead <- circular(b %>% filter(survive == FALSE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")

    # Plot histograms for 1st trip survival
    p1 <- plot_circular_histogram(circ_live, "Survived Trip 1", paste0("1st trip survived, n = ", length(circ_live)), paste0(output_path, "survived_trip1_histogram.png"))
    p2 <- plot_circular_histogram(circ_dead, "Died Trip 1", paste0("1st trip died, n = ", length(circ_dead)), paste0(output_path, "died_trip1_histogram.png"))

    # Step 4: Analyze survival for the first year
    b <- b %>% mutate(survive = surviveYear1)

    # Create circular data for seals that survived and didn't survive the first year
    circ_live_y1 <- circular(b %>% filter(survive == TRUE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")
    circ_dead_y1 <- circular(b %>% filter(survive == FALSE) %>% pull(bearing) %% 360, units = "degrees", zero = pi / 2, rotation = "clock")

    # Plot histograms for 1st year survival
    p3 <- plot_circular_histogram(circ_live_y1, "Survived Year 1", paste0("1st year survived, n = ", length(circ_live_y1)), paste0(output_path, "survived_year1_histogram.png"))
    p4 <- plot_circular_histogram(circ_dead_y1, "Died Year 1", paste0("1st year died, n = ", length(circ_dead_y1)), paste0(output_path, "died_year1_histogram.png"))

    # Step 5: Create and save a combined plot of all survival analyses
    combined_plot <- ggarrange(
        p1,
        p2,
        p3,
        p4,
        ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), align = "hv"
    )
    ggsave(paste0(output_path, "dispersal_direction_by_survival.png"), plot = combined_plot, width = 20, height = 20, units = "cm", dpi = 500)

    # Interpretation: These circular histograms show the distribution of dispersal directions for different groups.
    # Compare the patterns between survivors and non-survivors to see if there are any preferred directions associated with survival.
    # Look for peaks in the histograms that might indicate favorable directions for survival.

    # Step 6: Perform statistical tests to compare dispersal directions
    # Return the results of the statistical tests
    list(
        watson_live_dead_trip1 = watson.two.test(circ_live, circ_dead), watson_p_live_dead_trip1 = watson.two.test(circ_p, circ_dead),
        watson_p_female = watson.two.test(circ_p, circ_f),
        HR_live = HR_test(circ_live),
        HR_dead = HR_test(circ_dead)
    )
}
