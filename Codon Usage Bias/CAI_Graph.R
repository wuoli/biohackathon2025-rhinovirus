
rhinovirus_values <- c(0.582252, 0.577522, 0.605716, 0.588947, 0.570363, 0.588305, 
                       0.620437, 0.572831, 0.592268, 0.603954, 0.585898, 0.582957,
                       0.585307, 0.599864, 0.576769, 0.570005)
west_nile_values <- c(0.626476, 0.625656, 0.626079, 0.625656, 0.625915, 0.625771, 
                      0.626051, 0.626657, 0.627439, 0.627065, 0.625563, 0.627509)


data <- data.frame(
  virus = c(rep("Rhinovirus", length(rhinovirus_values)), rep("West Nile", length(west_nile_values))),
  cai_values = c(rhinovirus_values, west_nile_values)
)

summary_data <- data %>%
  group_by(virus) %>%
  summarise(
    mean_cai = mean(cai_values),
    sd_cai = sd(cai_values),
    .groups = "drop"
  )


ggplot(summary_data, aes(x = virus, y = mean_cai, fill = virus)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = mean_cai - sd_cai, ymax = mean_cai + sd_cai), 
                width = 0.2) +
  scale_fill_manual(values = c(
    "Rhinovirus" = "#66c2a5",
    "West Nile" = "#8da0cb"
  )) +
  labs(title = "Mean CAI by Virus",
       y = "Mean CAI",
       x = "",
       fill = "Virus") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

wilcox_test_result <- wilcox.test(cai_values ~ virus, data = data, exact = FALSE)
print(wilcox_test_result)

print(summary_data)
