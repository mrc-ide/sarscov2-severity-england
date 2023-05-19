version_check("sircovid", "0.14.11")
version_check("spimalot", "0.8.23")

dir.create("outputs", FALSE, TRUE)

# Read scenario inputs
national <- load_national_estimates("input/")

# Plot sensitivity analyses
png("outputs/national_ve_scenarios.png", units = "in", width = 14, height = 10, res = 300)
plot_national(national, which = "ve_scenarios")
dev.off()

png("outputs/national_other_scenarios.png", units = "in", width = 14, height = 10, res = 300)
plot_national(national, which = "other_scenarios")
dev.off()

png("outputs/national_crim_scenarios.png", units = "in", width = 14, height = 10, res = 300)
plot_national(national, which = "crim_scenarios")
dev.off()
