
all_variant_colours <- c("#8a8f9c", "#1ea1eb", "#f76002", "#eb2f97")


plot_national <- function (national, period = "Emergence3",
                           which = "ve_scenarios") {
  
  if (which == "ve_scenarios") {
    which <- c("central", "alpha_ve_high", "alpha_ve_low",
               "delta_ve_high", "delta_ve_low",
               "booster_ve_high", "booster_ve_low")
  } else if (which == "crim_scenarios") {
    which <- c("central", "crim_infect_high", "crim_infect_low",
               "crim_hospi_high", "crim_hospi_low",
               "crim_death_high", "crim_death_low")
  } else if (which == "other_scenarios") {
    which <- c("central", "fixed_si_high", "fixed_si_low",
               "mu_d_summer", "mu_d_winter")
  } else {
    stop("Scenarios not supported")
  }
  
  national <- national[which]
  national <- switch_levels(national)[[period]]
  scenarios <- names(national)
  
  re_arrange <- function(x, w) {
    
    tmp <- switch_levels(x)[[w]]
    data.frame(matrix(unlist(tmp), nrow = length(tmp), byrow = TRUE)) %>%
      `colnames<-`(., c("mean", "lb", "ub")) %>%
      mutate(variant = factor(names(x), levels = names(x)))
  }
  
  
  ret <- NULL
  for (w in names(national[[1]][[1]])) {
    tmp <- purrr::map(names(national), function(x) re_arrange(national[[x]], w))
    names(tmp) <- scenarios
    
    tmp <- bind_rows(tmp, .id = "scenario") %>%
      pivot_longer(!c(scenario, variant)) %>%
      pivot_wider() %>%
      mutate(source = w)
    
    ret <- rbind(ret, tmp)
  }
  
  ret <- ret %>%
    mutate(source = factor(source,
                           levels = c("r0", "IHR", "HFR", "IFR"),
                           labels = c("Intrinsic R_0", "Basic IHR", "Baisc HFR", "Basic IFR"))) %>%
    mutate(scenario = factor(scenario))
  
  ret$scenario <- relevel(ret$scenario, ref = "central")

  
  ggplot(ret, aes(x = NA, y = mean, col = variant, fill = variant)) +
    geom_point(shape = 18, size = 2, alpha = 0.7, position = position_dodge(width = 0.5)) +
    geom_crossbar(aes(ymin = lb, ymax = ub), position = "dodge", width = 0.5, alpha = 0.15) +
    scale_color_manual(values = all_variant_colours) +
    scale_fill_manual(values = all_variant_colours) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, NA)) +
    scale_x_discrete(breaks = NULL) +
    labs(x = "", y = "") +
    facet_grid(cols = vars(scenario), rows = vars(source), scales = "free_y",
               switch = "y") +
    theme_bw() +
    theme(axis.line = element_line(),
          axis.title.y = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.position = "bottom",
          axis.text.x = element_blank(),
          strip.text.x = element_text(size = 10))

}


load_national_estimates <- function(dir) {
  
  scenarios <- grep("national_", list.files(dir), value = TRUE)
  
  pars_to_get <- function(file) {
    readRDS(file = paste0(dir, file))
  }
  
  out <- purrr::map(scenarios, function(x)  pars_to_get(x))
  
  nms <- scenarios %>%
    str_replace(., "national_sev_", "") %>%
    str_replace(., ".rds", "")
  names(out) <- nms
  out
  
}
