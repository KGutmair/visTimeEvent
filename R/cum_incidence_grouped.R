###############################
# Function for cumulative incidence plots of two groups with competing events
# Plotted will be only the first cometing event
##################################

#' Function for cumulative incidence plots of two groups with competing events
#'
#' @description
#' This function plots the cumualative incidence curves of two groups with competing
#' events. Plotted will only be the first competing event (first level). One can optinally display the
#' median incidence time or incidence after a certain time.
#'
#' @inheritParams km_grouped
#' @param time A character string specifying the column name of the numeric time
#'             variable for the time-to-event endpoint.
#' @param event A character string specifying the column name of the factor event
#'              variable for the time-to-event endpoint.
#' @param time_vec_prob numeric vector containing time points, on which one wants to display
#'                      the probability of the incidence.
#' @param timepoint_label character describing the timepoint + unit for which the survival probabilities are calculated
#' @param change_factor number for changing one time unit in another unit e.g. months in time
#' @param risk_table_element options: "n.risk", "cum.censor", "cum.event", which elements should be displayed in the risk tables
#' @param x_unit unit of the time (e.g. years)
#'
#' @return a list with the cumulative incidence plot, a table with the median incidence
#' and a table with the probability of an event after certain timepoints
#' @export
#'
#' @importFrom dplyr rename %>%
#' @import tidycmprsk
#' @importFrom colorRamps primary.colors
#' @import ggplot2
#' @importFrom ggsurvfit ggcuminc
#' @importFrom stats as.formula
#' @import grid

comp_risk_grouped <- function(data,
                              time,
                              event,
                              group,
                              title = "",
                              x_title = waiver(),
                              y_title = waiver(),
                              main_outcome,
                              time_survival = 1,
                              timepoint_label = "3-yr Prob.",
                              x_lim = NULL,
                              y_lim = c(0, 1),
                              x_breaks = waiver(),
                              y_breaks = seq(0, 1, by = 0.2),
                              colors = FALSE,
                              # options: median, probability
                              show_label = "none",
                              text_size_title = 15,
                              text_size = 12,
                              p_placement = c(0.05, 0.15),
                              legend_placement = c(0.6, 0.8),
                              time_vec_prob = c(1, 2, 3),
                              change_xscale = TRUE,
                              change_factor = 12,
                              risk_table_elements = c("n.risk", "cum.censor", "cum.event"),
                              x_unit = "years") {
  #----------------------------------------------------------
  # Checking input variables
  #----------------------------------------------------------

  assert_numeric(data[[time]], lower = 0, any.missing = FALSE)
  assert_factor(data[[event]], any.missing = FALSE)
  assert_factor(data[[group]], any.missing = FALSE, max.levels = 2, min.levels = 2)

  #------------------------------------------------------
  # Creation of survival object
  #-------------------------------------------------------

  formula <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", group))

  surv_object <- tidycmprsk::cuminc(
    formula = formula,
    data = data, conf.type = "arcsin"
  )

  # specifying colors, if they were not specified in the parameters
  if (is.logical(colors)) {
    n <- length(unique(data[[group]]))
    colors <- colorRamps::primary.colors(n + 1)[-1]
  }
  cr <- levels(data[[event]])[-1]
  surv_object$tidy$outcome <- factor(surv_object$tidy$outcome, levels = cr, labels = cr)


  #-----------------------------------------------------
  # Survival probability after a certain timepoint
  #----------------------------------------------------
  labels <- paste0("{time} ", x_unit)

  tbl <-
    surv_object %>%
    tbl_cuminc(
      times = time_survival, label_header = labels,
      outcomes = cr
    )

  table_surv_prob <- tbl$tidy

  # b <- tbl$tidy
  # # restructure the probability output
  # b$statistic <- gsub("\\((\\d+)%?,\\s*(\\d+)%?\\)", "(\\1-\\2)", b$statistic)

  #---------------------------------------------------
  # Median survival probability
  #-------------------------------------------------
  median_table <- surv_object$tidy %>%
    filter(.data$estimate >= 0.5) %>%
    group_by(.data$outcome) %>%
    arrange(.data$outcome, .data$time) %>%
    slice_head() %>%
    select(.data$time)


  #--------------------------------------------------------------
  # Output table for survival probabilities at multiple timepoints
  #-----------------------------------------------------------

  tab <- surv_object %>%
    tbl_cuminc(
      times = time_vec_prob,
      label_header = labels,
      outcomes = cr,
      estimate_fun = function(x) style_number(x, scale = 100, digits = 0)
    ) %>%
    add_nevent() %>%
    add_n() %>%
    add_nevent(location = c("label", "level")) %>%
    add_n(location = c("label", "level"))

  tab <- as.data.frame(tab)

  if ("**Group**" %in% names(tab)) {
    tab <- tab %>%
      select(-`**Characteristic**`) %>%
      rename(
        `Competing events` = `**Group**`,
        N = `**N**`,
        `N Event` = `**N Event**`
      )
    cols <- 4:ncol(tab)
    tab[cols] <- lapply(tab[cols], function(x) gsub("\\((\\d+)%?,\\s*(\\d+)%?\\)", "(\\1-\\2)", x))
  } else {
    tab <- tab %>%
      select(-`**Characteristic**`) %>%
      rename(
        N = `**N**`,
        `N Event` = `**N Event**`
      )
    cols <- 3:ncol(tab)
    tab[cols] <- lapply(tab[cols], function(x) gsub("\\((\\d+)%?,\\s*(\\d+)%?\\)", "(\\1-\\2)", x))
  }


  tab1 <- flextable(tab) %>%
    bold(j = 1:ncol(tab), part = "header") %>%
    align(j = 2:ncol(tab), part = "all", align = "center")


  #-------------------------------------------------
  # Grey´s Test for difference
  #------------------------------------------------
  p_value <- surv_object$cmprsk$Tests[1, 2]
  p_value <- ifelse(p_value >= 0.1, round(p_value, 2), round(p_value, 3))
  p_value <- ifelse(p_value < 0.001, "< 0.001", p_value)


  #-----------------------------------------------------
  # Calculation of the plot
  #---------------------------------------------------
  cuminc_plot <- cuminc(formula = formula, data = data) %>%
    ggcuminc() +
    labs(
      x = x_title,
      y = y_title
    ) +
    add_risktable(risktable_stats = risk_table_elements) +
    ggtitle(title) +
    annotation_custom(
      grob = textGrob(paste0("p = ", p_value),
        x = p_placement[1], y = p_placement[2],
        hjust = 0, gp = gpar(col = "black")
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    scale_y_continuous(
      labels = scales::percent,
      breaks = y_breaks
    ) +
    scale_x_continuous(breaks = x_breaks) +
    # theme_minimal() +
    theme(
      plot.title = element_text(
        face = "bold", # makes it bold
        hjust = 0.5,
        size = text_size_title
      ),
      legend.position = legend_placement, legend.key.size = unit(0.7, "cm"),
      legend.background = element_rect(fill = alpha("blue", 0)),
      axis.text.x = element_text(
        color = "black",
        size = text_size
      ),
      axis.text.y = element_text(
        color = "black",
        size = text_size
      ),
      axis.title.x = element_text(color = "black", size = text_size),
      axis.title.y = element_text(color = "black", size = text_size),
      legend.text = element_text(size = text_size)
    ) +
    add_censor_mark()

  if (change_xscale == TRUE) {
    cuminc_plot <- cuminc_plot +
      scale_x_continuous(
        breaks = x_breaks,
        labels = function(x) round(x / change_factor)
      )
  }


  # adding the labels (median cumulative incidence or probability, if desired)

  if (show_label == "probability") {
    cuminc_plot <- cuminc_plot +
      scale_color_manual(
        values = colors,
        # labeling both competing risks
        labels = c(
          paste0(table_surv_prob$strata[1], " ", timepoint_label, ": ", table_surv_prob[1, 5]),
          paste0(table_surv_prob$strata[2], " ", timepoint_label, ": ", table_surv_prob[2, 5])
        )
      )
  } else if (show_label == "median") {
    cuminc_plot <- cuminc_plot +
      scale_color_manual(
        values = colors,
        # labeling both competing risks
        labels = c(
          paste0(median_table$outcome[1], ": Median ", median_table$time[1]),
          paste0(median_table$outcome[2], ": Median ", median_table$time[2])
        )
      )
  } else {
    cuminc_plot <- cuminc_plot +
      scale_color_manual(
        values = colors
      )
  }


  #------------------------------------------------------------------------------------
  # adding the ARD: cave: this is just for the main outcome, not for the competing event
  #------------------------------------------------------------------------------------
  cumprob_tidy <- surv_object$tidy

  arm_levels <- levels(data[[group]])

  cumprob_tidy <- cumprob_tidy %>%
    filter(outcome == main_outcome) %>%
    select(time, strata, estimate, std.error) %>%
    mutate(strata = factor(strata, levels = arm_levels))

  # I need common timepoints for both groups
  # --- 2. Helper: build step functions per stratum ---
  # KM curves start at S(0) = 1, so we prepend time = 0 with estimate = 1, se = 0
  build_stepfun <- function(df, value_col) {
    df <- df %>% arrange(time)
    times <- c(0, df$time)
    values <- c(1, df[[value_col]])
    if (value_col == "std.error") values[1] <- 0
    stepfun(times[-1], values, right = FALSE)
  }

  # --- 3. Split by strata ---
  df1 <- filter(cumprob_tidy, strata == arm_levels[1])
  df2 <- filter(cumprob_tidy, strata == arm_levels[2])

  est_fun1 <- build_stepfun(df1, "estimate")
  est_fun2 <- build_stepfun(df2, "estimate")
  se_fun1 <- build_stepfun(df1, "std.error")
  se_fun2 <- build_stepfun(df2, "std.error")

  # --- 4. Union of all time points (this "unifies" the time axis) ---
  all_times <- sort(unique(c(0, df1$time, df2$time)))

  # --- 5. Build the wide comparison table ---
  ard_df <- tibble(time = all_times) %>%
    mutate(
      estimate_1  = est_fun1(time),
      estimate_2  = est_fun2(time),
      se_1        = se_fun1(time),
      se_2        = se_fun2(time),
      diff        = estimate_1 - estimate_2,
      se_diff     = sqrt(se_1^2 + se_2^2), # independent strata -> variances add
      ci_lower    = diff - qnorm(0.975) * se_diff,
      ci_upper    = diff + qnorm(0.975) * se_diff
    ) %>%
    rename_with(~ paste0(., "_", arm_levels[1]), c(estimate_1, se_1)) %>%
    rename_with(~ paste0(., "_", arm_levels[2]), c(estimate_2, se_2))


  #------------------------------------------------
  # Extract the difference at xy time
  #-----------------------------------------------
  ard_diff_table <- map_dfr(time_vec_prob, function(t) {
    ard_df %>%
      arrange(time) %>%
      filter(time <= t) %>%
      filter(time == max(time)) %>%
      slice(1) %>% # in case of ties at the same max time
      mutate(
        timepoint = t,
        surv_prob = paste0(
          round(diff * 100), "% (",
          round(ci_lower * 100), "-",
          round(ci_upper * 100), ")"
        )
      ) %>%
      select(timepoint, time_used = time, diff, ci_lower, ci_upper, surv_prob)
  })

  ard_diff_table <- ard_diff_table %>%
    mutate(
      NNT = 1 / diff,
      NNT_LCI = min(1 / ci_upper, 1 / ci_lower),
      NNT_UCI = max(1 / ci_upper, 1 / ci_lower),
      NNT_CI = paste0(round(NNT), " (", round(NNT_LCI), "-", round(NNT_UCI), ")")
    )


  res_list <- list(
    "plot" = cuminc_plot,
    "probabilities" = table_surv_prob,
    "median" = median_table,
    "all_probabilities" = tab1,
    "ARD" = ard_diff_table
  )
  res_list
}
