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
                              time_survival = 1,
                              unit = "months",
                              endpoint = "",
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
                              time_vec_prob = c(1, 2, 3)) {

  #########################
  # Checking input variables
  ###########################

  assert_numeric(data[[time]], lower = 0, any.missing = FALSE)
  assert_factor(data[[event]], any.missing = FALSE)
  assert_factor(data[[group]], any.missing = FALSE, max.levels = 2, min.levels = 2)

  # Creation of survival object and the cumulative incidence for a given timepoint

  formula <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", group))
print(formula)
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

  # Survival probability after a certain amount of time
  tbl <-
    surv_object %>%
    tbl_cuminc(
      times = time_survival, label_header = unit,
      outcomes = cr
    )

  table_surv_prob <- tbl$tidy

  # Median survival probability

  # Median probability

  median_table <- surv_object$tidy %>%
    filter(.data$estimate >= 0.5) %>%
    group_by(.data$outcome) %>%
    arrange(.data$outcome, .data$time) %>%
    slice_head() %>%
    select(.data$time)



  # Output table for survival probabilities at multiple timepoints
  tab <- surv_object %>%
    tbl_cuminc(
      times = time_vec_prob, label_header = "**Months {time}**",
      outcomes = cr
    ) %>%
    add_nevent() %>%
    add_n() %>%
    add_nevent(location = c("label", "level")) %>%
    add_n(location = c("label", "level"))



  # Grey´s Test for difference
  p_value <- surv_object$cmprsk$Tests[1, 2]
  p_value <- ifelse(p_value >= 0.1, round(p_value, 2), round(p_value, 3))
  p_value <- ifelse(p_value < 0.001, "< 0.001", p_value)



  cuminc_plot <- cuminc(formula = formula, data = data) %>%
    ggcuminc() +
    labs(
      x = x_title,
      y = y_title
    ) +
    add_risktable() +
    ggtitle(title) +
    annotation_custom(
      grob = textGrob(paste0("p = ", p_value),
                      x = p_placement[1], y = p_placement[2],
                      hjust = 0, gp = gpar(col = "black")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    scale_y_continuous(
      labels = scales::percent,
      breaks = y_breaks
    ) +
    scale_x_continuous(breaks = x_breaks) +
    theme(
      plot.title = element_text(size = text_size_title),
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


  # adding the labels (median cumulative incidence or probability, if desired)

  if (show_label == "probability") {
    cuminc_plot <- cuminc_plot +
      scale_color_manual(
        values = colors,
        # labeling both competing risks
        labels = c(
          paste0(table_surv_prob$strata[1], " ", endpoint, ": ", time_survival, " ", unit, ": ", table_surv_prob[1, 5]),
          paste0(table_surv_prob$strata[2], " ", endpoint, ": ", time_survival, " ", unit, ": ", table_surv_prob[2, 5])
        )
      )
  } else if (show_label == "median") {
    cuminc_plot <- cuminc_plot +
      scale_color_manual(
        values = colors,
        # labeling both competing risks
        labels = c(
          paste0(median_table$outcome[1], ": median survival ", median_table$time[1]),
          paste0(median_table$outcome[2], ": median survival ", median_table$time[2])
        )
      )
  } else {
    cuminc_plot <- cuminc_plot +
      scale_color_manual(
        values = colors
      )
  }



  res_list <- list("plot" = cuminc_plot, "probabilities" = table_surv_prob, "median" = median_table)
  res_list
}
