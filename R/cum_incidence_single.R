#' Cumulative incidence plot for comepeting risks for one group
#'
#' @description
#' This function plots the cumulative incidence curves for two competing events.
#' Dispalyed are both curves (for both competing events) as well as optionally
#' median or probability for a certain time points plotted within the plot.
#'
#'
#' @param time A character string specifying the column name of the numeric time
#'             variable for the time-to-event endpoint. 0 = censored.
#' @param comp_risk_labels a character vector of the same length as the number of competing
#'                         risks with its names (e.g. death due to cancer, death due to corona)
#' @inheritParams km_grouped
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
#' @importFrom survival Surv
#'
#'
#'
comp_risk_single <-
  function(data,
           time,
           event,
           title = "",
           time_survival = 1,
           unit = "months",
           x_title = waiver(),
           y_title = waiver(),
           x_lim = NULL,
           y_lim = c(0, 1),
           x_breaks = waiver(),
           y_breaks = seq(0, 1, by = 0.2),
           colors = FALSE,
           # otpions: median, probability
           show_label = "none",
           text_size_title = 15,
           text_size = 12,
           legend_placement = c(0.48, 0.8),
           time_vec_prob = c(1, 2, 3),
           comp_risk_labels = NULL) {

    # Checking parameter input
    assert_numeric(data[[time]], lower = 0, any.missing = FALSE)
    assert_factor(data[[event]], any.missing = FALSE)

    time_sym <- sym(time)
    event_sym <- sym(event)
    data1 <- data %>%
      rename(
        time1 = !!time_sym,
        event1 = !!event_sym
      )
if(!is.null(comp_risk_labels) == TRUE) {
  data1$event1 <- factor(data1$event1, labels = comp_risk_labels)
}

    cr <- levels(data1$event1)[-1]
    print(cr)

    helper_df <- data.frame(outcome = cr)



    # Creation of survival object and the cumulative incidence for a given timepoint
    surv_object <- tidycmprsk::cuminc(Surv(time1, event1) ~ 1,
      data = data1, conf.type = "arcsin"
    )
    # adding the outcome (competing events) as factor again, becaus the factor coding is
    # deleted when computing the cuminc object
    surv_object$tidy$outcome <- factor(surv_object$tidy$outcome, levels = cr, labels = cr)

    # probability to survive (for both competing events) after a certain time point
    df_surv_prob <-
      surv_object %>%
      tbl_cuminc(
        times = time_survival, label_header = "**Months {time}**",
        outcomes = cr
      )

    df_surv_prob_tidy <- df_surv_prob$tidy

    # median probability
    median_table <- surv_object$tidy %>%
      group_by(.data$outcome) %>%
      filter(.data$estimate >= 0.5) %>%
      group_by(.data$outcome) %>%
      slice_head() %>%
      select(.data$outcome, .data$time) %>%
      arrange(factor(.data$outcome, levels = cr))

median_table <- merge(helper_df, median_table, by = "outcome", all.x = TRUE)



    options("ggsurvfit.switch-color-linetype" = TRUE)
    # specifying colors, if they were not specified in the parameters
    if (is.logical(colors)) {
      colors <- colorRamps::primary.colors(1)
    }

    plot <- surv_object %>%
      ggcuminc(outcome = cr) +
      labs(
        x = x_title,
        y = y_title
      ) +
      ggtitle(title) +
      add_risktable() +
      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      scale_y_continuous(
        labels = scales::percent,
        breaks = y_breaks
      ) +
      scale_x_continuous(breaks = x_breaks) +
      theme(
        plot.title = element_text(size = text_size_title, hjust = 0.5, face = "bold"),
        legend.position = legend_placement,
        legend.key.size = unit(0.7, "cm"),
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
      plot <- plot +
        scale_color_manual(
          values = colors,
          # labeling both competing risks
          labels = c(
            paste0(cr[1], ": ", time_survival, " ", df_surv_prob_tidy[1, 4]),
            paste0(cr[2], ": ", time_survival, " ", df_surv_prob_tidy[2, 4])
          )
        )
    } else if (show_label == "median") {
      plot <- plot +
        scale_color_manual(
          values = colors,
          # labeling both competing risks
          labels = c(
            paste0(median_table$outcome[1], ": median survival ", round(median_table$time[1], 2)),
            paste0(median_table$outcome[2], ": median survival ", round(median_table$time[2], 2))
          )
        )
    } else {
      plot <- plot +
        scale_color_manual(
          values = colors
        )
    }

    tab_cum_incidence_prob <- surv_object %>%
      tbl_cuminc(
        times = time_vec_prob, label_header = "**Months {time}**",
        outcomes = cr
      ) %>%
      add_nevent() %>%
      add_n() %>%
      add_nevent(location = c("label", "level")) %>%
      add_n(location = c("label", "level"))

    result <- list(plot, tab_cum_incidence_prob, median_table)
    names(result) <- c("plot", "cumulative incidence probability", "median table")
    options("ggsurvfit.switch-color-linetype" = FALSE)
    result
  }
