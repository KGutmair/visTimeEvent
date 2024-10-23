#' Plotting Kaplan-Meier curves for one group
#'
#'@description
#'This function is plotting Kaplan-Meier curves stratified by groups.
#'            the `survfit2()` function is used to plot the Kaplan-Meier curves.
#'            Additionally, p-values, median survival and probability of survival
#'            after a certain time is computed and added as information to the plot,
#'            if desired.
#'
#' @inheritParams km_grouped

#' @return list containing
#' a, the ggplot object of the Kaplan-Meier plot
#' b, a data frame containing median survival + CI,
#' c, a data frame containing the event probability after a certain time point + CI

#' @export
#'
#' @import checkmate
#' @importFrom dplyr %>% rename
#' @importFrom colorRamps primary.colors
#' @import  ggsurvfit
#' @importFrom survival Surv
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom scales percent
#' @import grid


km_single <-
  function(data,
           time,
           event,
           title = "",
           time_survival = 1,
           unit = "months",
           x_title = waiver(),
           y_title = waiver(),
           endpoint = "",
           x_lim = NULL,
           y_lim = c(0, 1),
           x_breaks = waiver(),
           y_breaks = seq(0, 1, by = 0.2),
           colors = FALSE,
           # otpions: median, probability
           show_label = "none",
           text_size_title = 15,
           text_size = 12,
           legend_placement = c(0.5, 0.8)) {


    # Assert data frame structure
    assert_data_frame(data, min.rows = 1, min.cols = 1)

    # Check that time, event, and group are valid columns in the data
    variables <- c(time, event)
    lapply(variables, function(var) assert_vector(var, any.missing = FALSE, len = 1))
    if (any(!variables %in% colnames(data))) {
      stop("At least one of the variable names cannot be found in the data.")
    }

    # Assertions for strings and numeric values
    assert_string(title)
    assert_string(endpoint)
    assertNumeric(data[[time]], lower = 0, any.missing = TRUE)
    assertNumeric(data[[event]], any.missing = TRUE)
    assertSubset(data[[event]], choices = c(0, 1), empty.ok = TRUE)
    assertNumeric(time_survival, lower = 0, any.missing = TRUE, len = 1)
    assertNumeric(legend_placement, any.missing = FALSE, len = 2)


    # Assert x_lim and y_lim must have exactly two values
    if(!is.null(x_lim) == TRUE) {

      assert_numeric(x_lim, len = 2, any.missing = FALSE)
    }

    if(!is.null(y_lim) == TRUE) {
      assert_numeric(y_lim, len = 2, any.missing = FALSE)
    }

    is_waiver <- function(x) {
      inherits(x, "waiver")
    }

    # Assertions for numeric x_breaks, y_breaks
    assertNumeric(y_breaks, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 2)

    if(is_waiver (x_breaks) == FALSE) {
      assertNumeric(x_breaks, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 2)
    }

    if(is_waiver(x_title) == FALSE) {
      assert_string(x_title)
    }

    if(is_waiver(y_title) == FALSE) {
      assert_string(y_title)
    }

    # Assert show_label is a valid option
    assert_character(show_label, any.missing = FALSE, max.len = 1)
    assertSubset(show_label, choices = c("none", "median", "probability"))

    # Assertions for text sizes
    lapply(list(text_size_title, text_size),
           function(x) assertNumeric(x, lower = 0, len = 1, any.missing = FALSE))

    print("Following options were chosen: ")
    print(paste0("- time-to-event probability of ", time_survival, " ", unit))

    time_sym <- sym(time)
    event_sym <- sym(event)

    data1 <- data %>%
      rename(
        time1 = !!time_sym,
        event1 = !!event_sym
      )

    # specifying colors, if they were not specified in the parameters
    if (is.logical(colors)) {
      colors <- colorRamps::primary.colors(1)
    }

    # calling the function whose output are a data frame with
    # the median survival times and survival probabilities + its CIs
    km_table1 <-
      median_probability_single(
        data = data,
        time = time,
        event = event,
        time_survival = time_survival )

    # table containing the information for median for every treatment arm
    median_table <- km_table1[[1]]
    # table containing the information to xy% probability
    surv_prob_table <- km_table1[[2]]



    # plotting the survival curves
    km_plot <-
      survfit2(Surv(time = time1, event = event1) ~ 1,
               data = data1
      ) %>%
      ggsurvfit(color = colors) +
      labs(
        x = x_title,
        y = y_title
      ) +
      ggtitle(title) +
      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      scale_y_continuous(
        labels = scales::percent,
        breaks = y_breaks
      ) +
      scale_x_continuous(breaks = x_breaks) +
      theme(
        plot.title = element_text(size = text_size_title),
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
      add_censor_mark(color = colors) +
      add_risktable()


    # setting labels, if desired
    if (show_label == "probability") {
        label_vec <- paste(time_survival, unit, endpoint, " probability  = ",
                           surv_prob_table$surv_prob, sep = " ")

      km_plot <- km_plot +
        annotation_custom(
          grob = textGrob(label_vec,
                          legend_placement[1], y = legend_placement[2],
                          hjust = 0, gp = gpar(col = "black")),
          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
        )
    } else if (show_label == "median") {
      label_vec <- paste( " median ", endpoint, ": ", median_table$median_CI, sep = " ")

      km_plot <- km_plot +
      annotation_custom(
        grob = textGrob(label_vec,
                        legend_placement[1], y = legend_placement[2],
                        hjust = 0, gp = gpar(col = "black")),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      )
    } else {
      km_plot
    }

    km_result <- list(km_plot, median_table, surv_prob_table)
    names(km_result) <- c("plot", "table of median survival", "table of survival probability")
    km_result
  }








#' Calculation of median survival time and survival probabiliy of time-to-event endpoints
#'
#' @inheritParams km_grouped
#' @importFrom dplyr %>% mutate select relocate rename mutate_at
#' @importFrom survival Surv survfit survdiff
#' @importFrom rlang .data
#'
#' @noRd
#'
#' @return a list with two data frames, one containing the information for median event time
#'        and one for the probability of an event after a certain time
#'
#'
#'
#'
median_probability_single <-
  function(data,
           time,
           event,
           # 5 year- probbaility of survival
           time_survival = 1) {


    time_sym <- sym(time)
    event_sym <- sym(event)
    data1 <- data %>%
      rename(
        time1 = !!time_sym,
        event1 = !!event_sym
      )

    res <-
      summary(survfit(Surv(time = time1, event = event1) ~ 1, data = data1),
              times = time_survival, extend = TRUE, conf.type = "arcsin"
      )



    # table of survival probabilities
    survival_prob <- data.frame(res$time, res$n.risk, res$n.event, res$surv, res$lower, res$upper)
    survival_prob <- survival_prob %>%
      mutate(names = "whole group") %>%
      mutate(
        survival_prop = round(.data$res.surv * 100, 0),
        LCI = round(.data$res.lower * 100, 0),
        UCI = round(.data$res.upper * 100, 0)
      ) %>%
      mutate(surv_prob = paste0(.data$survival_prop, "% (", .data$LCI, " - ", .data$UCI, ")")) %>%
      select(-c(.data$res.surv, .data$res.lower, .data$res.upper)) %>%
      relocate(.data$names, .data$res.time, .data$res.n.risk, .data$res.n.event, .data$surv_prob)


    # Table of median
    median_tab <- as.data.frame(t(res$table))
    median_tab <- median_tab %>%
      mutate_at(c("median", "0.95LCL", "0.95UCL"), function(x) format(round(x, 2), nsmall = 2)) %>%
      mutate(
        time = rep(as.character(time_survival), times = 1),
        median_CI = paste0(.data$median, " (", .data$`0.95LCL`, " - ", .data$`0.95UCL`, ")")
      ) %>%
      select(-c(.data$rmean, .data$`se(rmean)`, .data$`0.95LCL`,
                .data$`0.95UCL`))


    res_list <- list(median_tab, survival_prob)
    names(res_list) <- c("table of median survival", "table of survival probability")
    res_list
  }
