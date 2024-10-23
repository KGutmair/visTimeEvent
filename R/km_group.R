#' Plotting Kaplan-Meier curves stratified by groups
#'
#'@description This function is plotting Kaplan-Meier curves stratified by groups.
#'            the `survfit2()` function is used to plot the Kaplan-Meier curves.
#'            Additionally, p-values, median survival and probability of survival
#'            after a certain time is computed and added as information to the plot,
#'            if desired.
#'
#'
#' @param data A `data.frame` containing the necessary columns: time, event, and group.
#' @param time A character string specifying the column name of the numeric time
#'             variable for the time-to-event endpoint.
#' @param event A character string specifying the column name of the event indicator.
#'              Should be a numeric variable where 1 represents an event, and 0
#'              represents no event or right-censoring.
#' @param group A character string specifying the grouping variable used for stratifying
#'              the Kaplan-Meier curves.
#' @param title A character string providing the title for the Kaplan-Meier plot.
#' @param time_survival A numeric value representing the time point at which the
#'                      event pobability is calculated. Default is 1.
#' @param unit A character string indicating the unit for `time_survival` (e.g., "months").
#'             Default is "months".
#' @param test A character string specifying the type of statistical test for
#'             comparing survival curves. Options are `log-rank` or `wilcoxon`. Default is `log-rank`.
#' @param x_title A character string specifying the label for the x-axis.
#' @param y_title A character string specifying the label for the y-axis.
#' @param endpoint A character string providing the name of the endpoint.
#' @param x_lim A numeric vector of length two defining the limits of the x-axis.
#' @param y_lim A numeric vector of length two defining the limits of the y-axis.
#'              Default is `c(0, 1)`.
#' @param x_breaks A numeric vector specifying the break points for the x-axis labels.
#' @param y_breaks A numeric vector specifying the break points for the y-axis labels.
#'                 Default is from 0 to 1 with steps of 0.2.
#' @param colors A character vector specifying the colors for the curves. If not provided,
#'               random colors will be assigned.
#' @param show_label A logical or character value indicating whether to display median survival time of event probabilites
#'                   labels within the plot. Options are `FALSE` (default), `median`, or `probability`.
#' @param text_size_title A numeric value specifying the font size of the title.
#' @param text_size A numeric value specifying the font size of the text elements.
#' @param show_p_values logical, indicating whether the p-values be displayed in the plot?. Default = `TRUE`

#' @param p_placement numeric vector of length 2 indicating the position of the the p-value.
#'                    The first number correspondend to the hjust placement, the second to the vjust placement
#' @param legend_placement numeric vector of length 2 indicating the legend position. Default is c(0.6, 0.2).
#'
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


km_grouped <-
  function(data,
           time,
           event,
           group,
           title = "",
           time_survival = 1,
           unit = "months",
           test = "log-rank",
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
           show_p_values = TRUE,
           p_placement = c(0.05, 0.15),
           legend_placement = c(0.6, 0.2)) {


    # Assert data frame structure
    assert_data_frame(data, min.rows = 1, min.cols = 1)

    # Check that time, event, and group are valid columns in the data
    variables <- c(time, event, group)
    lapply(variables, function(var) assert_vector(var, any.missing = FALSE, len = 1))
    if (any(!variables %in% colnames(data))) {
      stop("At least one of the variables time, event, group cannot be found in the data.")
    }

    # Assertions for strings and numeric values
    assert_string(title)
    assert_string(endpoint)
    assertNumeric(data[[time]], lower = 0, any.missing = FALSE)
    assert_integer(data[[event]], lower = 0, any.missing = FALSE)
    assertSubset(data[[event]], choices = c(0, 1), empty.ok = FALSE)
    assertNumeric(time_survival, lower = 0, any.missing = TRUE, len = 1)
    assertNumeric(p_placement, any.missing = FALSE, len = 2)
    assertNumeric(legend_placement, any.missing = FALSE, len = 2)


    # Check if the number of colors matches the number of groups
    if (!identical(colors, FALSE) && length(colors) != length(unique(data[[group]]))) {
      warning("Number of colors does not correspond to the number of groups.")
    }



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
    print(paste0("- significance test: ", test))
    print(paste0("- time-to-event probability of ", time_survival, " ", unit))

    group_sym <- sym(group)
    time_sym <- sym(time)
    event_sym <- sym(event)

    data1 <- data %>%
      rename(
        group1 = !!group_sym,
        time1 = !!time_sym,
        event1 = !!event_sym
      )

    # specifying colors, if they were not specified in the parameters
    if (is.logical(colors)) {
      n <- length(unique(data1$group1))
      colors <- colorRamps::primary.colors(n + 1)[-1]
    }

    # calling the function whose output are a data frame with
    # the median survival times and survival probabilities + its CIs
    km_table1 <-
      median_probability(
        data = data,
        time = time,
        event = event,
        group = group,
        time_survival = time_survival,
        test = test
      )

    # table containing the information for median for every treatment arm
    median_table <- km_table1[[1]]
    median_table$median <- format(round(as.numeric(median_table$median), 1), nsmall = 1)
    # table containing the information to xy% probability
    surv_prob_table <- km_table1[[2]]



    # plotting the survival curves
    km_plot <-
      survfit2(Surv(time = time1, event = event1) ~ group1,
               data = data1
      ) %>%
      ggsurvfit() +
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
        # legend position: between 0 and 1, relative position, not bound to
        # absolute values, so no dependence on the plot's scale
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
      add_censor_mark() +
      add_risktable()

    if(show_p_values == TRUE) {
      km_plot <- km_plot +
        annotation_custom(
          grob = textGrob(paste0("p = ", median_table$p_value[1]),
                          x = p_placement[1], y = p_placement[2],
                          hjust = 0, gp = gpar(col = "black")),
          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
        )
    }

    # setting labels, if desired
    if (show_label == "probability") {
      label_vec <- c()
      for (i in seq_len(nrow(surv_prob_table))) {
        label_vec[i] <- paste(surv_prob_table$names[i], time_survival, unit, endpoint, " probability  = ", surv_prob_table$surv_prob[i], sep = " ")
      }
      km_plot <- km_plot +
        scale_color_manual(
          values = colors,
          labels = label_vec
        )
    } else if (show_label == "median") {
      label_vec <- c()
      for (i in seq_len(nrow(median_table))) {
        label_vec[i] <- paste(surv_prob_table$names[i], " median  = ", median_table$median_CI[i], sep = " ")
      }
      km_plot <- km_plot +
        scale_color_manual(
          values = colors,
          labels = label_vec
        )
    } else {
      km_plot <- km_plot +
        scale_color_manual(
          values = colors
        )
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
median_probability <-
  function(data,
           time,
           event,
           group,
           # 5 year- probbaility of survival
           time_survival = 1,
           test = "log-rank") {


    group_sym <- sym(group)
    time_sym <- sym(time)
    event_sym <- sym(event)
    data1 <- data %>%
      rename(
        group1 = !!group_sym,
        time1 = !!time_sym,
        event1 = !!event_sym
      )

    res <-
      summary(survfit(Surv(time = time1, event = event1) ~ group1, data = data1),
              times = time_survival, extend = TRUE, conf.type = "arcsin"
      )

    # table of survival probabilities
    survival_prob <- data.frame(res$time, res$n.risk, res$n.event, res$surv, res$lower, res$upper)
    survival_prob <- survival_prob %>%
      mutate(names = unique(data1$group1)) %>%
      mutate(
        survival_prop = round(.data$res.surv * 100, 0),
        LCI = round(.data$res.lower * 100, 0),
        UCI = round(.data$res.upper * 100, 0)
      ) %>%
      mutate(surv_prob = paste0(.data$survival_prop, "% (", .data$LCI, " - ", .data$UCI, ")")) %>%
      select(-c(.data$res.surv, .data$res.lower, .data$res.upper)) %>%
      relocate(.data$names, .data$res.time, .data$res.n.risk, .data$res.n.event, .data$surv_prob)

    type <- switch(test,
                   "log-rank" = 0,
                   "wilcoxon" = 1
    )
    res_test <-
      survdiff(Surv(time1, event1) ~ group1, data = data1, rho = type)

    # Table of median
    median_tab <- as.data.frame(res$table)
    median_tab <- median_tab %>%
      mutate_at(c("median", "0.95LCL", "0.95UCL"), function(x) format(round(x, 2), nsmall = 2)) %>%
      mutate(
        time = rep(as.character(time_survival), times = 2),
        pvalue = ifelse(res_test$pvalue >= 0.1,
                        format(round(res_test$pvalue, 2), nsmall = 2),
                        format(round(res_test$pvalue, 3), nsmall = 3)
        ),
        p_value = ifelse(.data$pvalue < 0.001, "< 0.001", .data$pvalue),
        median_CI = paste0(.data$median, " (", .data$`0.95LCL`, " - ", .data$`0.95UCL`, ")")
      ) %>%
      select(-c(.data$rmean, .data$`se(rmean)`, .data$`0.95LCL`,
                .data$`0.95UCL`, .data$pvalue))


    res_list <- list(median_tab, survival_prob)
    names(res_list) <- c("table of median survival", "table of survival probability")
    res_list
  }
