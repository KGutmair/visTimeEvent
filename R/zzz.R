#' Kaplan Meier plots for two groups with IPTW
#'
#' @description
#' This function plots Kaplan-Meier curves and calculates median survival time,
#' the probability of getting an event after a certain period and p values of
#' either the log-rank test or Wilcoxon test for two strata
#' Note, that the log.rank test is only possible for two groups
#'
#' @inheritParams km_grouped
#'
#' @param weight_col Character string indicating the column name containing the weights.
#' @param group_name Names of the groups to be compared (only used for the labels of the risk table).
#'                   The names are plotted in the corresponding order: first element of the vector
#'                   is printed in the risk table in first position, the second element of the
#'                   vector is printed above. In the risl table, the order of the risk table
#'                   corresponds to the order of the levels of the group variable: the first level
#'                   is printed at first, the second level above
#' @param risk_table Logical, default: TRUE. Should a risk table be displayed?
#' @param custom_colors Character string specifying the colors for the risk table. If not provided, default colors will be automatically selected.
#' @param color_curves Character string specifying the colors for the survival curves.
#' @param show_weighted_n Logical, default: TRUE. Should the weighted or unweighted number of observations at risk be displayed?
#'
#'
#' @return
#' # list containing
#' a, the ggplot object of the Kaplan-Meier plot
#' b, a data frame containing median survival, the probability to get an event after
#' the prespecified period and its confidence interval, the p value of the log rank or wilcoxon test
#'
#' @export
#'
#' @import adjustedCurves
#' @importFrom dplyr %>% rename mutate select group_by filter relocate
#' @importFrom RISCA ipw.log.rank
#' @importFrom colorRamps primary.colors
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom scales percent
#' @import grid
#'
#' km_grouped_weighted <- function(data,
#'                                 time,
#'                                 event,
#'                                 group,
#'                                 weight_col,
#'                                 group_name,
#'                                 endpoint,
#'                                 title = "",
#'                                 x_title = "Time",
#'                                 y_title = "Adjusted Survival Probability",
#'                                 time_survival = 1,
#'                                 risk_table = TRUE,
#'                                 x_lim = NULL,
#'                                 y_lim = c(0, 1),
#'                                 x_breaks = waiver(),
#'                                 y_breaks = seq(0, 1, by = 0.2),
#'                                 custom_colors = NULL,
#'                                 text_size_title = 15,
#'                                 text_size = 12,
#'                                 show_label = FALSE,
#'                                 unit = "months",
#'                                 color_curves = FALSE,
#'                                 # options: n_weighted, n_unweighted
#'                                 show_weighted_n = TRUE,
#'                                 show_p_values = TRUE,
#'                                 p_placement = c(0.05, 0.15),
#'                                 legend_placement = c(0.6, 0.25)) {
#'
#'   ##############
#'   # Checking input variables
#'   ##############
#'   assert_numeric(data[[time]], lower = 0, any.missing = FALSE)
#'   assert_numeric(data[[event]], lower = 0, any.missing = FALSE)
#'   assertSubset(data[[event]], choices = c(0, 1), empty.ok = FALSE)
#'   assert_factor(data[[group]], any.missing = FALSE, min.levels = 2)
#'
#'   group_sym <- sym(group)
#'   time_sym <- sym(time)
#'   event_sym <- sym(event)
#'
#'
#'   data1 <- data %>%
#'     rename(
#'       group1 = !!group_sym,
#'       time1 = !!time_sym,
#'       event1 = !!event_sym
#'     )
#'
#'
#'
#'   adjsurv <- adjustedsurv(
#'     data = data1,
#'     variable = "group1",
#'     ev_time = "time1",
#'     event = "event1",
#'     method = "iptw_km",
#'     treatment_model = data1[[weight_col]],
#'     conf_int = TRUE
#'   )
#'
#'
#'   time_table <- adjsurv$adj
#'
#'   ### Median survival and weighted logrank test
#'   median_surv <- adjusted_surv_quantile(adjsurv, p = 0.5, conf_int = TRUE)
#'
#'   # the log-rank test can only handle group variables with numeric 0/1
#'   # or factor 0/1
#'
#'   data1$group_num <- as.numeric(data1$group1) -1
#'   test1 <- ipw.log.rank(data1$time1, data1$event1, data1$group_num,
#'                         weights = data1[[weight_col]])
#'
#'
#'
#'   median_surv$p_value <- c(test1$p.value, rep(NA, times = length(unique(data1$group1)) - 1))
#'   median_surv <- median_surv %>%
#'     mutate(
#'       surv_time = format(round(.data$q_surv, 2), nsmall = 2),
#'       ci_lower = format(round(.data$ci_lower, 2), nsmall = 2),
#'       ci_upper = format(round(.data$ci_upper, 2), nsmall = 2),
#'       CI_95 = paste0(.data$ci_lower, " - ", .data$ci_upper),
#'       p_value = ifelse(.data$p_value >= 0.1, format(round(.data$p_value, 2),
#'                                                     nsmall = 2
#'       ), format(round(.data$p_value, 3), nsmall = 3)),
#'       p_value = ifelse(.data$p_value < 0.001, "< 0.001", .data$p_value)
#'     ) %>%
#'     select(.data$group, .data$surv_time, .data$CI_95, .data$p_value)
#'
#'
#'
#'
#'   # 5 year survival probability
#'
#'   surv_probability <- time_table %>%
#'     group_by(.data$group) %>%
#'     filter(.data$time < time_survival) %>%
#'     filter(.data$time == max(.data$time))
#'
#'   surv_probability <- surv_probability %>%
#'     mutate(
#'       surv = round(.data$surv * 100, 0),
#'       ci_lower = round(.data$ci_lower * 100, 0),
#'       ci_upper = round(.data$ci_upper * 100, 0)
#'     ) %>%
#'     mutate(surv_prob = paste0(.data$surv, "% (", .data$ci_lower, " - ", .data$ci_upper, ")")) %>%
#'     select(-c(.data$se, .data$ci_lower, .data$ci_upper)) %>%
#'     relocate(.data$group, .data$time, .data$surv)
#'
#'
#'   ### Survival curve
#'
#'   # specifying colors, if they were not specified in the parameters
#'   if (is.logical(color_curves)) {
#'     n <- length(unique(data1$group1))
#'     color_curves <- primary.colors(n + 1)[-1]
#'   }
#'
#'   if(is.null(custom_colors)) {
#'     custom_colors <- rep("black", times = 2*length(unique(data1$group1)))
#'   }
#'
#'
#'   ### Survival curve
#'   #### calcualting the survival curve with the normal surv object, because
#'   # the adjusted curv object does only calcuate survival probabilities till
#'   # the last event. Meaning, from the curve, I cannot see the whole
#'   # observation period
#'   surv_object <- survfit2(Surv(time1, event1) ~ group1,
#'                           weights = data1[[weight_col]], data = data1)
#'
#'
#'   survival.data <- data.frame("time" = surv_object$time,
#'                               "surv" = surv_object$surv,
#'                               "group" = c(rep(surv_probability$group[1], times = surv_object$strata[1]),
#'                                           rep(surv_probability$group[2], times = surv_object$strata[2]),
#'                                           rep(surv_probability$group[3], times = surv_object$strata[3])),
#'                               "ncensor" = surv_object$counts[, 3])
#'
#'   km_plot <- surv_object %>%
#'     ggsurvfit() +
#'     # adding censor marks, becaus the ggsurfit does not want to plot them
#'     geom_point(data = survival.data[survival.data$ncensor > 0, ],
#'                aes(x = time, y = surv, color = group),
#'                shape = 3, size = 3) +
#'     labs(
#'       x = x_title,
#'       y = y_title
#'     ) +
#'     ggtitle(title) +
#'     coord_cartesian(xlim = x_lim, ylim = y_lim) +
#'     scale_y_continuous(
#'       labels = scales::percent,
#'       breaks = y_breaks
#'     ) +
#'     scale_x_continuous(breaks = x_breaks) +
#'     theme(
#'       plot.title = element_text(size = text_size_title),
#'       legend.position = legend_placement, legend.key.size = unit(0.7, "cm"),
#'       legend.background = element_rect(fill = alpha("blue", 0)),
#'       axis.text.x = element_text(
#'         color = "black",
#'         size = text_size
#'       ),
#'       axis.text.y = element_text(
#'         color = "black",
#'         size = text_size
#'       ),
#'       axis.title.x = element_text(color = "black", size = text_size),
#'       axis.title.y = element_text(color = "black", size = text_size),
#'       legend.text = element_text(size = 12)
#'     )
#'   if (show_p_values == TRUE) {
#'     km_plot <- km_plot +
#'       # annotate(
#'       #   geom = "text",
#'       #   x = -Inf,
#'       #   y = -Inf,
#'       #   hjust = p_placement[1],
#'       #   vjust = p_placement[2],
#'       #   label = paste0("p = ", median_surv$p_value[1]),
#'       #   color = "black", size = 5
#'       # )
#'       annotation_custom(
#'         grob = textGrob(paste0("p = ", median_surv$p_value[1]),
#'                         x = p_placement[1], y = p_placement[2],
#'                         hjust = 0, gp = gpar(col = "black")),
#'         xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
#'       )
#'   }
#'
#'
#'   # # setting labels, if desired
#'   if (show_label == "probability") {
#'     label_vec <- c()
#'     for (i in seq_len(nrow(surv_probability))) {
#'       label_vec[i] <- paste(surv_probability$group[i], time_survival, unit, endpoint, " probability  = ", surv_probability$surv_prob[i], sep = " ")
#'     }
#'     km_plot <- km_plot +
#'       scale_color_manual(
#'         values = color_curves,
#'         labels = label_vec
#'       )
#'   } else if (show_label == "median") {
#'     label_vec <- c()
#'     for (i in seq_len(nrow(median_surv))) {
#'       label_vec[i] <- paste(median_surv$group[i], " median  = ", median_surv$median_CI[i], sep = " ")
#'     }
#'     km_plot <- km_plot +
#'       scale_color_manual(
#'         values = color_curves,
#'         labels = label_vec
#'       )
#'   } else {
#'     km_plot <- km_plot +
#'       scale_color_manual(
#'         values = color_curves
#'       )
#'   }
#'
#'   # if x_lim is nor defined, then this will be calculated by waiver()
#'   # I need to exreact these values for plotting the risk table
#'   x_lim_present <- TRUE
#'   if(is.null(x_lim)) {
#'     plot_data <- ggplot_build(km_plot)
#'
#'     x_breaks <- plot_data$layout$panel_params[[1]]$x$breaks
#'     x_lim <- plot_data$layout$panel_params[[1]]$x.range
#'     x_lim_present <- FALSE
#'     print(c("km plot", x_breaks))
#'     print(c("km plot", x_lim))
#'
#'   }
#'
#'
#'
#'
#'   if (risk_table) {
#'     risk_tab <- plot_risk_table.groups(
#'       data = data,
#'       x_breaks = x_breaks,
#'       x_lim = x_lim,
#'       variable = group,
#'       ev_time = time,
#'       weight_col = weight_col,
#'       event = event,
#'       custom_colors = custom_colors,
#'       group_name = group_name,
#'       show_weighted_n = show_weighted_n,
#'       x_lim_present = x_lim_present
#'     )
#'
#'     km_plot <- cowplot::plot_grid(km_plot, risk_tab,
#'                                   ncol = 1, align = "v",
#'                                   rel_heights = c(3, 1.2), axis = "l"
#'     )
#'   }
#'
#'
#'   result <- list(km_plot, median_surv, surv_probability)
#'   result
#' }
#'
#'
#'
#' #' This function plots with ggplot a risk table
#' #' @inheritParams km_grouped
#'
#' #' @param variable character string indicating the stratification variable
#' #' @param ev_time A character string specifying the column name of the numeric time
#' #'             variable for the time-to-event endpoint.
#' #' @param group_name Names of the groups to be compared (used for the risk table).
#' #' @param show_weighted_n Logical, default: TRUE. Should the weighted or unweighted
#' #' number of observations at risk be displayed?
#' #' @param weight_col Character string indicating the column name containing the weights.
#' #' @param text_size a number inidicating the font size of the text
#' #' @param text_alpha a number inidicating the font size of the text
#' #' @param text_color a character string indicating the color of the text in the risk table
#' #' @param text_family a character string specifying the text family
#' #' @param text_fontface a character string specifying the fontface
#' #' @param color_groups logical TRUE/FALSE, default is TRUE
#' #' @param reverse_order should the risk table be printed in the reversed order,
#' #' logical TRUE/FALSE, default is TRUE
#' #' @param custom_colors character string specifying the colors
#' #' @param vjust numeric value
#' #'
#' #' @return a ggplot object with the risk table
#' #' @noRd
#' #'
#' #' @import ggplot2
#' #'
#' plot_risk_table.groups <- function(data,
#'                                    x_breaks,
#'                                    x_lim,
#'                                    variable,
#'                                    ev_time,
#'                                    group_name,
#'                                    show_weighted_n = TRUE,
#'                                    weight_col,
#'                                    event,
#'                                    text_size = 3.5,
#'                                    text_alpha = 1,
#'                                    text_color = "black",
#'                                    text_family = "sans",
#'                                    text_fontface = "plain",
#'                                    color_groups = TRUE,
#'                                    reverse_order = TRUE,
#'                                    custom_colors,
#'                                    vjust = 5,
#'                                    x_lim_present) {
#'   plotdata <- get_risk_table.groups(
#'     times = x_breaks, data = data,
#'     variable = variable,
#'     ev_time = ev_time,
#'     weights = data[[weight_col]],
#'     event = event,
#'     show_weighted_n = show_weighted_n
#'   )
#'
#'   mapping <- ggplot2::aes(
#'     x = .data$time,
#'     y = .data$group,
#'     color = .data$group1,
#'     label = .data$value
#'   )
#'
#'   if (!color_groups) {
#'     mapping$colour <- NULL
#'   }
#'
#'   if (color_groups) {
#'     main_geom <- ggplot2::geom_text(
#'       size = text_size, alpha = text_alpha,
#'       family = text_family, fontface = text_fontface
#'     )
#'   } else {
#'     main_geom <- ggplot2::geom_text(
#'       size = text_size, alpha = text_alpha,
#'       family = text_family, fontface = text_fontface,
#'       color = text_color
#'     )
#'   }
#'
#'   p <- ggplot2::ggplot(data = plotdata, mapping) +
#'     main_geom +
#'     ggplot2::scale_y_continuous(breaks = c(0, 0.3, 1, 1.3), labels = c(
#'       "Event", "At risk",
#'       "Event", "At risk"
#'     )) +
#'     annotate(
#'       geom = "text",
#'       hjust = -Inf,
#'       x = 3.8,
#'       y = 0.6,
#'       label = group_name[2],
#'       color = "black", size = 3.5
#'     ) +
#'     annotate(
#'       geom = "text",
#'       hjust = -Inf,
#'       x = 5.2,
#'       y = 1.6,
#'       label = group_name[1],
#'       color = "black", size = 3.5
#'     ) +
#'
#'     # remove axis labels and title
#'     ggplot2::theme_classic() +
#'     ggplot2::theme(
#'       legend.position = "none",
#'       axis.text.x = element_blank(),
#'       axis.ticks.x = element_blank(),
#'       # axis.text.y=element_blank(),
#'       axis.ticks.y = element_blank(),
#'       axis.title = element_blank(),
#'       axis.text.y = element_text(color = "black", size = 10),
#'       axis.line = element_line(size = 0.1, colour = "white", linetype = 2)
#'     )
#'
#'   if (!is.null(custom_colors)) {
#'     p <- p + ggplot2::scale_color_manual(values = custom_colors)
#'   }
#'
#'   if (reverse_order) {
#'     # p <- p + ggplot2::scale_y_continous(limits=rev)
#'   }
#'   # when I have not defined a x_limit value, than I need here to switch the
#'   # padding of, so that the y limits of the km plot and the risk table are the same
#'   # and I do not have a shift in the breaks
#'   # If I have x_limits defined, than the paddings are the same in the km plot
#'   # and in the risk table
#'   if (x_lim_present == FALSE) {
#'     p <- p + ggplot2::scale_x_continuous(limits = x_lim, breaks = x_breaks,
#'                                          expand = c(0, 0))
#'   } else {
#'     ggplot2::scale_x_continuous(limits = x_lim, breaks = x_breaks)
#'   }
#'   p_data <- ggplot_build(p)
#'
#'   x_breaks1 <- p_data$layout$panel_params[[1]]$x$breaks
#'   x_lim1 <- p_data$layout$panel_params[[1]]$x.range
#'   print(c("Risk table", x_breaks1))
#'   print(c("Risk table", x_lim1))
#'
#'   return(p)
#' }
#'
#'
#' ###########################################################
#' # risk table for multiple groups
#' # this function calls internally the get_risk_table function
#' #############################################################
#'
#' #' This function calculated a risk table for two groups
#' #'
#' #' @description
#' #' this function calls internally the `get_risk_table`
#' #'
#' #'
#' #' @inheritParams km_grouped
#' #' @param times numeric vector indicating the breaking point on which the number at risks
#' #'              and numer of events should be displayed
#' #'
#' #' @return a data frame containing the number at risk and number of events
#' #'         for the presepcified timepoints stratified after groups
#' #' @noRd
#' #' @importFrom dplyr bind_rows mutate_at mutate
#' #' @importFrom tidyr gather
#' #' @importFrom rlang .data
#' #'
#' get_risk_table.groups <- function(times, data, variable, ev_time, event = NULL,
#'                                   weights, show_weighted_n = TRUE) {
#'   levs <- levels(data[, variable])
#'   out <- vector(mode = "list", length = length(levs))
#'   for (i in seq_len(length(levs))) {
#'     dat_temp <- data[data[, variable] == levs[i], ]
#'
#'     if (show_weighted_n == TRUE) {
#'       weights_i <- weights[data[, variable] == levs[i]]
#'     } else {
#'       weights_i <- NULL
#'     }
#'
#'
#'
#'     out_i <- get_risk_table.all(
#'       times = times, data = dat_temp, ev_time = ev_time,
#'       event = event, weights = weights_i, show_weighted_n = show_weighted_n
#'     )
#'
#'     out_i$group <- levs[i]
#'     out[[i]] <- out_i
#'   }
#'   out <- bind_rows(out)
#'   out$group <- factor(out$group, levels = levs)
#'
#'
#'
#'   # doing some customazion for easier plotting of the risk table
#'   out <- out %>%
#'     mutate_at(c("n_risk", "n_event"), function(x) round(x, 1))
#'   #out <- gather(out, .data$measure, .data$value, c(.data$n_risk, .data$n_event), factor_key = TRUE)
#'   out<- gather(out, measure, value, c(n_risk, n_event), factor_key=TRUE)
#'   out <- out %>%
#'     mutate(
#'       group = (as.numeric(.data$group) - 1),
#'       # I want it the other way round: the reference level of the
#'       # group should be 1, because this is printed above the other level
#'       group = ifelse(.data$group == 0, 1, 0),
#'       group = ifelse(.data$measure == "n_risk", .data$group + 0.3, .data$group),
#'       group1 = factor(.data$group, levels = c(0, 0.3, 1, 1.3))
#'     )
#'
#'   return(out)
#' }
#'
#'
#' #' This function calculates the risk table for one group
#' #'
#' #' @inheritParams get_risk_table.groups
#' #'
#' #' @return a data frame containing the number at risk and number of events
#' #'         for the presepcified timepoints
#' #'
#' #' @noRd
#' #'
#' get_risk_table.all <- function(times, data, ev_time, event = NULL,
#'                                weights, show_weighted_n = TRUE) {
#'   # set weights to all 1 if show_weighted_n = FALSE
#'   if (show_weighted_n == FALSE) {
#'     weights <- rep(1, nrow(data))
#'   } else {
#'     weights <- weights
#'   }
#'
#'   # number at risk
#'
#'   n_risk <- vapply(times, function(x) {
#'     sum(weights[data[, ev_time] >= x])
#'   },
#'   FUN.VALUE = numeric(1)
#'   )
#'
#'   n_event <- vapply(times,
#'                     FUN = function(x) {
#'                       sum(weights[data[, ev_time] < x &
#'                                     data[, event] == 1])
#'                     },
#'                     FUN.VALUE = numeric(1)
#'   )
#'   # cumulative number of right-censored observations
#'
#'   out <- data.frame(time = times, n_risk = n_risk, n_event = n_event)
#'
#'   return(out)
#' }
#'
