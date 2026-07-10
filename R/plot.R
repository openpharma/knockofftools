
#' Heatmap of multiple variable selections ordered by importance
#'
#' @param x data.frame of variable selections from multiple knockoffs
#' (each entry is either 1 if variable is selected and 0 otherwise).
#' Columns correspond to different knockoffs and rows correspond to the
#' underlying variables. row.names(x) records the variable names.
#'
#' @param ... Additional arguments passed to other plot methods (currently ignored).
#'
#' @param nbcocluster bivariate vector c(number of variable clusters, number of selection clusters).
#' The former number must be specified less than nrow(x) and the latter must be less than ncol(x).
#'
#' @details To help visualize most important variables we perform clustering both selections and variables.
#'
#' @return plot of heatmap
#'
#' @method plot variable.selections
#' @export
#'
#' @examples
#' library(knockofftools)
#' set.seed(1)
#'
#' # Simulate 8 Gaussian covariate predictors and 2 binary factors:
#' X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)
#'
#' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
#' lp <- generate_lp(X, p_nn = 5, a=1)
#'
#' # Simulate response:
#' y <- lp + rnorm(100)
#'
#' # Calculate knockoff statistics:
#' W <- knockoff.statistics(y=y, X=X, type="regression", M=5)
#'
#' S <- variable.selections(W, error.type = "pfer", level = 1)
#'
#' # plot heatmap of knockoff selections:
#' plot(S)
plot.variable.selections <- function(x, ..., nbcocluster=c(7,7)) {

  S <- x

  # Check x is of class variable.selections to use components 'selected' and 'stable.variables'
  if (!inherits(x, "variable.selections")) {
    stop("Input x must be of class \'variable.selections\'. Please see ?variable.selections.")
  }

  stable.vars <- S$stable.variables
  S <- S$selected
  if (ncol(S) < 2 || nrow(S) < 2) {
    stop("Input x must have at least 2 knockoff repetitions and 2 variables to plot the clustered heatmap.")
  }

  # Basic checks on nbcocluster ensuring clustering doesn't fail
  stopifnot(
    is.numeric(nbcocluster),
    length(nbcocluster) == 2,
    all(nbcocluster > 0)
  )
  # Make sure user doesn't specify more clusters than data to support:
  nbcocluster <- pmin(nbcocluster, dim(S))

  axis.text.y.color <- list(
    selected     = "firebrick2",
    not_selected = "black"
  )

  selections <- data.frame(draw = factor(rep(1:ncol(S),each=nrow(S))),
                           variable = factor(rownames(S)),
                           selected = as.numeric(as.matrix(S)))

  sel.mat <- matrix(selections$selected,nrow=nrow(S))
  hclust.row <- hclust(dist(sel.mat, method="binary"), method="ward.D")
  hclust.col <- hclust(dist(t(sel.mat), method="binary"), method="ward.D")

  selections$varclass = factor(cutree(hclust.row, k=nbcocluster[1]), ordered=TRUE)
  selections$drawclass = factor(rep(cutree(hclust.col, k=nbcocluster[2]), each=nrow(S)), ordered=TRUE)

  # Calculate means per draw cluster (to order heatmap)
  meta.draw <- selections %>%
    dplyr::group_by(drawclass) %>%
    dplyr::summarise(selected = mean(selected)) %>%
    dplyr::arrange(dplyr::desc(selected))

  # Calculate means per variable block (to order heatmap)
  meta.var <- selections %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(selected = mean(selected))%>%
    dplyr::arrange(selected)

  # Order selections according to block cluster means:
  selections <- selections %>%
    dplyr::mutate(drawclass = factor(drawclass, levels = as.character(meta.draw$drawclass)),
                  variable = factor(variable, levels = as.character(meta.var$variable))) %>%
    dplyr::arrange(drawclass) %>%
    dplyr::mutate(draw = factor(draw, levels = unique(draw))) %>%
    dplyr::mutate(
      draw_plot = match(draw, levels(draw)),
      selected = factor(selected, levels = 0:1, labels = c("no", "yes")),
      variable_colored_label = ifelse(
        variable %in% stable.vars,
        paste0("<b><span style='color:", axis.text.y.color$selected,     ";font-weight:bold;'>", variable, "</span></b>"),
        paste0(   "<span style='color:", axis.text.y.color$not_selected, ";font-weight:bold;'>", variable, "</span>")
      ) # forcats::fct_reorder(as.numeric(variable))
      ) %>%
    dplyr::arrange(variable) %>%
    dplyr::mutate(
      variable_colored_label = factor(
        variable_colored_label,
        levels = unique(variable_colored_label)
        )
    )

  if (utils::packageVersion("ggplot2") >= "3.5.0") {

    ggplot2::ggplot(
      data = selections,
      mapping = ggplot2::aes(x = draw_plot, y = as.numeric(variable_colored_label))
    ) +
      ggplot2::geom_tile(ggplot2::aes(fill = selected)) +
      ggplot2::scale_fill_manual(
        values = c("no" = "#132B43", "yes" = "#56B1F7")
      ) +
      ggplot2::scale_x_continuous(
        expand = 0,
        breaks = seq(1, ncol(S), by = 1),
        minor_breaks = seq(1.5, ncol(S) - 0.5, by = 1)
      ) +
      ggplot2::scale_y_continuous(
        expand = 0,
        labels = levels(selections$variable_colored_label),
        breaks = seq(1, nrow(S), by = 1),
        minor_breaks = seq(1.5, nrow(S) - 0.5, by = 1)
      ) +
      ggplot2::labs(
        x = "knockoff repetition (reordered)",
        y = "variable (stable selections colored in red)"
      ) +
      ggplot2::guides(
        x = ggplot2::guide_axis(minor.ticks = TRUE),
        y = ggplot2::guide_axis(minor.ticks = TRUE)
      ) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_blank(),
        axis.text.y  = ggtext::element_markdown(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.minor.ticks.x.bottom = ggplot2::element_line(linewidth = ggplot2::rel(1)),
        axis.minor.ticks.y.left = ggplot2::element_line(linewidth = ggplot2::rel(1)),
      )

  } else {

    message(
      "You are currently using ggplot2 version ",
      utils::packageVersion("ggplot2"), ".\n",
      "Consider upgrading to at least version 3.5.0. for an updated plot version."
    )

    lbl_var_marked <- levels(selections$variable) %>%
      paste0(ifelse(levels(selections$variable) %in% stable.vars, "*", ""))

    selections %>%
      dplyr::mutate(variable = factor(variable, labels = lbl_var_marked)) %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = draw, y = variable)) +
      ggplot2::geom_tile(ggplot2::aes(fill = selected)) +
      ggplot2::xlab("knockoff repetition (reordered)") +
      ggplot2::ylab("variable (stable selections marked with *)") +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::scale_fill_manual(values = c("no" = "#132B43","yes" = "#56B1F7"))

  }

}


suppress_warning_vectorized_input_element_text <- function(expr) {

  # Vectorized input to `element_text()` is not officially supported.
  # ℹ Results may be unexpected or may change in future versions of ggplot2.
  pattern <- paste0(
    "Vectorized.input.*element_text.*not.officially.supported.",
    "Results.*unexpected.*change.*ggplot2."
  )

  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(pattern, conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
