# Purpose: To plot relative error metrics
# Creator: Matthew LH. Cheng
# Date 1/12/22


#' Title To plot relative error time series
#'
#' @param data datframe
#' @param x character of x axis variable
#' @param y character of y axis variable
#' @param ylim vector of y limits
#' @param lwr_1 character of lwr CI 1
#' @param upr_1 character of upr CI 1
#' @param lwr_2 character of lwr CI 2
#' @param upr2 character of upr CI 2
#' @param facet_name character name of variable we want to facet by
#'
#' @return
#' @export
#'
#' @examples
plot_RE_ts <- function(data, x, y, ylim = c(NULL, NULL), lwr_1, upr_1,
                    lwr_2, upr_2, facet_name) {
  
  require(rlang)
  
  plot_RE <- ggplot(data, aes(x = {{x}}, y = {{y}}) ) +
    geom_ribbon(aes(ymin = {{lwr_1}}, ymax = {{upr_1}}), alpha = 0.5, fill = "grey4") +
    geom_ribbon(aes(ymin = {{lwr_2}}, ymax = {{upr_2}}), alpha = 0.3, fill = "grey2") +
    geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.8, alpha = 0.85) +
    geom_hline(aes(yintercept = 0), col = "black", lty = 2, size = 1, alpha = 1) +
    facet_wrap(enquo(facet_name), scales = "free", ncol = 2) +
    coord_cartesian(ylim = c(ylim)) +
    labs(x = "Year", y = "Relative Error") +
    theme_bw() +
    theme(strip.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13, color = "black"))
    
    return(plot_RE)
}

