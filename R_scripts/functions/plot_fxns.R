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
plot_RE_ts_ggplot <- function(data, x, y, ylim = c(NULL, NULL), lwr_1, upr_1,
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

plot_RE_ts_base <- function(data, par_name, ylim) {
  
  # Filter for a particular parameter name
  data<- data[data$par_name == par_name, ]
  
  # Empty plot here
  plot(1, type="n", ylim = ylim, xlim = c(min(data$year), max(data$year)),
       xlab = "Year", ylab = "Relative Error", main = unique(data$par_name))
  
  # Percentile intervals here
  # 100% percentiles
  polygon(x = c(data$year, rev(data$year)),  y = c(data$lwr_100, data$upr_100), 
          col =  adjustcolor("blue", alpha.f = 0.2), border = NA,
          lty = 1, lwd = 3) 
  
  # 95% percentiles
  polygon(x = c(data$year, rev(data$year)),  y = c(data$lwr_95, data$upr_95), 
          col =  adjustcolor("blue", alpha.f = 0.4), border = NA,
          lty = 1, lwd = 3) 
  
  # 75% percentiles
  polygon(x = c(data$year, rev(data$year)),  y = c(data$lwr_75, data$upr_75), 
          col =  adjustcolor("blue", alpha.f = 0.6), border = NA,
          lty = 1, lwd = 3) 
  
  # Add median points
  points(data$year, y = data$median, ylim = ylim, col = "black",
         lwd = 1, pch = 21, cex = 2.5, bg = "white")
  abline(h = 0, col = "black", lty = 2, lwd= 3) # 0 bias line
  
  
}

