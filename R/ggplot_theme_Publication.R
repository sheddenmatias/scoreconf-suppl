# Based on https://rpubs.com/koundy/71792

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "plain",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "plain",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            #legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="plain")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
# Shared scales for the four statistics. The colours are the Okabe-Ito
# barrier-free palette, which separates under the common forms of colour
# blindness and in greyscale; linetype repeats the identity so it never rests on
# colour alone. 01, 05 and 07 currently define their own colours, in which the
# green for SCR_C and the red for LRT are 4.1 apart under deuteranopia; adopting
# these scales makes the figures consistent and legible.
stat_colours <- c(WLD = "#0072B2", LRT = "#D55E00",
                  SCR_C = "#009E73", SCR_U = "#CC79A7")
stat_linetypes <- c(WLD = "dotdash", LRT = "solid",
                    SCR_C = "dotted", SCR_U = "dashed")
scale_colour_stat <- function(...)
  ggplot2::scale_colour_manual(values = stat_colours, ...)
scale_linetype_stat <- function(...)
  ggplot2::scale_linetype_manual(values = stat_linetypes, ...)

# Display names. "Score, extended" is the interval the paper proposes: the same
# statistic as "Score", with the nuisance parameters maximized over the extended
# parameter set.
# Wald and the likelihood ratio are restricted too, so only the extended
# statistic carries a qualifier.
stat_labels <- c(WLD = "Wald", LRT = "Likelihood ratio",
                 SCR_C = "Score", SCR_U = "Score, extended")
scale_colour_stat <- function(...)
  ggplot2::scale_colour_manual(values = stat_colours, labels = stat_labels, ...)
scale_linetype_stat <- function(...)
  ggplot2::scale_linetype_manual(values = stat_linetypes, labels = stat_labels, ...)
