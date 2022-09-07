#' Remove the x-axis and all related items for sharing the x-axis in a plot
#' patchwork.
#' 
remove_xaxis <- function() {
  theme(
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank()
  )
}

#' Plot test set performance across subgroups as a bar chart.
#' 
#' @param data A data frame of test set performance.
#' @param y_var The name of a column to plot on the y-axis (subgroups).
#' @param fill_var The name of a column to colour the bars by.
#' @param metric The name of a metric to plot. If `NULL`, will plot all.
#' 
plot_performance_bars <- function(data, y_var, fill_var=NULL, metric=NULL) {
  
  if (! is.null(metric)) {
    data <- filter(data, .metric == metric)
  }
  
  data[[y_var]] <- reorder(data[[y_var]], data$.estimate)

  p <- ggplot(data, mapping=aes(x=.estimate, y=.data[[y_var]]))
  
  if (! is.null(fill_var)) {
    p <- p + geom_col(mapping=aes(fill=.data[[fill_var]]), position=position_dodge())
  } else {
    p <- p + geom_col(position=position_dodge())
  }
  
  if ("n" %in% colnames(data)) {
    p <- p + geom_text(mapping=aes(label=n, group=.data[[fill_var]]), hjust=-0.1, 
                       position=position_dodge(width=1))
  }
    
  p <- p +
    scale_x_continuous(limits=c(0, 1), expand=expansion(add=0.1)) +
    labs(x="", y="")
  
  if (is.null(metric)) {
    p <- p + facet_wrap(~.metric)
  }
  
  p
}

#' Plot the predicted and observed proportion of threatened species across groups.
#' 
#' If the predictions are probabilistic (i.e. a `.draw` column is present), a summary
#' of their distribution will be plotted. Otherwise, it will be a joined dot-plot.
#' 
#' @param data A data frame with columns for the predicted threat status of a species
#'   (`threatened` or `.pred_class`) and the observed threat status (`obs`).
#' @param y_var The column to group the predictions by, which will form the items on
#'  the y-axis of the plot.
#'
plot_threat <- function(data, y_var) {
  is_samples <- ".draw" %in% colnames(data)
  
  if (is_samples) {
    dist <-
      data |>
      group_by({{ y_var }}, .draw) |>
      summarise(threatened=mean(threatened),
                .groups="drop") |>
      mutate({{ y_var }} := reorder({{ y_var }}, threatened))
    data <-
      data |>
      filter(.draw == 1)
  } else {
    data <- 
      data |>
      mutate(threatened=ifelse(.pred_class == "threatened", 1, 0))
  }
  
  obs <-
    data |>
    group_by({{ y_var }}) |>
    summarise(threatened=mean(obs == "threatened", na.rm=TRUE),
              status="observed",
              .groups="drop") |>
    replace_na(list(threatened=0))
  
  if (is_samples) {
    props <-
      dist |>
      group_by({{ y_var }})
  } else {
    props <-
      data |>
      group_by({{ y_var }})
  }
  
  props <-
    props |>
    summarise(threatened=mean(threatened),
              status="predicted",
              .groups="drop") |>
    bind_rows(obs) |>
    pivot_wider(names_from=status, values_from=threatened) |>
    mutate({{ y_var }} := reorder({{ y_var }}, predicted)) |>
    pivot_longer(c(observed, predicted), names_to="status", values_to="threatened")
  
  p <- ggplot(data=props, mapping=aes(y={{ y_var }}, x=threatened))
  
  if (is_samples) {
    p <- p + stat_gradientinterval(data=dist, point_interval="mean_qi", fill_type="segments")
  } else {
    p <- p + geom_line(mapping=aes(group={{ y_var }}), size=2, colour="grey80")
  }
  
  p +
    geom_point(mapping=aes(colour=status), size=3) +
    scale_x_continuous(limits=c(0, 1), labels=scales::label_percent()) +
    scale_colour_manual(values=c("observed"="red", "predicted"="black"), name="") +
    labs(x="Threatened plant species", y="")
}

#' Plot a choropleth map, where the WGSRPD level 3 regions are coloured by some quantity.
#'
#' @param data A data frame with some quantity for each of the WGSRPD3 level 3 regions.
#' @param fill_var The column to colour the regions by.
#'
plot_map <- function(data, fill_var) {
  filled_regions <- 
    rWCVPdata::wgsprd3 |>
    left_join(data, by=c("LEVEL3_COD"="area_code_l3")) |>
    st_wrap_dateline(options=c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) |>
    st_transform(st_crs("+proj=igh"))
  
  ggplot() +
    geom_sf(data=filled_regions, mapping=aes(fill={{ fill_var }}), colour="grey90", size=0.5/.pt) +
    geom_sf(data=goode_outline(mask=TRUE), colour=NA, fill="#ffffff", size=0.5/.pt) +
    geom_sf(data=goode_outline(), fill=NA, colour="grey50", size=0.5) +
    scico::scale_fill_scico(
      palette="batlowK",
      na.value="grey80"
    ) +
    coord_sf(crs="+proj=igh") +
    guides(fill=
      guide_colourbar(
        direction="horizontal", 
        title.position="top", 
        label.position="bottom",
        label.hjust=0.5,
        barwidth=15,
        barheight=1
      )
    ) +
    theme_map(legend.position="bottom")
}

#' A (hopefully) nice looking theme for maps, derived from 
#' https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
#'
theme_map <- function(background_colour="#ffffff", grid_colour="#dbdbd9", ...) {
  theme_minimal() +
    theme(
      # remove axes
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      # put in a subtle grid
      panel.grid.major=element_line(colour=grid_colour, size=0.2),
      panel.grid.minor=element_blank(),
      # shade background
      plot.background=element_rect(fill=background_colour, colour=NA),
      panel.background=element_rect(fill=background_colour, colour=NA),
      legend.background=element_rect(fill=background_colour, colour=NA),
      panel.border=element_blank(),
      # set margins
      plot.margin=unit(c(0.5, 0.5, 0.2, 0.5), "cm"),
      panel.spacing=unit(c(-0.1, 0.2, 0.2, 0.2), "cm"),
      # define specific text settings
      legend.title=element_text(size=11),
      legend.text=element_text(size=9, hjust=0),
      plot.title=element_text(size=15, hjust=0),
      plot.subtitle=element_text(size=10, hjust=0.5,
                                 margin=margin(b=-0.1, t=-0.1, l=2, unit="cm"),
                                 debug=FALSE),
      plot.caption=element_text(size=7, hjust=0.5, margin=margin(t=0.2, b=0, unit="cm")),
      ...
    )
}

#' a handy function for making an outline or mask of the goode map projection
#' adapted from https://wilkelab.org/practicalgg/articles/goode.html
#' 
goode_outline <- function(mask=FALSE) {
  lats <- c(
    90:-90,
    -90:0, 0:-90,
    -90:0, 0:-90,
    -90:0, 0:-90,
    -90:90,
    90:0, 0:90,
    90
  )
  
  longs <- c(
    rep(180, 181),
    rep(c(80.01, 79.99), each=91),
    rep(c(-19.99, -20.01), each=91),
    rep(c(-99.99, -100.01), each=91),
    rep(-180, 181),
    rep(c(-40.01, -39.99), each=91),
    180
  )
  
  outline <-
    list(cbind(longs, lats)) %>%
    sf::st_polygon() %>%
    sf::st_sfc(
      crs="+proj=longlat +ellps=WGS84 +no_defs"
    )
  
  if (mask) {
    outline <- sf::st_transform(outline, crs="+proj=igh")
    
    xlim <- sf::st_bbox(outline)[c("xmin", "xmax")]*1.1
    ylim <- sf::st_bbox(outline)[c("ymin", "ymax")]*1.1
    
    rectangle <- list(
      cbind(
        c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
        c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
      )
    )
    
    rectangle <- 
      rectangle %>%
      sf::st_polygon() %>%
      sf::st_sfc(crs="+proj=igh")
    
    outline <- sf::st_difference(rectangle, outline)
  }
  
  outline
}
