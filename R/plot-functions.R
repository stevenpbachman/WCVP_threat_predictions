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
    p <- p + geom_linerange(mapping=aes(xmin=0, xmax=.estimate, group=.data[[fill_var]]), colour="grey80",
                          position=position_dodge(width=0.75), size=1)
    p <- p + geom_point(mapping=aes(colour=.data[[fill_var]]), position=position_dodge(width=0.75), size=3)
    #p <- p + geom_col(mapping=aes(fill=.data[[fill_var]]), position=position_dodge(width=0.75), width=0.75)
  } else {
    p <- p + geom_point(position=position_dodge())
  }
  
  if ("n" %in% colnames(data)) {
    labels <- 
      data |>
      group_by(.data[[y_var]], n) |>
      summarise(.estimate=max(.estimate, na.rm=TRUE), .groups="drop") |>
      mutate(n=scales::label_comma()(n))
      
    p <- p + geom_text(data=labels, mapping=aes(label=n), hjust=-0.1, size=3, colour="grey50")
  }
    
  p <- p +
    scale_x_continuous(limits=c(0, 1), expand=expansion(add=c(0, 0.1))) +
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
plot_threat <- function(data, y_var, draws=NULL) {
  is_samples <- ! is.null(draws)
  is_conformal <- sum(str_detect(colnames(data), ".set_")) > 0
  
  if (is_samples) {
    dist <-
      draws |>
      left_join(
        data |> select(plant_name_id, set, {{ y_var }}),
        by=c("plant_name_id", "set")
      ) |>
      group_by({{ y_var }}, .draw) |>
      summarise(threatened=mean(threatened),
                .groups="drop") |>
      mutate({{ y_var }} := reorder({{ y_var }}, threatened))
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
  
  if (is_conformal) {
    props <-
      props |>
      summarise(threatened=mean(threatened),
                .lower=mean(.set_threatened),
                .upper=mean(ifelse(`.set_not threatened`, FALSE, .set_threatened)),
                status="predicted",
                .groups="drop")
  } else if (is_samples) { 
    props <-
      props |>
      median_qi(threatened) |>
      mutate(status="predicted")
  } else {
    props <-
      props |>
      summarise(threatened=mean(threatened),
                status="predicted",
                .groups="drop")
  }
  
  props <- mutate(props, {{ y_var }} := reorder({{ y_var }}, threatened))
  y_levels <- levels(pull(props, {{ y_var }}))
  obs <- mutate(obs, {{ y_var }} := factor({{ y_var }}, levels=y_levels))
    
  p <- ggplot(data=props, mapping=aes(y={{ y_var }}, x=threatened)) +
    geom_col(data=obs, mapping=aes(fill=status))
  
  if (is_conformal | is_samples) {
    p <- p + geom_pointinterval(mapping=aes(xmin=.lower, xmax=.upper, colour=threatened),
                                show.legend=FALSE)
  } else {
    p <- p + geom_line(mapping=aes(group={{ y_var }}), size=2, colour="grey80")
  }
  
  p +
    scale_x_continuous(limits=c(0, 1), labels=scales::label_percent()) +
    scale_fill_manual(values=c("observed"="grey80", "predicted"="black"), name="") +
    labs(x="Threatened plant species", y="")
}


#' Plot a choropleth map, where the WGSRPD level 3 regions are coloured by some quantity.
#'
#' @param data A data frame with some quantity for each of the WGSRPD3 level 3 regions.
#' @param fill_var The column to colour the regions by.
#'
plot_map <- function(data, fill_var, .proj=c("moll", "igh"), .points=FALSE, .shapes=NULL) {
  .proj <- match.arg(.proj)
  
  if (is.null(.shapes)) {
    filled_regions <- 
      rWCVPdata::wgsprd3 |>
      left_join(data, by=c("LEVEL3_COD"="area_code_l3"))
  } else {
    filled_regions <-
      .shapes |>
      left_join(data)
  }
  
  
  if (.proj == "goode") {
    filled_regions <- 
      filled_regions |>
      st_wrap_dateline(options=c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) |>
      st_transform(st_crs("+proj=igh"))
  } else {
    filled_regions <- 
      filled_regions |>
      st_wrap_dateline(options=c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) |>
      st_transform(st_crs("+proj=moll"))
  }
    
  
  p <- ggplot() +
    geom_sf(data=filled_regions, mapping=aes(fill={{ fill_var }}), colour="grey90", size=0.5/.pt)
  
  if (.points) {
    points_sf <- 
      filled_regions |>
      st_make_valid() |>
      mutate(island=sapply(st_intersects(filled_regions, filled_regions), function(x) (length(x) == 1)),
             area=st_area(geometry)) |>
      filter(island, area <= units::set_units(2e10, "m^2")) |>
      st_centroid()
    
    # ALU is split at dateline, so need to shift centroid
    points_sf[points_sf$LEVEL3_COD == "ALU",]$geometry <- points_sf[points_sf$LEVEL3_COD == "ALU",]$geometry - c(60, 0)
    
    p <-
      p +
      geom_sf(data=points_sf, mapping=aes(colour={{ fill_var }}), 
              show.legend=FALSE, size=.pt / 6)
  }
  
  if (.proj == "goode") {
    p <- 
      p +
      geom_sf(data=goode_outline(mask=TRUE), colour=NA, fill="#ffffff", size=0.5/.pt) +
      geom_sf(data=goode_outline(), fill=NA, colour="grey50", size=0.5) +
      coord_sf(crs="+proj=igh")
  } else {
    p <- 
      p + 
      geom_sf(data=moll_outline(), fill=NA, colour="grey50", size=0.5) +
      coord_sf(crs="+proj=moll")
  }
  p <- 
    p +
    scico::scale_fill_scico(
      palette="batlowK",
      na.value="grey80"
    ) +
    guides(fill=
      guide_colourbar(
        direction="horizontal", 
        title.position="top", 
        label.position="bottom",
        label.hjust=0.5,
        barwidth=15,
        barheight=1
      )
    ) 
  if (.points) {
    p <- 
      p +
      scico::scale_color_scico(
        palette="batlowK",
        na.value="grey80"
      )
  }
  p +
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

moll_outline <- function() {
  lats <- c(
    90:-90,
    -90:90,
    90
  )
  
  longs <- c(
    rep(180, 181),
    rep(-180, 181),
    180
  )
  outline <-
    list(cbind(longs, lats)) %>%
    sf::st_polygon() %>%
    sf::st_sfc(
      crs="+proj=longlat +ellps=WGS84 +no_defs"
    )
  
  outline
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
