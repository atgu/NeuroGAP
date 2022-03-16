## some helper functions

# functions to convert between spherical lat-long and 3D XYZ coordinates
sphere_to_coord_mtrx <- function(d) {
  if (length(dim(d) > 1)) {
    d1 <- d[, 1] + 90
    d2 <- d[, 2]
  } else {
    d1 <- d[1] + 90
    d2 <- d[2]
  }
    
  return(cbind(
    sinpi(d1/180) * cospi(d2/180),
    sinpi(d1/180) * sinpi(d2/180),
    cospi(d1/180)
  ))
}

cartesian_to_spherical_mtrx <- function(d) {
  if (length(dim(d)) == 2) {
    d1 <- d[, 1]
    d2 <- d[, 2]
    d3 <- d[, 3]
    dsum2 <- rowSums(d ** 2)
  } else {
    d1 <- d[1]
    d2 <- d[2]
    d3 <- d[3]
    dsum2 <- sum(d ** 2)
  }
  return(cbind(
    acos(d3 / sqrt(dsum2)) / pi * 180 - 90,
    atan2(d2, d1) / pi * 180
  ))
}


# pca used for procrustes: pca_filtered <- f(pca, row selection, languages used)
procrustes_pca_fxn <- function(pca, rows, langs_used) {
  # Each individual is associated with some (possibly empty) set of
  #  languages, and each language has a position in PC-space. This
  #  fxn finds the mean in PC-space of languages associated with each
  #  individual.
  
  # pca :  matrix, containing the location in PCA-space for each language
  #                based on the linguistic data (phoneme presence-absence)
  # rows : list, with each entry containing the 
  # langs_used :  vector, containing the codes for languages which are used here
  
  rows_used <- c()
  return_pca <- matrix(NA, nrow = 0, ncol = ncol(pca))
  for (i in 1:length(rows)) {
    langs <- rows[[i]]
    
    # langs: set of available languages
    # lang_codes_upd$Code[which(! is.na(lang_codes_upd$lang_family))] # for classification PCA
    # langs_to_use  
    if (! any(langs_used %in% langs)) {
      next
    }
    
    rows_used <- c(rows_used, i)
    
    if (length(which(langs_used %in% langs)) == 1) {
      return_pca <- rbind(return_pca,
                          pca[which(langs_used %in% langs),])
    } else {
      
      pca_col <- colMeans(
        pca[which(langs_used %in% langs), ]
      )
      return_pca <- rbind(
        return_pca,
        pca_col
      )
    }
  }
  # rownames used to store the indiv's with data
  rownames(return_pca) <- rows_used
  
  return(return_pca)
}
 
# pca used for procrustes <- f(phoneme_mtrx, row selection, languages used)
procrustes_phoneme_pca_fxn <- function(phoible_phoneme_mtrx, rows, langs_used,
                                       rows_gen2=NA, rows_gen3=NA,
                                       indiv_subset=1:length(rows)) {
  require(ggplot2)
  require(vegan)
  require(maps)

  indiv_phoneme_matrix <- matrix(NA, nrow = 0, ncol = ncol(phoible_phoneme_mtrx))
  
  rows_used <- c()
  for (i in indiv_subset) {
    langs <- rows[[i]]
    
    if (! all(is.na(rows_gen2))) {
        langs_2 <- rows_gen2[[i]]
    } else { langs_2 <- c() }
    
    if (! all(is.na(rows_gen3))) {
      langs_3 <- rows_gen3[[i]]
    } else { langs_3 <- c() }
    
    langs_all <- union(langs, union(langs_2, langs_3))
    
    if (! any(langs_used %in% langs_all)) {
      next
    } else {
      rows_used <- c(rows_used, i)
    }
    
    if (sum(langs_used %in% langs_all) == 1) {
      indiv_phoneme_matrix <- rbind(
        indiv_phoneme_matrix,
        phoible_phoneme_mtrx[
          which(langs_used %in% langs_all),
        ])
    } else {
      # weighted sum of indiv.'s, parents', and grandparents' languages
      mtrx_col <- matrix(0, nrow = 1, ncol = ncol(phoible_phoneme_mtrx))

      multiplier <- 1
      mean_divisor <- 0
      for (langs_set in list(langs, langs_2, langs_3)) {
        langs_set <- langs_set[which(langs_set %in% langs_used)]
        
        for (lang in langs_set) {
          mtrx_col <- mtrx_col +
            multiplier *
            phoible_phoneme_mtrx[which(langs_used %in% lang), ]
          
          mean_divisor <- mean_divisor + multiplier
        }
        
        multiplier <- multiplier / 2
      }
      
      mtrx_col <- mtrx_col / mean_divisor
      
      indiv_phoneme_matrix <- rbind(
        indiv_phoneme_matrix,
        mtrx_col
      )
    }
  }
  
  return_PCA <- prcomp(as.matrix(dist(
    indiv_phoneme_matrix,
    method = "canberra", diag = T, upper = T
    )))$x
  
  # rownames used to store the indiv's with data
  rownames(return_PCA) <- as.character(rows_used)
  
  return(return_PCA)
} 

# calculate lat/long locations for individuals based on a weighted avg
#  of the languages spoken by themselves and their family
lat_long_calc_weighted <- function() {
  
  longitudes <- c()
  latitudes <- c()
  
  for (indiv in 1:nrow(neuroGAP_lang)) {
    
    indiv_value <- rep(0, 3)
    
    for (i in 0:2) {
      if (i==0) {langs <- neuroGAP_lang$lang_self_all[[indiv]]}
      if (i==1) {langs <- neuroGAP_lang$parents[[indiv]]}
      if (i==2) {langs <- neuroGAP_lang$grandparents[[indiv]]}
      
      langs <- langs[langs %in% langs_to_use]
      if (length(langs) == 0) next
      
      for (lang in langs) {
        # mean_divisor <- mean_divisor + 0.5 ** i
        
        indiv_value <- indiv_value + (0.5 ** i) *
          sphere_to_coord_mtrx(lang_codes_upd[
            which(lang_codes_upd$Code == lang),
            c("longitude", "latitude")])
      }
    }
    
    indiv_long_lat <- cartesian_to_spherical_mtrx(indiv_value)
    longitudes <- c(longitudes, indiv_long_lat[1])
    latitudes <- c(latitudes, indiv_long_lat[2])
  }
  
  return(list(latitudes = latitudes, longitudes = longitudes))
}

# function to run a procrustes analysis, plot the results,
#  and print some potentially useful output
procrust_fxn_short <- function(data1, data2, geographic = T,
                               nperm = 19999, title = NA,
                               pdf_title = NA, EAfr = F,
                               just_return_t = F) {
  require(vegan)
  require(ggplot2)
  require(maps)
  
  rn_d1 <- rownames(data1)
  if (typeof(data1) == "list") {
    data1 <- as.matrix(data1)
  }
  
  if (geographic) {
    # convert longitude/latitude to points on a unit sphere, 
    #   so that euclidean distances are meaningful when
    #   running procrustes.
    data1 <- sphere_to_coord_mtrx(data1)
  }

  # run procrustes
  p_test <- protest(data1, data2, permutations = how(nperm = nperm))
  
  if (just_return_t) {
    bootstrapped_proc <- sapply(
      1:500,
      function (d) {
        samples <- sample(
          1:nrow(data1),
          size = nrow(data1),
          replace = T)
        return(protest(data1[samples, ],
          data2[samples, ],
          permutations = 9)$t0)
        })
    
    print(quantile(bootstrapped_proc,
                   probs = c(0.025, 0.5, 0.975)))
    
    return(list(t0 = p_test$t0, t0vec = bootstrapped_proc))
  }

  
  if (!is.na(pdf_title)) {
    p <- procrustes(data1, data2)

    if (geographic) {
      
      if (EAfr) {
        limits = list(xlim = c(27, 47), ylim = c(-7, 15))
      } else {
        limits = list(xlim = c(15, 50), ylim = c(-35, 20))
      }
      
      # original points, for plotting on a map
      orig_points <- cartesian_to_spherical_mtrx(
        p$X + t(matrix(rep(p$xmean, nrow(p$X)), ncol = nrow(p$X)))
      )
      
      # procrustes-tranformed points, for plotting on a map
      proc_points <- cartesian_to_spherical_mtrx(
        p$Yrot + t(matrix(rep(p$xmean, nrow(p$X)), ncol = nrow(p$X)))
      )
      
    } else {
      orig_points <- p$X + t(matrix(rep(p$xmean, nrow(p$X)), ncol = nrow(p$X)))
      proc_points <- p$Yrot + t(matrix(rep(p$xmean, nrow(p$X)), ncol = nrow(p$X)))
    }
    
    # extract three-letter codes that represent study of origin
    point_qual <- factor(c(sub(
      "[0-9]*$", "",
      neuroGAP_data$IID[as.integer(rn_d1)]
    )))
    color_list <- c(
      "AAP"=rgb(0, 0, 0),
      "CTP"=rgb(.9, .4, .5),
      "KWP"=rgb(.2, .6, .9),
      "MAP"=rgb(1, 0, 0),
      "MOP"=rgb(0, .9, 0))
    shape_list <- c("AAP"=8, "CTP"=10, "KWP"=7,
                    "MAP"=12, "MOP"=9)
    color_list <- color_list[names(color_list) %in% point_qual]
    
    shape_list <- shape_list[names(shape_list) %in% point_qual]    
    
    plotting_list <- list(x=c(), y=c(), numind=c(), point_qual=c())

    for (upoint1 in unique(proc_points[, 1])) {
      rows_up1 <- which(proc_points[, 1] == upoint1)

      for (upoint2 in unique(proc_points[rows_up1, 2])) {
        rows_up <- intersect(rows_up1, which(proc_points[, 2] == upoint2))

        for (qual in unique(point_qual[rows_up])) {
          rows_used <- rows_up[point_qual[rows_up] == qual]
          plotting_list$x <- append(plotting_list$x, orig_points[rows_used[1], 1])
          plotting_list$y <- append(plotting_list$y, orig_points[rows_used[1], 2])
          plotting_list$xend <- append(plotting_list$xend, upoint1)
          plotting_list$yend <- append(plotting_list$yend, upoint2)
          plotting_list$numind <- append(plotting_list$numind, length(rows_used))
          plotting_list$point_qual <- append(plotting_list$point_qual, qual)
        }
      }
    }
    plotting_df <- as.data.frame(plotting_list)
    plotting_df <- plotting_df[sample(1:nrow(plotting_df)),]
    
    if (geographic) {
      world <- map_data("world")
      ggp <- ggplot(world, mapping=aes(long, lat, group = group)) +
        geom_polygon(fill="#f0e8d8", color='#00000010') +
        coord_fixed(xlim = limits$xlim, ylim = limits$ylim) +
        theme_void() + theme(
            panel.background = element_rect(fill='#A0C0F0',
            linetype=0)
        )
    } else {
      ggp <- ggplot() + theme_classic() +
        theme(
          rect = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(color = '#505050', size = .5),
          panel.grid = element_blank(),
          # panel.grid.minor.y = element_line(color = '#CCCCCC', size = .1),
          # panel.grid.major.y = element_line(color = '#CCCCCC', size = .2),
          plot.title = element_text(size = 8),
          axis.ticks.length = unit(-0.05, "in"),
          axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
          axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
          axis.ticks.x = element_blank(),
          aspect.ratio = 0.8,
          legend.background = element_rect(color = "white", fill = "white"),
          text = element_text(size = 20)
        )
    }

    # setup plotting theme
    ggp <- ggp + ggtitle(label=title)

    # add mean points representing the center-points of geographic sites
    mean_aggr <- aggregate(x=orig_points, by = list(point_qual), mean)
    colnames(mean_aggr)[1:3] <- c("site", "x", "y")

    ggp <- ggp + geom_point(
      data = mean_aggr,
      mapping = aes(x=x, y=y, color=site, shape=site),
      alpha = 1, inherit.aes = FALSE, size = 8) +
      scale_color_manual(values = color_list) +
      scale_shape_manual(values = shape_list) +
      guides(shape = "none")

    # add procr-rotated points
    range_max <- min(8, log(2 + 0.5 * max(plotting_df$numind)))
    ggp <- ggp + geom_point(
        plotting_df,
        mapping = aes(x=xend, y=yend,
                      color=point_qual, size=numind),
        inherit.aes = FALSE, # color = rgb(1, 1, 1, 0), 
        alpha = 0.6, shape = 16) +
      scale_color_manual(values = color_list) +
      scale_size_continuous(range = c(4, range_max)) +
      labs(color = "Site", size = "Number of individuals")

    if (max(plotting_list$numind) == 1) {
      ggp <- ggp + guides(size = "none")
    }

    # # # Draw arrows
    # # ggp <- ggp + geom_segment(
    # #   data = plotting_df,
    # #   mapping = aes(x=x, y=y, xend=xend, yend=yend,
    # #                 size=numind, color=point_qual),
    # #   arrow = arrow(), inherit.aes = FALSE, alpha = 0.5,
    # #   lineend = "butt", linejoin = "mitre"
    # # )

    ggsave(ggp, filename = pdf_title, device = cairo_pdf, 
           width = 8, height = 10, units = "in")   
    
    
  }
  print(title)
  print(p_test)
}

