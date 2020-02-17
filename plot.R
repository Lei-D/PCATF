library(tidyverse)
library(tidyr)
library(RColorBrewer)

th = theme_minimal() + theme(
  # legend.title = element_blank(),
  legend.position = 'bottom',
  strip.background = element_rect(fill='gray70',color='gray70'),
  text = element_text(size=12,family='Palatino')
)

plot.one = function(n, p, d, mu1, gamma1, 
                    outlier.index = NULL, n.outlier = NULL,
                    SNR = NULL, sigma0 = NULL, group, 
                    V, mu0, gamma0, lambda, seed = NULL, 
                    par = c("SNR", "mu1", "gamma1", "n.outlier")){
  # place holder for the output
  out = data.frame(Parameter = character(),
                   Method = character(),
                   Index = integer(),
                   L2norm = double(),
                   Outlier = logical())
  if (par == "SNR") {
    for (i in SNR) {
      # simulate data
      data = dataModel(n = n, p = p, d = d, 
                       mu1 = mu1, gamma1 = gamma1, 
                       outlier.index = outlier.index, 
                       SNR = i, group = group, 
                       V = V[, 1:d, drop = FALSE], 
                       mu0 = mu0, gamma0 = gamma0, 
                       seed = seed)
      # compute estimated U for various methods
      df = compare.methods(data = data, lambda = lambda, d = d)
      # combine output
      out = rbind(out, 
                  data.frame(Parameter = rep(paste(par, "=", i, sep = ""),  
                                             nrow(df)),
                             Method = df$method,
                             Index = df$index,
                             L2norm = df$l2norm,
                             Outlier = df$outlier))
    } 
  } else if (par == "mu1") {
    for (i in mu1) {
      # simulate data
      data = dataModel(n = n, p = p, d = d, 
                       mu1 = i, gamma1 = gamma1, 
                       outlier.index = outlier.index, 
                       sigma0 = sigma0, group = group, 
                       V = V[, 1:d, drop = FALSE], 
                       mu0 = mu0, gamma0 = gamma0, 
                       seed = seed)
      # compute estimated U for various methods
      df = compare.methods(data = data, lambda = lambda, d = d)
      # combine output
      out = rbind(out, 
                  data.frame(Parameter = rep(paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                             nrow(df)),
                             Method = df$method,
                             Index = df$index,
                             L2norm = df$l2norm,
                             Outlier = df$outlier))
    } 
  } else if (par == "gamma1") {
    for (i in gamma1) {
      # simulate data
      data = dataModel(n = n, p = p, d = d, 
                       mu1 = mu1, gamma1 = i, 
                       outlier.index = outlier.index, 
                       sigma0 = sigma0, group = group, 
                       V = V[, 1:d, drop = FALSE], 
                       mu0 = mu0, gamma0 = gamma0, 
                       seed = seed)
      # compute estimated U for various methods
      df = compare.methods(data = data, lambda = lambda, d = d)
      # combine output
      out = rbind(out, 
                  data.frame(Parameter = rep(paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                             nrow(df)),
                             Method = df$method,
                             Index = df$index,
                             L2norm = df$l2norm,
                             Outlier = df$outlier))
    } 
  } else if (par == "n.outlier") {
    for (i in n.outlier) {
      temp = 10+i
      # simulate data
      data = dataModel(n = n, p = p, d = d, 
                       mu1 = mu1, gamma1 = gamma1, 
                       outlier.index = 11:temp, 
                       sigma0 = sigma0, group = group, 
                       V = V[, 1:d, drop = FALSE], 
                       mu0 = mu0, gamma0 = gamma0, 
                       seed = seed)
      # compute estimated U for various methods
      df = compare.methods(data = data, lambda = lambda, d = d)
      # combine output
      out = rbind(out, 
                  data.frame(Parameter = rep(paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                             nrow(df)),
                             Method = df$method,
                             Index = df$index,
                             L2norm = df$l2norm,
                             Outlier = df$outlier))
    } 
  }

  out$Method = as.character(out$Method)
  out$Method[which(out$Method == paste("lambda", round(lambda, 4)))] = "PCATF"
  out$Method = factor(out$Method, 
                      levels = c("Population Mean", "PCA Leverage", "PCATF"))
  
  # plot
  p = ggplot(out, aes(x = Index, y = L2norm, colour = Outlier)) + 
    geom_point() + 
    facet_grid(vars(Method), vars(Parameter)) +
    xlab("Index (Time Point)") +
    ylab("Leverage") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(2,1)]) +
    th
  return(p)
}


plot.multiple = function(n, p, d, mu1, gamma1, 
                         outlier.index = NULL, n.outlier = NULL, 
                         SNR = NULL, sigma0 = NULL,
                         niter, V, group, multiplier, mu0, gamma0,
                         lambda, par = c("SNR", "mu1", "gamma1", "n.outlier")){
  # place holder for the output
  out = data.frame(Parameter = character(),
                   Method = character(),
                   Measure = character(),
                   Value = double())
  if (par == "SNR") {
    for (i in SNR) {
      for (a in 1:niter){
        # simulate data
        data = dataModel(n = n, p = p, d = d, mu1 = mu1, gamma1 = gamma1, 
                         outlier.index = outlier.index, 
                         SNR = i, group = group, 
                         V = V[, 1:d, drop = FALSE], 
                         mu0 = mu0, gamma0 = gamma0)
        df = compare.methods(data = data, lambda = lambda, d = d)
        df$method[which(df$method == paste("lambda", round(lambda, 4)))] = "PCATF"
        df$method = factor(df$method, 
                           levels = c("Population Mean", "PCA Leverage", "PCATF"))
        
        # PCA leverage
        PCAleverage.median = median(df[df$method == "PCA Leverage", "l2norm"])
        PCAleverage.outlier = df[df$method == "PCA Leverage", "l2norm"] >= PCAleverage.median*(1+1.4826*multiplier)
        PCAleverage.precision = ifelse(sum(PCAleverage.outlier) != 0,
                                       sum(PCAleverage.outlier[outlier.index])/sum(PCAleverage.outlier), 0) 
        out = rbind(out, data.frame(Parameter = paste(par, "=", i, sep = ""), 
                                    Method = "PCA Leverage", 
                                    Measure = "Precision", 
                                    Value = PCAleverage.precision))
        PCAleverage.recall = sum(PCAleverage.outlier[outlier.index])/length(outlier.index)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i, sep = ""), 
                                    Method = "PCA Leverage", 
                                    Measure = "Recall", 
                                    Value = PCAleverage.recall))
        
        # our
        our.median = median(df[df$method == "PCATF", "l2norm"])
        our.outlier = df[df$method == "PCATF", "l2norm"] >= our.median*(1+1.4826*multiplier)
        our.precision = ifelse(sum(our.outlier) != 0,
                               sum(our.outlier[outlier.index])/sum(our.outlier), 0)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i, sep = ""), 
                                    Method = "PCATF", 
                                    Measure = "Precision", 
                                    Value = our.precision))
        our.recall = sum(our.outlier[outlier.index])/length(outlier.index)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i, sep = ""), 
                                    Method = "PCATF", 
                                    Measure = "Recall", 
                                    Value = our.recall))
      }
    }
  } else if (par == "mu1") {
    for (i in mu1) {
      for (a in 1:niter){
        # simulate data
        data = dataModel(n = n, p = p, d = d, mu1 = i, gamma1 = gamma1, 
                         outlier.index = outlier.index, 
                         sigma0 = sigma0, group = group, 
                         V = V[, 1:d, drop = FALSE], 
                         mu0 = mu0, gamma0 = gamma0)
        df = compare.methods(data = data, lambda = lambda, d = d)
        df$method[which(df$method == paste("lambda", round(lambda, 4)))] = "PCATF"
        df$method = factor(df$method, 
                           levels = c("Population Mean", "PCA Leverage", "PCATF"))
        
        # PCA leverage
        PCAleverage.median = median(df[df$method == "PCA Leverage", "l2norm"])
        PCAleverage.outlier = df[df$method == "PCA Leverage", "l2norm"] >= PCAleverage.median*(1+1.4826*multiplier)
        PCAleverage.precision = ifelse(sum(PCAleverage.outlier) != 0,
                                       sum(PCAleverage.outlier[outlier.index])/sum(PCAleverage.outlier), 0) 
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCA Leverage", 
                                    Measure = "Precision", 
                                    Value = PCAleverage.precision))
        PCAleverage.recall = sum(PCAleverage.outlier[outlier.index])/length(outlier.index)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCA Leverage", 
                                    Measure = "Recall", 
                                    Value = PCAleverage.recall))
        
        # our
        our.median = median(df[df$method == "PCATF", "l2norm"])
        our.outlier = df[df$method == "PCATF", "l2norm"] >= our.median*(1+1.4826*multiplier)
        our.precision = ifelse(sum(our.outlier) != 0,
                               sum(our.outlier[outlier.index])/sum(our.outlier), 0)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCATF", 
                                    Measure = "Precision", 
                                    Value = our.precision))
        our.recall = sum(our.outlier[outlier.index])/length(outlier.index)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCATF", 
                                    Measure = "Recall", 
                                    Value = our.recall))
      }
    }
  } else if (par == "gamma1") {
    for (i in gamma1) {
      for (a in 1:niter){
        # simulate data
        data = dataModel(n = n, p = p, d = d, mu1 = mu1, gamma1 = i, 
                         outlier.index = outlier.index, 
                         sigma0 = sigma0, group = group, 
                         V = V[, 1:d, drop = FALSE], 
                         mu0 = mu0, gamma0 = gamma0)
        df = compare.methods(data = data, lambda = lambda, d = d)
        df$method[which(df$method == paste("lambda", round(lambda, 4)))] = "PCATF"
        df$method = factor(df$method, 
                           levels = c("Population Mean", "PCA Leverage", "PCATF"))
        
        # PCA leverage
        PCAleverage.median = median(df[df$method == "PCA Leverage", "l2norm"])
        PCAleverage.outlier = df[df$method == "PCA Leverage", "l2norm"] >= PCAleverage.median*(1+1.4826*multiplier)
        PCAleverage.precision = ifelse(sum(PCAleverage.outlier) != 0,
                                       sum(PCAleverage.outlier[outlier.index])/sum(PCAleverage.outlier), 0) 
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCA Leverage", 
                                    Measure = "Precision", 
                                    Value = PCAleverage.precision))
        PCAleverage.recall = sum(PCAleverage.outlier[outlier.index])/length(outlier.index)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCA Leverage", 
                                    Measure = "Recall", 
                                    Value = PCAleverage.recall))
        
        # our
        our.median = median(df[df$method == "PCATF", "l2norm"])
        our.outlier = df[df$method == "PCATF", "l2norm"] >= our.median*(1+1.4826*multiplier)
        our.precision = ifelse(sum(our.outlier) != 0,
                               sum(our.outlier[outlier.index])/sum(our.outlier), 0)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCATF", 
                                    Measure = "Precision", 
                                    Value = our.precision))
        our.recall = sum(our.outlier[outlier.index])/length(outlier.index)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCATF", 
                                    Measure = "Recall", 
                                    Value = our.recall))
      }
    }
  } else if (par == "n.outlier") {
    for (i in n.outlier) {
      temp = 10 + i
      outlier.index = 11:temp
      for (a in 1:niter){
        # simulate data
        data = dataModel(n = n, p = p, d = d, mu1 = mu1, gamma1 = gamma1, 
                         outlier.index = outlier.index, 
                         sigma0 = sigma0, group = group, 
                         V = V[, 1:d, drop = FALSE], 
                         mu0 = mu0, gamma0 = gamma0)
        df = compare.methods(data = data, lambda = lambda, d = d)
        df$method[which(df$method == paste("lambda", round(lambda, 4)))] = "PCATF"
        df$method = factor(df$method, 
                           levels = c("Population Mean", "PCA Leverage", "PCATF"))
        
        # PCA leverage
        PCAleverage.median = median(df[df$method == "PCA Leverage", "l2norm"])
        PCAleverage.outlier = df[df$method == "PCA Leverage", "l2norm"] >= PCAleverage.median*(1+1.4826*multiplier)
        PCAleverage.precision = ifelse(sum(PCAleverage.outlier) != 0,
                                       sum(PCAleverage.outlier[outlier.index])/sum(PCAleverage.outlier), 0) 
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCA Leverage", 
                                    Measure = "Precision", 
                                    Value = PCAleverage.precision))
        PCAleverage.recall = sum(PCAleverage.outlier[outlier.index])/length(outlier.index)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCA Leverage", 
                                    Measure = "Recall", 
                                    Value = PCAleverage.recall))
        
        # our
        our.median = median(df[df$method == "PCATF", "l2norm"])
        our.outlier = df[df$method == "PCATF", "l2norm"] >= our.median*(1+1.4826*multiplier)
        our.precision = ifelse(sum(our.outlier) != 0,
                               sum(our.outlier[outlier.index])/sum(our.outlier), 0)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCATF", 
                                    Measure = "Precision", 
                                    Value = our.precision))
        our.recall = sum(our.outlier[outlier.index])/length(outlier.index)
        out = rbind(out, data.frame(Parameter = paste(par, "=", i,"\n", "SNR=", round(data$SNR, 4), sep = ""), 
                                    Method = "PCATF", 
                                    Measure = "Recall", 
                                    Value = our.recall))
      }
    }
  }

  
  # plot
  p = ggplot(out, aes(x = Method, y = Value, colour = Method, fill = Method)) + 
    geom_boxplot() + 
    facet_grid(vars(Measure), vars(Parameter)) + 
    xlab("") + 
    ylab("") + 
    scale_fill_manual(values = brewer.pal(5, "Set1")[c(3,5)]) +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(3,5)]) +
    th +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank()) +
    stat_summary(geom="crossbar",color="gray70",fatten=1,
                 fun.data = function(x){
                   c(y=median(x), ymin=median(x),ymax=median(x))
                 }) 
  return(p)
}

library(ggrepel)
plot.real <- function(clev, ...){
  choosePCs <- clev$params$choosePCs
  choosePCs_formatted <- switch(choosePCs,
                                mean='Mean',
                                kurtosis='Kurtosis')
  method <- clev$params$method
  method_formatted <- switch(method,
                             leverage='Leverage',
                             robdist='Robust Distance',
                             robdist_subset='Robust Distance Subset')
  measure <- switch(method,
                    leverage=clev$leverage,
                    robdist=clev$robdist,
                    robdist_subset=clev$robdist)
  outliers <- clev$outliers
  cutoffs <- clev$cutoffs
  args <- list(...)
  
  #Log the y-axis if the measurement is robust distance.
  log_measure <- switch(method,
                        leverage=FALSE,
                        robdist=TRUE,
                        robdist_subset=TRUE)
  
  # Identify outliers and their levels of outlyingness.
  index <- 1:length(measure)
  outlier_level_num <- apply(outliers, 1, sum)  # get outlier levels as a single factor
  outlier_level_names <- c('not an outlier', colnames(outliers))
  outlier_level <- factor(outlier_level_names[outlier_level_num + 1], levels=outlier_level_names)
  d <- data.frame(index, measure, outlier_level)
  if(method %in% c('robdist','robdist_subset')){
    d$inMCD <- ifelse(clev$inMCD, 'In MCD', 'Not In MCD')
  }
  
  # The plot will have lines extending downward from outliers
  #  to the x-axis.
  is_outlier <- d$outlier_level != 'not an outlier'
  any_outliers <- any(is_outlier)
  if(any_outliers){
    # Obtain the coordinates of the outliers' lines' vertices.
    drop_line <- d[is_outlier,]
    drop_line$outlier_level <- factor(colnames(outliers)[outlier_level_num[is_outlier]], levels=colnames(outliers)) # remove 'not an outlier' level
    drop_line$xmin <- drop_line$index - .5
    drop_line$xmax <- drop_line$index + .5
    drop_line$ymin <- 0
    drop_line$ymax <- drop_line$measure
    # ggplot will draw rows from top to bottom.
    # Ordering by increasing outlier level ensures lines for the
    #  most outlying observations are drawn last, i.e. on the top layer.
    drop_line <- drop_line[order(drop_line$outlier_level),]
  }
  
  if(log_measure){
    # Add 1 before log transforming to ensure a positive range.
    method_formatted <- paste0('log10(', method_formatted, ' + 1)')
    d$measure <- log(d$measure + 1, base = 10)
    if(any_outliers){
      drop_line$ymax <- log(drop_line$ymax + 1, base = 10)
    }
    cutoffs <- log(cutoffs + 1, base = 10)
  }
  
  # The lowest, middle, and highest outlier levels are colored
  #  yellow, orange, and red, respectively.
  cols <- hsv(h=c(.1,.05,1), s=c(.6,.8,1))
  #if(any_outliers){
  #	cols <- cols[sort(unique(
  #		outlier_level_num[outlier_level_num!=0]))]
  #}
  
  main <- ifelse('main' %in% names(args), args$main,
                 paste0('Outlier Distribution',
                        ifelse(any_outliers, '', ' (None Identified)')))
  sub <- ifelse('sub' %in% names(args), args$sub,
                paste0(choosePCs_formatted,', ',method_formatted))
  xlab <- ifelse('xlab' %in% names(args), args$xlab, 'Index (Time Point)')
  ylab <- ifelse('ylab' %in% names(args), args$ylab, method_formatted)
  legend.position <- ifelse('show.legend' %in% names(args),
                            ifelse(args$show.legend, 'bottom', 'none'),
                            'none')
  if(method=='leverage'){ ylim_max <- 1 }
  else { ylim_max <- max(d$measure) }
  
  plt <- ggplot(d, aes(x=index,y=measure, color=outlier_level))
  if(any_outliers){
    if(method=='leverage'){ nudge_y <- ylim_max * .08 }
    else { nudge_y <- ylim_max * .12 }
    plt <- plt +
      geom_rect(data=drop_line, inherit.aes=FALSE,
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                    fill=outlier_level), alpha=.9) +
      #geom_text(aes(label=ifelse(is_outlier, as.character(index) ,'')), size=4,
      #	nudge_y=nudge_y, check_overlap=TRUE, show.legend=FALSE)
      geom_text_repel(aes(label=ifelse(is_outlier, as.character(index) ,'')), size=4,
                      nudge_y=nudge_y, show.legend=FALSE)
  }
  plt <- plt + geom_point(show.legend=FALSE) +
    scale_color_manual(values=c('grey','black','black','black')) +
    scale_fill_manual(values=cols, labels=colnames(outliers), drop=FALSE) +
    geom_hline(yintercept=cutoffs, linetype='dashed', color='gray') +
    labs(x=xlab, y=ylab, fill='Outlier Level') +
    coord_cartesian(xlim=c(0, floor(max(d$index)*1.02)),
                    ylim=c(0, ylim_max*4)) +
    theme_classic() +
    theme(legend.position=legend.position, panel.spacing.y=unit(1.5, 'lines')) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(main, subtitle=sub) #+ geom_rug(sides='l', col=rgb(.5,0,0,alpha=.2))
  
  if(method %in% c('robdist','robdist_subset')){
    plt <- plt + facet_grid(inMCD~.)
  }
  
  if('type' %in% names(args)){
    if(args$type == 'n'){ return(plt) }
  }
  return(plt)
}
