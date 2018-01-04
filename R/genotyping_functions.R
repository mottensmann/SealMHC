#' genotype using an modified doc approach (Lighten et al. 2014)
#' 
#' @description 
#' (I) Calculate cumulative sequencing depth
#' (II) Compute rate of change (ROC, first derivative)
#' (III) Compute degree of change (DOC) to get inflection points
#' 
#' @details 
#' Mathematically, inflections points are identical for 1 or 2 alleles, therefore
#' the two are discriminated based on the gain in sequencing depth from 1 to 2 
#' (i.e. stay with 1 allele, if step to 2 is below 'gain')
#' 
#' @param x vector of reads per variant
#' @param n numbers of most abundant variants to consider. 
#' @param names names of variants
#' @param gain relative increase in sequencing depth standardised over total depth
#' @param doc_min minimal required degree of change giving the significane level of inflection points
#' @param depth minimum relative sequencing depth of alleles. Value between 0 and 1

get_genotypes <- function(x = NULL, n = 8, names = NULL, plot = F, gain = 0.05, doc_min = 40, depth_min = 0.7) {
  # 1. Calculate degree of changes per for n variants
  # 1.1 Cumulative sequence depth
  
  # ensure that n <= nmax
  if (n > length(x)) n <- length(x)
  cum_dep <- cumsum(sort(x, decreasing = T)[1:n])
  
  # 1.2. Rate of change i.e. first derivative
  roc <- diff(cum_dep)
  
  # 1.3 Degree of change i.e. second derivative
  doc <- rep(0, length(roc) - 1)
  for (i in 1:length(doc)) doc[i] <-  roc[i]/roc[i + 1]
  
  # 1.4 Standardise doc
  doc <- doc/sum(doc)*100
  
  # 2. Get inflection point of curve
  infl <- which(doc == max(doc))[1] + 1
  
  # 3. Distinguish between 1 and 2 alleles
  # Assume one allele, when raise in cumsum is less than 5 %
  
  if (infl == 2) {
    if ((cum_dep[2]/max(cum_dep)) - (cum_dep[1]/max(cum_dep)) < gain)  {
      infl <- 1 
    }
  }
  loop <- "Start"
  if (infl > 2) {
    while (loop != "Stop") {
      if ((cum_dep[infl]/max(cum_dep)) - (cum_dep[infl - 1]/max(cum_dep)) < gain) {
        infl <- infl - 1 
      } else {
        loop <- "Stop"
      }
      if (infl <= 1) {
        loop <- "Stop"
      }
    }
  }
  # 3.1 Set quality score to high initially
  quality_score <- "High"
  
  # data frame for plotting
  df2 <- data.frame(x = 0:length(cum_dep),
                    y = c(0, cum_dep/max(cum_dep)*100),
                    quality = "High")
  
  # Linear models
  lm_fit <- with(df2, lm(y~x))
  slope <- lm_fit$coefficients[[2]]
  intercept <- lm_fit$coefficients[[1]]
  
  # correlation between depth and order of variants
  cor_val <- with(df2, cor(x,y))
  if (cor_val > 0.90) {
    df2$quality <- "Low" 
    quality_score <- "Low"
  }
  
  # Check Depth of variants
  if (cum_dep[[infl]]/max(cum_dep) < depth_min) {
    df2$quality <- "Low"
    quality_score <- "Low"
  }
  
  # Require doc to be significant i.e. > 0.4
  if (max(doc) < doc_min & infl %in% 2:3) {
    df2$quality <- "Low"
    quality_score <- "Low"
  }
  
  
  # If amplicon quality is low, then:
  # a) Cumsum is straight line i.e. r > 0.85
  # b) Gain to n-1 is small AND
  # c) increase in depth marginal i.e < 10% of n-1
  
  # 4. Score alleles
  alleles <- names(x)[order(x, decreasing = T)][1:infl]
  df <- 
    data.frame(alleles = paste0(alleles, collapse = ";"),
               total_depth = max(cum_dep),
               allele_depth = cum_dep[[infl]]/max(cum_dep),
               doc = ifelse(infl == 1, NA, max(doc)),
               n_alleles = infl,
               quality = quality_score)
  
  if (plot == TRUE) {
    
    
    plot <- ggplot(df2, aes(x = x,y = y, shape = quality)) +
      geom_point() +
      geom_line() +
      theme_classic() +
      xlab("Number of variants") +
      ylab("Cumulative sequencing depth [%]") +
      scale_y_continuous(expand = c(0,0),
                         breaks = seq(0,100,10)) +
      scale_x_continuous(expand = c(0,0),
                         breaks = 0:16) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14)) +
      geom_vline(xintercept = infl,
                 linetype = "dashed", 
                 color = "darkorange",
                 size = 1.2) +
      geom_abline(slope = slope,
                  intercept = intercept,
                  linetype = "dotted",
                  color = "blue") +
      annotate("text", x = n - 1, y = 20,
               label = paste0("R: ", round(cor_val, 2)),
               color = "blue")
    
    return(plot)
  } else {
    return(list(df = df,
                alleles = alleles,
                coord = 
                  data.frame(x = 1:n,
                             y = cum_dep/max(cum_dep)*100,
                             group = infl,
                             quality = df$quality)))
  }
  
  
}

