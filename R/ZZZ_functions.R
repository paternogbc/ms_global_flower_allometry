# Script to make alternative trees using V.PhyloMaker
make_trees <- function(x, scenarios = c("S1", "S2", "S3"), r = 1) {
  cl <- x
  dl <- data.frame(
    species =  cl$tpl_name,
    genus   =  cl$tpl_genus,
    family  =  cl$family
  )
  
  # Make tree
  tree_out <- V.PhyloMaker::phylo.maker(sp.list = dl, scenarios = scenarios, r = r)
  
  return(tree_out)
}

### Function to summarise a modl estimates from a list of RMA models.
tab.rma <- function(mod.list, title = "RMA", model.names){
  ### Arguments check
  ####
  if (class(mod.list) != "list") {
    mod.list <- list(mod.list)
  }
  nm        <- length(mod.list)
  if (length(model.names) != nm) {stop("Check model names")}
  
  tab <- setNames(data.frame(matrix(ncol = 13, nrow = nm)),
                  c("Title", "Subtitle", "n", "R2",
                    "int", "int.se","int.lo", "int.hi",
                    "slo", "slo.se", "slo.lo", "slo.hi",
                    "p_b1"))
  
  ### loop across models
  for (i in 1:nm) {
    
    m <- mod.list[[i]]
    
    title    <- title
    subtitle <- model.names[i]
    n        <- m$n
    r2       <- m$r2[[1]]
    int      <- m$groupsummary$Int
    int.se   <- NA
    int.lo   <- m$groupsummary$Int_lowCI
    int.hi   <- m$groupsummary$Int_highCI
    slo      <- m$groupsummary$Slope
    slo.se   <- NA
    slo.lo   <- m$groupsummary$Slope_lowCI
    slo.hi   <- m$groupsummary$Slope_highCI
    pv.slo   <- m$groupsummary$Slope_test_p
    tab[i,]  <- data.frame(title, subtitle, n, r2, int, int.se, int.lo, int.hi,
                           slo, slo.se, slo.lo, slo.hi, pv.slo, stringsAsFactors = F)
  }
  return(tab)
}

### Summary function for phyrma
summary.phyl.RMA <- function(x) {
  data.frame(
    Slope = x$RMA.beta[2],
    Int = x$RMA.beta[1], 
    model = x$model,
    lambda = x$lambda, 
    r2 = x$test[[1]], 
    Slope_test = x$h0,
    Slope_test_stat = x$test[[2]],
    Slope_test_p = x$test[[4]]
  )
}

### Function to summarise a modl estimates from a list of rma.phy models.
tab.phy.rma <- function(mod.list, title = "RMA.phy", model.names){
  ### Arguments check
  if (class(mod.list) != "list") {
    mod.list <- list(mod.list)
  }
  nm        <- length(mod.list)
  if (length(model.names) != nm) {stop("Check model names")}
  
  tab <- setNames(data.frame(matrix(ncol = 13, nrow = nm)),
                  c("Title", "Subtitle", "n", "R2",
                    "int", "int.se","int.lo", "int.hi",
                    "slo", "slo.se", "slo.lo", "slo.hi",
                    "p_b1"))
  
  ### loop across models
  for (i in 1:nm) {
    m <- mod.list[[i]]
    title    <- title
    subtitle <- model.names[i]
    n        <- length(m$resid)
    r2       <- m$test[[1]]
    int      <- m$RMA.beta[[1]]
    int.se   <- NA
    int.lo   <- NA
    int.hi   <- NA
    slo      <- m$RMA.beta[[2]]
    slo.se   <- NA
    slo.lo   <- NA
    slo.hi   <- NA
    pv.slo   <- m$test[[4]]
    tab[i,]  <- data.frame(title, subtitle, n, r2, int, int.se, int.lo, int.hi,
                           slo, slo.se, slo.lo, slo.hi, pv.slo, stringsAsFactors = F)
  }
  return(tab)
}

### Function to summarise a modl estimates from a list of lm models.
tab.ols <- function(mod.list, title = "OLS", model.names){
  ### Arguments check
  if (class(mod.list) != "list") {
    mod.list <- list(mod.list)
  }
  nm        <- length(mod.list)
  if (length(model.names) != nm) {stop("Check model names")}
  
  tab <- setNames(data.frame(matrix(ncol = 13, nrow = nm)),
                  c("Title", "Subtitle", "n", "R2",
                    "int", "int.se","int.lo", "int.hi",
                    "slo", "slo.se", "slo.lo", "slo.hi",
                    "p_b1"))
  
  ### loop across models
  for (i in 1:nm) {
    m <- mod.list[[i]]
    
    title    <- title
    subtitle <- model.names[i]
    n        <- length(m$residuals)
    r2       <- summary(m)[[8]]
    int      <- m$coefficients[1]
    int.se   <- summary(m)[[4]][1,2]
    int.lo   <- confint(m)[1,1]
    int.hi   <- confint(m)[1,2]
    slo      <- m$coefficients[2]
    slo.se   <- summary(m)[[4]][2,2]
    slo.lo   <- confint(m)[2,1]
    slo.hi   <- confint(m)[2,2]
    pv.slo   <- ttest(x = m, h0 = 1)$p.value
    tab[i,]  <- data.frame(title, subtitle, n, r2, int, int.se, int.lo, int.hi,
                           slo, slo.se, slo.lo, slo.hi, pv.slo, stringsAsFactors = F)
  }
  return(tab)
}

### FUnction to summarise a modl estimates from a list of phylolm models.
tab.phylolm <- function(mod.list, phy, title = "PGLS", model.names){
  ### Arguments check
  tree <- phy
  if (class(tree) != "phylo") {stop("Check tree")}
  if (class(tree) != "phylo") {stop("Check tree")}
  if (class(mod.list) != "list") {
    mod.list <- list(mod.list)
  }
  nm        <- length(mod.list)
  if (length(model.names) != nm) {stop("Check model names")}
  
  tab <- setNames(data.frame(matrix(ncol = 13, nrow = nm)),
                  c("Title", "Subtitle", "n", "R2",
                    "int", "int.se","int.lo", "int.hi",
                    "slo", "slo.se", "slo.lo", "slo.hi",
                    "p_b1"))
  
  ### loop across models
  for (i in 1:nm) {
    m <- mod.list[[i]]
    
    title    <- title
    subtitle <- model.names[i]
    n        <- m$n
    r2       <- R2.pred(mod = m, phy = tree)
    int      <- m$coefficients[1]
    int.se   <- summary(m)[[2]][1,2]
    int.lo   <- m$bootconfint95[,1][1]
    int.hi   <- m$bootconfint95[,1][2]
    slo      <- m$coefficients[2]
    slo.se   <- summary(m)[[2]][2,2]
    slo.lo   <- m$bootconfint95[,2][1]
    slo.hi   <- m$bootconfint95[,2][2]
    pv.slo   <- ttest(x = m, h0 = 1)$p.value
    tab[i,]  <- data.frame(title, subtitle, n, r2, int, int.se, int.lo, int.hi,
                           slo, slo.se, slo.lo, slo.hi, pv.slo, stringsAsFactors = F)
  }
  return(tab)
}

# Function to plot log10 x log10 allometric regressions from model estimates  .
# Intercept is converted to original scale by 10^intercept
plot_phylolm <- function(y,x, data,estimates, subtitle, labx, laby,
                         bsize  = 18,
                         c.ponto = "steelblue",
                         c.line  = "tomato",
                         alp     = .5,
                         s.ponto = 2,
                         s.line  = .8,
                         xa,
                         ya,
                         anota.size = 4.5,
                         axis_size  = 18,
                         intx,
                         inty,
                         limx = NULL,
                         limy = NULL){
  ### Model estimates:
  es     <- estimates
  es     <- filter(es, Subtitle == subtitle)
  if (nrow(es) == 0) {stop("Check subtitle")}
  al     <- (es$int)
  al.lo  <- (es$int.lo)
  al.hi  <- (es$int.hi)
  be     <- es$slo
  be.lo  <- es$slo.lo
  be.hi  <- es$slo.hi
  alr    <- round(10^al, digits = 3)
  ber    <- round(be, digits = 3)
  pv     <- round(es$p_b1, digits = 5) 
  if (pv < 0.00001) {pv <- "0.00001"}
  r2     <- round(es$R2, 2)
  
  ### Plot labels
  lr2 <- paste("~R^2==~", r2, sep = "")
  leq <- paste("~italic(y) ==~", alr, "~italic(x)^", ber)
  #leq <- paste("~beta ==~", deparse(ber), "~','~", "alpha ==~", alr)
  lpv <- paste("P[beta == 1] ==~", deparse(pv), "~ '' ")
  
  ### Plot axis limits
  if (is.null(limx)) {limx <- c(min(data[, x]), max(data[, x]))}
  if (is.null(limy)) {limy <- c(min(data[, y]), max(data[, y]))}
  
  ### Estimate CI intevals
  feq <- function(x, a, b) {a * x^b}
  data$yme <- feq(x = data[[x]], a = 10^al, b = be)
  data$ylo <- feq(x = data[[x]], a = 10^al.lo, b = be.lo)
  data$yhi <- feq(x = data[[x]], a = 10^al.hi, b = be.hi)
  
  g1 <-
    ggplot(data, aes_string(y = y, x = x)) +
    geom_point(size = s.ponto, alpha = alp, fill = c.ponto, shape = 21, color = "black") +
    geom_ribbon(aes(ymin = ylo, ymax = yhi), fill = c.line, alpha = .35) +
    geom_line(aes_string(y = "yme", x = x),
              color = c.line, size = s.line) +
    scale_x_log10(breaks = intx,
                  labels = intx,
                  limits = limx) +
    scale_y_log10(breaks = inty,
                  labels = inty,
                  limits = limy) +
    labs(y = laby, 
         x = labx) +
    annotate("text", x = xa, y = ya, label = leq, parse = T, fontface = 'italic', 
             vjust = 2, hjust = 0, color = gray(.2), size = anota.size) + 
    annotate("text", x = xa, y = ya, label = lr2, parse = T,
             vjust = 4, hjust = 0, size = anota.size,
             color = gray(.2)) +
    annotate("text", x = xa, y = ya, label = lpv, parse = T,
             vjust = 5.5, hjust = -0.03, fontface = 2,
             color = gray(.2), size = anota.size) +
    tema(base_size = bsize) +
    theme(axis.title = element_text(size = axis_size))
  
  return(g1)
}
  
### Function to extract named vectors from data.frame:
extract <- function(x, column){
  setNames(x[,column], nm = rownames(x))
}

### Function to calculate t-test for regression slope with different H0.
ttest <- function(x, h0 = 0){
  if(class(x) == "lm") { 
    beta    = coef(x)[2]
    se.beta = summary(x)[[4]][2,2]
    df      = length(x$residuals) - 2 # df = N - 2 
    tstat <- (beta - h0) / se.beta
    p.val = 2 * pt(abs(tstat), df, lower.tail = FALSE)
    res <- data.frame(t.value = tstat, df = df, p.value = p.val,
                      type = "two-tailed-test")
    
  }
  if(class(x) == "phylolm") { 
    beta    = coef(x)[2]
    se.beta = summary(x)[[2]][2,2]
    df      = length(x$residuals) - 2 # df = N - 2 
    tstat <- (beta - h0) / se.beta
    p.val = 2 * pt(abs(tstat), df, lower.tail = FALSE)
    res <- data.frame(t.value = tstat, df = df, p.value = p.val,
                      type = "two-tailed-test")
  }
  return(res)
}

### Sensitivity analysis for sma models (sampling)
samp_sma <- function(formula, data, breaks = seq(0.05,0.5,0.05), times = 100) {
  d <- data
  n <- nrow(d)
  
  ### LOOP---
  estim <- 
    data.frame(estimate    = numeric(times * length(breaks)),
               slope.test  = numeric(times * length(breaks)),
               intercepts  = numeric(times * length(breaks)),
               n.remove    = numeric(times * length(breaks)),
               n.percent   = numeric(times * length(breaks)),
               repetition  = numeric(times * length(breaks)))
  cc <- 1
  for (i in 1:length(breaks)) {
    for (j in 1:times) {
      nr     <- n - round(n*breaks[i], digits = 0)
      sa.sp  <- sample(x = 1:n, size = nr, replace = F)
      dc     <- d[sa.sp, ]
      m1     <- sma(formula, slope.test = 1, dc)
      slo    <- coef(m1)[[2]]
      int    <- coef(m1)[[1]]
      pslo   <- m1$slopetest[[1]]$p
      
      ### store results
      estim[cc, ] <- c(slo, pslo, int, nr, breaks[i] * 100, j)
      cc <- cc + 1
    }
  }
  estim$n.percent <- as.factor(estim$n.percent)
  return(estim)
}

### Function to predict flower biomass allocation from allometric regressions
pred.flower <- function(estimates, names = c("male", "female", "petals", "sepals"),
                        xint = c(10^seq(-5,0,0.5))) {
  e <- estimates
  feq <- function(x, a, b) {a * x^b}
  ny <- nrow(e)
  nr <- length(xint)
  stopifnot(ny == length(names))
  
  ### Predict biomass between components
  res <- matrix(nrow = nr, ncol = ny)
  
  for(i in 1:ny){
    res[,i] <- feq(x = xint, a = e$int[i],       b =  e$slo[i])
  }
  
  res.raw <- as.data.frame(res)
  colnames(res.raw) <- names
  res.raw$xint <- xint  
  
  res <- as.data.frame(res/xint)
  colnames(res) <- names
  res$xint <- xint
  
  ### Unbiased estimation (from allocation)
  res %>%
    mutate(tot  = male + female + petals + sepals,
           bias = 1 - tot,
           male = male + (male/tot) * bias,
           female = female + (female/tot) * bias,
           petals = petals + (petals/tot) * bias,
           sepals = sepals + (sepals/tot) * bias,
           tot2 = male + female + petals + sepals) -> res.unbiased
  
  ### Unbiased estimation (from raw biomass)
  res.raw %>%
    mutate(tot  = male + female + petals + sepals,
           bias = xint - tot,
           male = male + (male/tot) * bias,
           female = female + (female/tot) * bias,
           petals = petals + (petals/tot) * bias,
           sepals = sepals + (sepals/tot) * bias,
           tot2 = male + female + petals + sepals) -> res.raw.unbiased
  r <- list("biased" = res, "unbiased" = res.unbiased, "raw.unbiased" = res.raw.unbiased)
  return(r)
}

# ggplot theme
tema <- function(base_size = 18){
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(size = 15, hjust = 0),
      plot.subtitle = element_text(size = 15, hjust = 0, face = "italic"),
      #plot.caption = element_text(size = 9, hjust = 1),
      #axis.title = element_text(size = 30),
      #axis.text = element_text(size = 25),
      #axis.ticks = element_line(colour = "red", lineend = "round"),
      #legend.key=element_rect(colour=NA, fill =NA),
      panel.grid  = element_blank(),   
      #panel.border = element_rect(fill = NA, colour = "black", size=1),
      #panel.background = element_rect(fill = "white", colour = "black"), 
      #strip.background = element_rect(fill = "green")
    )
}

# Function to perform permutation test between two allometric slopes (phySMA)
# Assuming log-log regressions between groupvar against xvar
slope_test_phySMA <- function(data, phy, groupvar, xvar, n_sim, method = "BM", 
                              type = c("two-tailed", "one-tailed")){
  # Start
  d0 <- data 
  
  # Full models
  # Original difference in slope
  m01 <- phyl.RMA(x = log10(extract(x = d0, xvar)), 
                  y = log10(extract(x = d0, groupvar[1])), 
                  tree = dp$phy, method = method)
  m02 <- phyl.RMA(x = log10(extract(x = d0, xvar)), 
                  y = log10(extract(x = d0, groupvar[2])), 
                  tree = dp$phy, method = method)
  diff0 <- m01$RMA.beta[2] - m02$RMA.beta[2]
  
  # Loop - null distribution of differences
  diff_null <- numeric(length = n_sim)
  
  for (i in seq_len(n_sim)) {
    
    x <- d0
    swap <- rbernoulli(n = nrow(x))
    
    # get values
    var1 <- x[swap, groupvar[1]]
    var2 <- x[swap, groupvar[2]]
    
    # Swap values
    x[swap, groupvar[1]] <- var2
    x[swap, groupvar[2]] <- var1
    
    m1 <- phyl.RMA(x = log10(extract(x = x, xvar)), 
                   y = log10(extract(x = x, groupvar[1])), 
                   tree = dp$phy, method = method)
    m2 <- phyl.RMA(x = log10(extract(x = x, xvar)), 
                   y = log10(extract(x = x, groupvar[2])), 
                   tree = dp$phy, method = method)
    
    diff <- m1$RMA.beta[2] - m2$RMA.beta[2]
    diff_null[i] <- diff
  }
  if (type == "one-tailed") {
    p.value <- sum((diff_null) >= diff0)/n_sim
    direction = paste(groupvar[1], ">", 
                      groupvar[2])
  }
  if (type == "two-tailed") {
    p.value <- sum(abs(diff_null) >= abs(diff0))/n_sim
    direction = paste(groupvar[1], "> | < ", 
                      groupvar[2])
  }
  
  est <- data.frame(diff = diff0, n_sim, 
                    null_diff = round(mean(diff_null), digits = 4),
                    p.value, type, direction)
  res <- list(slope_test = est,
              null_distribution = data.frame(diff_null))
  return(res)
}

### Sensitivity analysis for SMA/OLS models---------------------------
### Sampling uncertainty
samp_sma <- function(formula, data, breaks = seq(0.05,0.5,0.05), 
                     times = 100, method = "SMA") {
  d <- data
  n <- nrow(d)
  
  ### LOOP---
  estim <- 
    data.frame(estimate      = numeric(times * length(breaks)),
               slope.test  = numeric(times * length(breaks)),
               intercepts  = numeric(times * length(breaks)),
               n.remove    = numeric(times * length(breaks)),
               n.percent   = numeric(times * length(breaks)),
               repetition  = numeric(times * length(breaks)))
  cc <- 1
  for (i in 1:length(breaks)) {
    for (j in 1:times) {
      nr     <- n - round(n*breaks[i], digits = 0)
      sa.sp  <- sample(x = 1:n, size = nr, replace = F)
      dc     <- d[sa.sp, ]
      m1     <- sma(formula, dc, slope.test = 1, method = method)
      slo    <- coef(m1)[[2]]
      int    <- coef(m1)[[1]]
      pslo   <- m1$slopetest[[1]]$p
      ### store results
      estim[cc, ] <- c(slo, pslo, int, nr, breaks[i] * 100, j)
      cc <- cc + 1
    }
  }
  estim$n.percent <- as.factor(estim$n.percent)
  res <- list(estim = estim,
              slope_test = beta_prop_samp(x = estim))
  return(res)
}

### Sensitivity analysis for phySMA models---------------------------
### Sampling uncertainty
samp_physma <- function(y, x, cdata, breaks = seq(0.05,0.5,0.05), 
                     times = 100, method = "BM") {
  d <- cdata
  n <- nrow(d$data)
  
  ### LOOP---
  estim <- 
    data.frame(estimate      = numeric(times * length(breaks)),
               slope.test  = numeric(times * length(breaks)),
               intercepts  = numeric(times * length(breaks)),
               n.remove    = numeric(times * length(breaks)),
               n.percent   = numeric(times * length(breaks)),
               repetition  = numeric(times * length(breaks)))
  cc <- 1
  for (i in 1:length(breaks)) {
    for (j in 1:times) {
      nr     <- n - round(n*breaks[i], digits = 0)
      sa.sp  <- sample(x = 1:n, size = nr, replace = F)
      dr     <- d$data[sa.sp, ]
      dc     <- match_dataphy(get(y) ~ get(x), data = dr, phy = d$phy, verbose = F)
      m1     <- phyl.RMA(x = log10(extract(x = dc$data, x)),
                         y = log10(extract(x = dc$data, y)), 
                         tree = dc$phy, method = method)
      slo    <- m1$RMA.beta[[2]]
      int    <- m1$RMA.beta[[1]]
      pslo   <- m1$test[[4]]
      ### store results
      estim[cc, ] <- c(slo, pslo, int, nr, breaks[i] * 100, j)
      cc <- cc + 1
    }
  }
  estim$n.percent <- as.factor(estim$n.percent)
  res <- list(estim = estim,
              slope_test = beta_prop_samp(x = estim))
  return(res)
}

### Sensitivity analysis for PGLS models---------------------------
### Sampling uncertainty
samp_pgls <- function(formula, cdata, breaks = seq(0.05,0.5,0.05), 
                        times = 100, method = "BM") {
  d <- cdata
  n <- nrow(d$data)
  
  ### LOOP---
  estim <- 
    data.frame(estimate      = numeric(times * length(breaks)),
               slope.test  = numeric(times * length(breaks)),
               intercepts  = numeric(times * length(breaks)),
               n.remove    = numeric(times * length(breaks)),
               n.percent   = numeric(times * length(breaks)),
               repetition  = numeric(times * length(breaks)))
  cc <- 1
  for (i in 1:length(breaks)) {
    for (j in 1:times) {
      nr     <- n - round(n*breaks[i], digits = 0)
      sa.sp  <- sample(x = 1:n, size = nr, replace = F)
      dr     <- d$data[sa.sp, ]
      dc     <- match_dataphy(formula, data = dr, phy = d$phy, verbose = F)
      m1     <- phylolm(formula, data = dc$data, phy = dc$phy,
                        model = method)
      slo    <- coef(m1)[[2]]
      int    <- coef(m1)[[1]]
      pslo   <- ttest(x = m1, h0 = 1)$p.value
      ### store results
      estim[cc, ] <- c(slo, pslo, int, nr, breaks[i] * 100, j)
      cc <- cc + 1
    }
  }
  estim$n.percent <- as.factor(estim$n.percent)
  res <- list(estim = estim,
              slope_test = beta_prop_samp(x = estim))
  return(res)
}

# Proportion of P.values == 1 for samp functions
beta_prop_samp <- function(x) {
    x %>%
    group_by(n.percent) %>%
    summarise(n = n(),
              proportion_b1 = 1 - sum(slope.test < 0.05)/n) %>%
    as.data.frame() %>% 
    ungroup()
}
  
# Plot sampling uncertainty-----------------------------------------
plot_sensi_samp <- function(x, color, fill, st){
  ggplot(x$estim, aes(y = estimate, x = factor(n.percent))) +
    geom_hline(yintercept = es$slo[1], size = 1, color = color) +
    geom_jitter(width = 0.05, alpha = .10, show.legend = F, color = color) +
    geom_boxplot(width = 0.25, alpha = .55, outlier.shape = NA,
                 show.legend = F, fill = fill) +
    geom_hline(yintercept = 1, lty = "dashed") +
    scale_y_continuous(breaks = bre, limits = lims) +
    theme_light(base_size = 18) +
    labs(x = "% of species removed", y = "Estimated slope",
         subtitle = paste(st)) +
    tema()
} 

# Function for phylogenetic uncertainty (phySMA)-----------
tree_physma <- function(y, x, data, trees, method = "BM"){
  estim <- data.frame(
    estimate = numeric(length(trees)),
    slope.test  = numeric(length(trees))
    )
  
  for (i in seq_len(length(trees))) {
    m1 <- phyl.RMA(x = log10(extract(x = data, x)), 
                   y = log10(extract(x = data, y)), 
                   tree = trees[[i]], method = method)
    estim[i, ] <- c(coef(m1)[[2]], m1$test[[4]])
  }
  return(estim)
}
