library(tourr)
library(tidyverse)
library(GGally)
library(broom)
library(rpart)

n <- 1234
set.seed(404)
d <- tibble(x1 = runif(n, -1, 1),
            x2 = runif(n, -1, 1),
            x3 = runif(n, -1, 1),
            x4 = runif(n, -1, 1))
d <- d %>%
  mutate(y = x1 + x2^2/0.25 - sqrt(abs(x3))/0.15 + rnorm(n, 0, 0.25))
ggscatmat(d)

m1 <- lm(y~., d)
m1_a <- augment(m1)

m2 <- rpart(y~., d)
m2_a <- d %>%
  mutate(.fitted = predict(m2)) %>%
  mutate(.resid = y-.fitted) %>%
  select(y, x1:x4, .fitted, .resid)

m_12 <- bind_cols(m1_a, select(m2_a, .fitted, .resid)) %>%
  rename(m1_f = `.fitted...6`,
         m1_e = `.resid...7`,
         m2_f = `.fitted...12`,
         m2_e = `.resid...13`) %>%
  mutate(de = m1_e - m2_e)


record <- animate(m1_a[,c(2,3,4,5,7)],
                     dependence_tour(c(1,1,1,1,2)),
                     display_xy(axes="bottomleft"),
                     axes="bottomleft",
                     start = matrix(c(1,0,0,0,0,0,0,0,0,1),
                                    ncol=2, byrow=FALSE),
                     rescale=TRUE)

path <- save_history(m1_a[,c(2,3,4,5,7)],
             dependence_tour(c(1,1,1,1,2)),
             start = matrix(c(1,0,0,0,0,0,0,0,0,1),
                            ncol=2, byrow=FALSE),
             max = 50,
             rescale=TRUE)

animate_xy(m1_a[,c(2,3,4,5,7)],
        planned_tour(path),
        axes="bottomleft",
        rescale=TRUE)


# HERE: new script starts
# diagnostic plots
# following calculations are for
# data = (m1_a[,c(2,3,4,5,7)]
library(tidyr)
library(patchwork)

data1 <- as.matrix(m1_a[,c(2,3,4,5,7)])
data2 <- as.matrix(m2_a[,c(2,3,4,5,7)])
data12 <- as.matrix(m_12[,c(2,3,4,5,14)])

# supplementary functions

# get chisq test for independence
# if independent should be uniform on [0,1]x[0,1]
get_p_value_chisq <- function(ndata, k = 5) {
  conf_tab <- table(cut(rank(ndata[,1]), k),
                    cut(rank(ndata[,2]), k))
  chisq.test(conf_tab)$p.value
}

# get pvalues for independence test for the tour
get_pvalues <- function(data, path) {
  chisq_pval <- apply(path, 3, function(x) get_p_value_chisq(data %*% x))
  data.frame(step = 1:dim(path)[3],
             log10_pval = -log10(chisq_pval))
}

# extract projections recorded during the tour
get_projections <- function(path) {
  attr(path, "data") <- NULL
  projections <- t(sapply(1:dim(path)[3],
                          function(x)
                            c(path[-nrow(path),1,x])
  ))

  projections_df <- as.data.frame(projections)
  projections_df$step <- 1:nrow(projections_df)
  pivot_longer(projections_df, -step)
}

projections_long <- get_projections(path)

get_triplots <- function(data, path, datacolor = "orange", modelname = "") {
  chisq_pval_log <- get_pvalues(data, path)
  # where is max dependence
  p_max <- which.max(chisq_pval_log$log10_pval)
  ndata <- data.frame(data %*% matrix(path[,,p_max], ncol = 2))
  # very very ugly solution
  # do not look at it
  ndata1 <- data.frame(data1 %*% matrix(path[,,p_max], ncol = 2))
  ndata2 <- data.frame(data2 %*% matrix(path[,,p_max], ncol = 2))
  ndata12 <- data.frame(data12 %*% matrix(path[,,p_max], ncol = 2))

  # now plots
  pl_coef <- ggplot(projections_long, aes(step, value, color=name)) +
    geom_line() +
    theme_bw() +
    ggtitle("projection weights") +
    geom_vline(xintercept = p_max+0.5, color="red", lty=2) +
    theme(legend.position = "none")

  pl_pval <- ggplot(chisq_pval_log, aes(step, log10_pval)) +
    geom_point(color=datacolor) +
    geom_step(color=datacolor) +
    theme_bw() +
    geom_vline(xintercept = p_max+0.5, color="red", lty=2) +
    ggtitle("log10 p-values for chisq independence test", "large values suggest structure")

  pl_resid <- ggplot(ndata, aes(X1, X2)) +
# very very ugly solution
# do not look at it
    geom_point(data=ndata1, color="black", size=0.1, alpha-0.2) +
    geom_point(data=ndata2, color="orange", size=0.1, alpha-0.2) +
    geom_point(data=ndata12, color="blue", size=0.1, alpha-0.2) +
    geom_point(color="grey", size=0.5, alpha-0.2) +

    geom_smooth(data=ndata1, se = FALSE, color = "black", linewidth=1) +
    geom_smooth(data=ndata2, se = FALSE, color = "orange", linewidth=1) +
    geom_smooth(data=ndata12, se = FALSE, color = "blue", linewidth=1) +

    geom_smooth(se = FALSE, color = datacolor, linewidth=2) +
    theme_bw() +
    xlab("projection") +
    ylab("residuals") +
    ggtitle(paste0(modelname, " vs. projection plot"), paste0("for step ", p_max))

  pl_proj <- ggplot(data = data.frame(var = 1:4, coef =c(path[1:4,1,p_max])),
                    aes(var, coef, fill = factor(var))) +
    geom_col() +
    theme_bw() +
    theme(legend.position = "none") +
    coord_flip()

  list(pl_coef = pl_coef,
        pl_pval = pl_pval,
        pl_resid = pl_resid,
        pl_proj = pl_proj)
}

# calculate actual plots
plots1 <- get_triplots(data1, path, datacolor = "black", modelname = "lm residuals")
plots2 <- get_triplots(data2, path, datacolor = "orange", modelname = "rpart residuals")
plots12 <- get_triplots(data12, path, datacolor = "blue", modelname = "diff residuals")

(plots1$pl_resid + plots2$pl_resid + plots12$pl_resid) /
  (plots1$pl_proj + plots2$pl_proj + plots12$pl_proj) /
(plots1$pl_pval + plots2$pl_pval + plots12$pl_pval) +
  plot_layout(heights = c(3,1,1))



# rest




animate_xy(m1_a[,c(2,3,4,5,7)],
           dependence_tour(c(1,1,1,1,2)),
           axes="bottomleft",
           start = matrix(c(1,0,0,0,0,0,0,0,0,1),
                          ncol=2, byrow=FALSE),
           rescale=TRUE)
animate_xy(m2_a[,c(2,3,4,5,7)],
           dependence_tour(c(1,1,1,1,2)),
           axes="bottomleft",
           start = matrix(c(1,0,0,0,0,0,0,0,0,1),
                          ncol=2, byrow=FALSE),
           rescale=TRUE)

animate(m_12[,c(2,3,4,5,14)],
           dependence_tour(c(1,1,1,1,2)),
           display_xy(axes="bottomleft"),
           start = matrix(c(1,0,0,0,0,0,0,0,0,1),
                          ncol=2, byrow=FALSE),
           rescale=TRUE)


data <- as.matrix(m_12[,c(2,3,4,5,14)])
