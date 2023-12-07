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

# !!
# I had to modify the animate() function
# to keep record of consecutive frames
# added line 80: record[[i]] <- list(i = i, step = step)

# HERE: new script starts
# diagnostic plots
# following calculations are for
# data = (m1_a[,c(2,3,4,5,7)]
library(tidyr)
library(patchwork)

data <- as.matrix(m1_a[,c(2,3,4,5,7)])

# supplementary functions

# get chisq test for independence
# if independent should be uniform on [0,1]x[0,1]
get_p_value_chisq <- function(ndata, k = 5) {
  conf_tab <- table(cut(rank(ndata[,1]), k),
                    cut(rank(ndata[,2]), k))
  chisq.test(conf_tab)$p.value
}

# get pvalues for independence test for the tour
get_pvalues <- function(data, record) {
  chisq_pval <- sapply(record, function(x) get_p_value_chisq(data %*% x$step$proj))
  data.frame(step = seq_along(record),
             log10_pval = -log10(chisq_pval))
}

# extract projections recorded during the tour
get_projections <- function(record) {
  projections <- t(sapply(record,
                          function(x) {
                            proj <- x$step$proj
                            proj[-nrow(proj), 1] # first column without last row
                          }))

  projections_df <- as.data.frame(projections)
  projections_df$step <- 1:nrow(projections_df)
  pivot_longer(projections_df, -step)
}

# calculate projections and pvalues
projections_long <- get_projections(record)
chisq_pval_log <- get_pvalues(data, record)
# where is max dependence
p_max <- which.max(chisq_pval_log$log10_pval)
ndata <- data.frame(data %*% record[[p_max]]$step$proj)

# now plots
pl_coef <- ggplot(projections_long, aes(step, value, color=name)) +
  geom_line() +
  theme_bw() +
  ggtitle("projection weights") +
  geom_vline(xintercept = p_max, color="red", lty=2) +
  theme(legend.position = "none")


pl_pval <- ggplot(chisq_pval_log, aes(step, log10_pval)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_vline(xintercept = p_max, color="red", lty=2) +
  ggtitle("log10 p-values for chisq independence test", "large values suggest structure")

pl_resid <- ggplot(ndata, aes(X1, X2)) +
  geom_point() +
  geom_smooth(se = FALSE, color = "orange", linewidth=2) +
  theme_bw() +
  xlab("projection") +
  ylab("residuals") +
  ggtitle("residual vs. projection plot", paste0("for step ", p_max))

pl_proj <- ggplot(data = data.frame(var = 1:4, coef = record[[p_max]]$step$proj[1:4,1]),
                    aes(var, coef, fill = factor(var))) +
                   geom_col() +
  theme_bw() +
  theme(legend.position = "none") +
  coord_flip()


(pl_pval + pl_resid) /
  (pl_coef + pl_proj) +
  plot_layout(heights = c(3,1))



# rest




record <- animate_xy(m1_a[,c(2,3,4,5,7)],
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
