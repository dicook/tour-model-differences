library(tourr)
library(tidyverse)
library(GGally)
library(broom)
library(rpart)

n <- 123
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

m_12 <- bind_cols(m1_a, select(m2_a, .fitted, .resid)) %>%
  rename(m1_f = `.fitted...6`,
         m1_e = `.resid...7`,
         m2_f = `.fitted...12`,
         m2_e = `.resid...13`) %>%
  mutate(de = m1_e - m2_e)

animate_xy(m_12[,c(2,3,4,5,14)],
           dependence_tour(c(1,1,1,1,2)),
           axes="bottomleft",
           start = matrix(c(1,0,0,0,0,0,0,0,0,1),
                          ncol=2, byrow=FALSE),
           rescale=TRUE)

