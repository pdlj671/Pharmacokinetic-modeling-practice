# Model and simulation of 
# Minimal Physiologically-Based Pharmacokinetic (mPBPK) Metamodeling of Target Engagement in 
# Skin Informs Anti-IL17A Drug Development in Psoriasis

# Goals:
# 1. Simulate model
# 2. Mimic results from paper
# 3. Sensitivity analysis
# 4. Random forest?

library(mrgsolve)
library(tidyverse)
library(PKPDmisc)
library(PKNCA)
library(knitr)
library(MESS)
library(FME)
# Read model 
# Note that their paper includes random effects
mod <- mread("IL-17A")

# Simulation of various dosing regimens for secukinumab
# Phase II low dose-ranging
# 25 mg every week for four weeks
# 24 hours * 7 = one week
dose <- ev(amt = 25, ii = 4, addl = 5, cmt = "ABS")

# Simulate
out <- mod %>% 
  ev(dose) %>%
  mrgsim(end = 48, delta = 0.1) 

# plot output of interest
plot(out, Cserum + resp~.)


# IV dose with multiple regimens
e1 <- ev(amt = 1000, ii = 28/7, addl = 4, cmt="SERUM") # 4 times a week three times?
e2 <- ev(amt = 1000, cmt ="SERUM")
e3 <- ev(amt = 300, cmt = "SERUM")

# Combine dosing regimens into one event
e <- as_data_set(e1,e2,e3)

# Simulate serum concentration with 3 events

mod %>%
  ev(e) %>%
  mrgsim(end = 392/7, delta = 0.1) %>% 
  plot()


# Local Sensitivity Analysis

# Create function for sensitivity analysis

fun <- function(pars) {
  mod %>%
    param(pars) %>%
    mrgsim_d(e1, delta = 0.0025) %>%
    select(-ID)
}

# Set pars (input) as everything

pars <- as.numeric(param(mod))
pars <- pars[rownames(pars1)]

locSens <- FME::sensFun(
   func = fun,
   parms = pars,
   sensvar = "SERUM",
   tiny= 1e-5
 )

summary(locSens)

plot(locSens)

plot(summary(locSens))

summ <- 
  as_tibble(summary(locSens)) %>%
  mutate(parms = names(pars))


ggplot(data = summ, aes(x=reorder(parms, Mean), y = Mean)) +
  geom_col() +
  labs(x = "Parameter", y = "Coefficient") +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2)


## Calculate AUC and Cmax?

# Set an output for the simulation
# Use carry.out function to obtain concentrations over time
out1 <- mod %>%
  ev(e) %>%
  carry.out("Cserum") %>%
  mrgsim(end = 392/7, delta = 0.1) 


# Set a data frame for only ID, time, and concentration
# There is most likely any easier way to do this

df <- data.frame(out1$ID, out1$time, out1$Cserum)
colnames(df) <- c("ID", "time","Cserum")

# Calculate AUC of each ID
df %>%
  group_by(ID) %>%
  mutate(AUC = auc(time, Cserum, type = "spline"))

# Calculate cmax and tmax

df1<- df %>%
  group_by(ID) %>%
  mutate(Cmax = max(Cserum), Tmax = time[which.max(Cserum)])



