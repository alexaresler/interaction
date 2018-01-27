####################################
## Interactions using Cox PH Models
####################################


## Load packages
library(survival)
library(survsim)
library(car)
library(multcomp)


####################################
##  SIMULATE DATA
####################################
N = 1000
DATA = simple.surv.sim(n=N, foltime = 10, 
                       dist.ev="weibull", anc.ev=0.8, beta0.ev=2, 
                       dist.cens="weibull", anc.cens=1.2, beta0.cens=4)
DATA$TIME = DATA$stop
DATA$EVENT = DATA$status

DATA$VAR1 = ifelse(runif(N)<.5, 1, 0)
DATA$VAR2 = ifelse(runif(N)<.5, 1, 0)
DATA$VAR3 = factor(sample(1:3, N, replace=TRUE, prob=c(0.4, 0.3, 0.3)))


####################################
##  TWO BINARY VARIABLES
####################################

# Set-up model
surv.func = Surv(TIME, EVENT) ~ VAR1*VAR2

surv.mod = coxph(surv.func, data=DATA)

# Interaction p-value
type3.test = Anova(surv.mod, type="III")

int.p = type3.test$"Pr(>Chisq)"[3] 

# HRs, 95% CI, p-values: Compare VAR2 within VAR1 groups 
mat.mod.1 = matrix(c(0, 1, 0), 1) # VAR2 1 vs. 0 among VAR1 0
mat.mod.2 = matrix(c(0, 1, 1), 1) # VAR2 1 vs. 0 among VAR1 1

con.1 = glht(surv.mod, linfct = mat.mod.1)
con.2 = glht(surv.mod, linfct = mat.mod.2)

info = rbind(rep(c("ref", rep("-",3)), 2),
             c(c(exp(confint(con.1)$confint[1:3]), summary(con.1)$test$pvalues[1]), 
               c(exp(confint(con.2)$confint[1:3]), summary(con.2)$test$pvalues[1])))
info = data.frame(info)

colnames(info) = c(paste(c("HR", "l", "u", "p"), "VAR1.0", sep="."), 
                   paste(c("HR", "l", "u", "p"), "VAR1.1", sep="."))
rownames(info) = c("VAR2.0", "VAR2.1")

info

# HRs, 95% CI, p-values: Compare VAR1 within VAR2 groups 
mat.mod.1 = matrix(c(1, 0, 0), 1) # VAR1 1 vs. 0 among VAR2 0
mat.mod.2 = matrix(c(1, 0, 1), 1) # VAR1 1 vs. 0 among VAR2 1

con.1 = glht(surv.mod, linfct = mat.mod.1)
con.2 = glht(surv.mod, linfct = mat.mod.2)

info = cbind(rep("ref", 2),
             rbind(c(exp(confint(con.1)$confint[1:3]), summary(con.1)$test$pvalues[1]), 
                   c(exp(confint(con.2)$confint[1:3]), summary(con.2)$test$pvalues[1])))
info = data.frame(info)

colnames(info) = c("VAR1.0", paste(c("HR", "l", "u", "p"), "VAR1.1", sep="."))
rownames(info) = c("VAR2.0", "VAR2.1")

info


####################################
##  BINARY VARIABLE & CATEGORICAL VARIABLE
####################################

# Set-up model
surv.func = Surv(TIME, EVENT) ~ VAR1*VAR3

surv.mod = coxph(surv.func, data=DATA)

# Interaction p-value
type3.test = Anova(surv.mod, type="III")

int.p = type3.test$"Pr(>Chisq)"[3] 

# HRs, 95% CI, p-values: Compare VAR3 within VAR1 groups 
mat.mod.1 = matrix(c(0, 1, 0, 0, 0), 1) # VAR3 2 vs. 1 among VAR1 0
mat.mod.2 = matrix(c(0, 0, 1, 0, 0), 1) # VAR3 3 vs. 1 among VAR1 0
mat.mod.3 = matrix(c(0, 1, 0, 1, 0), 1) # VAR3 2 vs. 1 among VAR1 1
mat.mod.4 = matrix(c(0, 0, 1, 0, 1), 1) # VAR3 3 vs. 1 among VAR1 1

con.1 = glht(surv.mod, linfct = mat.mod.1)
con.2 = glht(surv.mod, linfct = mat.mod.2)
con.3 = glht(surv.mod, linfct = mat.mod.3)
con.4 = glht(surv.mod, linfct = mat.mod.4)

info = rbind(rep(c("ref", rep("-",3)), 2),
             rbind(c(c(exp(confint(con.1)$confint[1:3]), summary(con.1)$test$pvalues[1]), 
                     c(exp(confint(con.3)$confint[1:3]), summary(con.3)$test$pvalues[1])),
                   c(c(exp(confint(con.2)$confint[1:3]), summary(con.2)$test$pvalues[1]), 
                     c(exp(confint(con.4)$confint[1:3]), summary(con.4)$test$pvalues[1]))))
info = data.frame(info)

colnames(info) = c(paste(c("HR", "l", "u", "p"), "VAR1.0", sep="."), 
                   paste(c("HR", "l", "u", "p"), "VAR1.1", sep="."))
rownames(info) = c("VAR3.1", "VAR3.2", "VAR3.3")

info

# HRs, 95% CI, p-values: Compare VAR1 within VAR3 groups  
mat.mod.1 = matrix(c(1, 0, 0, 0, 0), 1) # VAR1 1 vs. 0 among VAR3 1
mat.mod.2 = matrix(c(1, 0, 0, 1, 0), 1) # VAR1 1 vs. 0 among VAR3 2
mat.mod.3 = matrix(c(1, 0, 0, 0, 1), 1) # VAR1 1 vs. 0 among VAR3 3

con.1 = glht(surv.mod, linfct = mat.mod.1)
con.2 = glht(surv.mod, linfct = mat.mod.2)
con.3 = glht(surv.mod, linfct = mat.mod.3)

info = cbind(rep("ref", 3),
             rbind(c(exp(confint(con.1)$confint[1:3]), summary(con.1)$test$pvalues[1]), 
                   c(exp(confint(con.2)$confint[1:3]), summary(con.2)$test$pvalues[1]),
                   c(exp(confint(con.3)$confint[1:3]), summary(con.3)$test$pvalues[1])))
info = data.frame(info)

colnames(info) = c("VAR1.0", paste(c("HR", "l", "u", "p"), "VAR1.1", sep="."))
rownames(info) = c("VAR3.1", "VAR3.2", "VAR3.3")

info
