## Interactions using Cox PH Models


## Load packages
library(survival)
library(car)
library(multcomp)


####################################
##  2 BINARY VARIABLES
####################################

# Key variables in model
DATA$VAR1 = ifelse(DATA$substrate_cleavage_class == "SC > 25 %", 1, 0)

DATA$VAR2 = ifelse(DATA$TNM_GRP==1, 1, 0)

# Set-up model
surv.func = Surv(TIME_OS, EVENT_OS) ~ VAR1*VAR2

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
               c(exp(confint(con.2)$confint[1:3]), summary(con.2)$test$pvalues[1])) )
info = data.frame(info)

colnames(info) = c(paste(c("HR", "l", "u", "p"), "VAR1.0", sep="."), 
                   paste(c("HR", "l", "u", "p"), "VAR1.1", sep="."))
rownames(info) = c("VAR2.0", "VAR2.1")


# HRs, 95% CI, p-values: Compare VAR1 within VAR2 groups 
mat.mod.1 = matrix(c(1, 0, 0), 1) # VAR1 1 vs. 0 among VAR2 0
mat.mod.2 = matrix(c(1, 0, 1), 1) # VAR1 1 vs. 0 among VAR2 1

con.1 = glht(surv.mod, linfct = mat.mod.1)
con.2 = glht(surv.mod, linfct = mat.mod.2)

info = cbind(rep("ref", 2),
             rbind(c(exp(confint(con.1)$confint[1:3]), summary(con.1)$test$pvalues[1]), 
                   c(exp(confint(con.2)$confint[1:3]), summary(con.2)$test$pvalues[1])) )
info = data.frame(info)

colnames(info) = c("VAR1.0", paste(c("HR", "l", "u", "p"), "VAR1.1", sep="."))
rownames(info) = c("VAR2.0", "VAR2.1")


####################################
##  BINARY VARIABLE & CATEGORICAL VARIABLE
####################################

# Key variables in model
DATA$VAR1 = ifelse(DATA$substrate_cleavage_class == "SC > 25 %", 1, 0)

DATA$VAR2 = factor(DATA$TNM)

# Set-up model
surv.func = Surv(TIME_OS, EVENT_OS) ~ VAR1*VAR2

surv.mod = coxph(surv.func, data=DATA)

# Interaction p-value
type3.test = Anova(surv.mod, type="III")

int.p = type3.test$"Pr(>Chisq)"[3] 

# HRs, 95% CI, p-values: Compare VAR2 within VAR1 groups 
mat.mod.1 = matrix(c(0, 1, 0, 0, 0), 1) # VAR2 2 vs. 1 among VAR1 0
mat.mod.2 = matrix(c(0, 0, 1, 0, 0), 1) # VAR2 3 vs. 1 among VAR1 0
mat.mod.3 = matrix(c(0, 1, 0, 1, 0), 1) # VAR2 2 vs. 1 among VAR1 1
mat.mod.4 = matrix(c(0, 0, 1, 0, 1), 1) # VAR2 3 vs. 1 among VAR1 1

con.1 = glht(surv.mod, linfct = mat.mod.1)
con.2 = glht(surv.mod, linfct = mat.mod.2)
con.3 = glht(surv.mod, linfct = mat.mod.3)
con.4 = glht(surv.mod, linfct = mat.mod.4)

info = rbind(rep(c("ref", rep("-",3)), 2),
             rbind(c(c(exp(confint(con.1)$confint[1:3]), summary(con.1)$test$pvalues[1]), 
                     c(exp(confint(con.3)$confint[1:3]), summary(con.3)$test$pvalues[1])),
                   c(c(exp(confint(con.2)$confint[1:3]), summary(con.2)$test$pvalues[1]), 
                     c(exp(confint(con.4)$confint[1:3]), summary(con.4)$test$pvalues[1])) ) )
info = data.frame(info)

colnames(info) = c(paste(c("HR", "l", "u", "p"), "VAR1.0", sep="."), 
                   paste(c("HR", "l", "u", "p"), "VAR1.1", sep="."))
rownames(info) = c("VAR2.1", "VAR2.2", "VAR2.3")


# HRs, 95% CI, p-values: Compare VAR1 within VAR2 groups  
mat.mod.1 = matrix(c(1, 0, 0, 0, 0), 1) # VAR1 1 vs. 0 among VAR2 1
mat.mod.2 = matrix(c(1, 0, 0, 1, 0), 1) # VAR1 1 vs. 0 among VAR2 2
mat.mod.3 = matrix(c(1, 0, 0, 0, 1), 1) # VAR1 1 vs. 0 among VAR2 3

con.1 = glht(surv.mod, linfct = mat.mod.1)
con.2 = glht(surv.mod, linfct = mat.mod.2)
con.3 = glht(surv.mod, linfct = mat.mod.3)

info = cbind(rep("ref", 3),
             rbind(c(exp(confint(con.1)$confint[1:3]), summary(con.1)$test$pvalues[1]), 
                   c(exp(confint(con.2)$confint[1:3]), summary(con.2)$test$pvalues[1]),
                   c(exp(confint(con.3)$confint[1:3]), summary(con.3)$test$pvalues[1])) )
info = data.frame(info)

colnames(info) = c("VAR1.0", paste(c("HR", "l", "u", "p"), "VAR1.1", sep="."))
rownames(info) = c("VAR2.1", "VAR2.2", "VAR2.3")





