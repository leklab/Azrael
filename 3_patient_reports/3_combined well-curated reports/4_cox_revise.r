#Rscript --vanilla 4_cox_revise.r report_and_SMuRF_high_confidence_dimsum_NEW

library(survival)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1],header=T,sep="\t")
#df <- read.delim('./report_and_SMuRF')

df = df[df$Onset_m != "na" & df$fs1 != "na" & df$fs2 != "na",]

df$log_functional <- log2(2^as.numeric(df$fs1) + 2^as.numeric(df$fs2))
fivenum(df$log_functional)

df$Onset_m = as.numeric(df$Onset_m)
df$event = ifelse(df$Onset_m > 0, 1, 0)

# Fit the Cox proportional hazards model
cox_model = coxph(Surv(Onset_m, event) ~ log_functional, data=df)
summary(cox_model)
