# ================================
# 1. Load Packages
# ================================
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# ================================
# 2. Load Data
# ================================
data(rotterdam)

# ================================
# 3. Prepare Variables
# ================================
rotterdam$ttr <- rotterdam$dtime
rotterdam$relapse <- rotterdam$death
rotterdam$grp <- as.factor(rotterdam$hormon)

rotterdam$failure <- factor(rotterdam$relapse,
                            levels = c(0,1),
                            labels = c("Censored","Event"))
rotterdam$meno   <- factor(rotterdam$meno,
                           levels=c(0,1),
                           labels=c("Pre","Post"))

rotterdam$chemo  <- factor(rotterdam$chemo,
                           levels=c(0,1),
                           labels=c("No","Yes"))
# ================================
# 4. Check Missing Values
# ================================
colSums(is.na(rotterdam))


rotterdam <- na.omit(rotterdam)

# ================================
# 5. Fix Variable Types
# ================================

# size already categorical ✔
rotterdam$size <- as.factor(rotterdam$size)

# nodes → نحوله category
rotterdam$nodes_cat <- cut(rotterdam$nodes,
                           breaks = c(-1, 3, Inf),
                           labels = c("Low","High"))

# grade → factor
rotterdam$grade <- as.factor(rotterdam$grade)

# ================================
# 6. Event Chart
# ================================
rotterdam$id <- seq_len(nrow(rotterdam))

ggplot(rotterdam, aes(x = ttr, y = id)) +
  geom_segment(aes(x = 0, xend = ttr, y = id, yend = id),
               color = "lightgrey") +
  geom_point(aes(shape = failure), size = 1) +
  labs(title = "Event Chart",
       x = "Time",
       y = "Observation") +
  theme_minimal()

# ================================
# 7. KM Overall
# ================================
km_fit <- survfit(Surv(ttr, relapse) ~ 1, data = rotterdam)

ggsurvplot(km_fit,
           data = rotterdam,
           title = "Overall Survival",
           xlab = "Time",
           surv.median.line = "hv")

# ================================
# 8. KM by Groups (المهم 🔥)
# ================================
group_vars <- c("grp", "size", "grade", "nodes_cat")

for (grp_var in group_vars) {
  
  cat("\n==============================\n")
  cat("Grouping variable:", grp_var, "\n")
  cat("==============================\n")
  
  km_fit_gr <- eval(parse(text = paste0(
    "survfit(Surv(ttr, relapse) ~ ", grp_var,
    ", data = rotterdam)"
  )))
  
  print(summary(km_fit_gr))
  
  cat("\nLog-rank test:\n")
  print(eval(parse(text = paste0(
    "survdiff(Surv(ttr, relapse) ~ ", grp_var,
    ", data = rotterdam)"
  ))))
  
  g <- ggsurvplot(km_fit_gr,
                  data = rotterdam,
                  pval = TRUE,
                  title = paste("Survival by", grp_var),
                  xlab = "Time",
                  surv.median.line = "hv",
                  legend.title = grp_var)
  
  print(g)
}

# ================================
# 9. Cox Model
# ================================

cox_full <- coxph(
  Surv(ttr, relapse) ~
    grp + age + meno + nodes +
    size + grade +
    er + pgr +
    chemo ,
  data = rotterdam
)
summary(cox_full)
cox_model <- coxph(Surv(ttr, relapse) ~  age + nodes + size + grade+pgr,
                   data = rotterdam)
summary(cox_model)
AIC(cox_full, cox_model)
# Tidy output with broom
library(broom)
tidy(cox_model, exponentiate=TRUE, conf.int=TRUE)

# Likelihood Ratio Test statistic
LR_stat <- 2 * (cox_model$loglik[2] - cox_model$loglik[1])

# p-value for LRT
LR_p <- 1 - pchisq(LR_stat, df = length(coef(cox_model)))

# Hazard Ratio table
hr_table <- data.frame(
  Variable = names(coef(cox_model)),
  HR       = exp(coef(cox_model)),
  LowerCI  = summary(cox_model)$conf.int[, "lower .95"],
  UpperCI  = summary(cox_model)$conf.int[, "upper .95"],
  pvalue   = summary(cox_model)$coefficients[, "Pr(>|z|)"]
)

# Final combined output
list(
  Likelihood_Ratio_Stat = LR_stat,
  Likelihood_Ratio_pvalue = LR_p,
  HR_Table = hr_table
)
# ================================
# 10. Hazard Ratios
# ================================
exp(coef(cox_model))   # HR فقط
hr_table <- data.frame(
  Variable = names(coef(cox_model)),
  HR       = exp(coef(cox_model)),
  LowerCI  = summary(cox_model)$conf.int[, "lower .95"],
  UpperCI  = summary(cox_model)$conf.int[, "upper .95"],
  pvalue   = summary(cox_model)$coefficients[, "Pr(>|z|)"]
)

print(hr_table)

# ================================
# 11. Forest Plot
# ================================
ggforest(cox_model,
         data = rotterdam,
         main = "Hazard Ratios")
# cox.zph() tests PH assumption using Schoenfeld residuals
ph_test <-cox.zph(cox_model)
print(ph_test)
# Plot Schoenfeld residuals over time
plot(ph_test)     # Should be flat (no trend) if PH assumption holds











