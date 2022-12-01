# Import libraries
library(data.table)
library(ggplot2)
library(patchwork)
library(ggsci)

paths <- Sys.glob("~/Stanford/CARE 2022/results/*/adjusted-mortality-rates.rds")
dt.list <- lapply(paths, function(x){
  return(readRDS(x))
})
dt <- rbindlist(dt.list)

dt[, SEX := factor(SEX, levels = c("M", "F"))]
dt[CAUSE == "Flu_pneumonia", CAUSE := "Flu/Pneumonia"]
cause.order <- c("Cancer", "Heart", "Cerebrovascular", "Accidents", "Flu/Pneumonia", "Suicide",
                 "Diabetes", "Respiratory", "Alzheimer", "Kidney", "Liver")
dt[, CAUSE := factor(CAUSE, levels = cause.order)]

ggplot(dt, aes(YEAR, M)) +
  geom_line(aes(color = SEX)) +
  geom_ribbon(aes(ymin = CL, ymax = CU, fill = SEX), alpha = 0.4) +
  geom_point(aes(y = AMR, color = SEX), size = 0.7, alpha = 0.4) +
  facet_wrap(~CAUSE, scales = "free_y") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,NA)) +
  scale_color_lancet() +
  scale_fill_lancet() +
  labs(x = "Year", y = "Adjusted Mortality Rates [per 100000]",
       color = "Sex", fill = "Sex") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.position = c(0.85,0.15),
    axis.title = element_text(size = 9)
  )
