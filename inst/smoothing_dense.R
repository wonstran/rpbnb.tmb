rm(list=ls())
require(ggplot2)
library(ggpubr)
library(zoo)

data <- read.csv("C:\\Users\\litabook\\repos\\truck\\data\\export_dense_all.csv")
data$aadt = exp(data$LNAADT)
data$cr = data$C_HV*1000000/(365*data$aadt*0.1)

# Urban Dense Segment
seg_dense = data
sorted = seg_dense[order(seg_dense$SR40_MI3),]
n = round(nrow(sorted)/4-1)
sorted$CMA = rollmean(sorted$cr, k = n, fill=NA, na.pad = TRUE, align="left")
sorted$CMA2 = rollmean(sorted$cr, k = n, fill=NA, na.pad = TRUE, align="center")
sorted$CMA3 = rollmean(sorted$cr, k = n, fill=NA, na.pad = TRUE, align="right")

write.csv(sorted, "C:\\Users\\litabook\\repos\\truck\\results\\smoothing_dense_results.csv", row.names = FALSE)

ylim.prim = c(0, 1)
ylim.sec = c(0, 0.15)
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

ggplot(sorted, aes(x=SR40_MI3)) +
  geom_histogram(aes(x = SR40_MI3, y = (..count..)/sum(..count..)*b+a), bins = 30, fill = "gray", alpha = 0.5) +
  scale_y_continuous(name = "Crash Rates (per 100 million VMT)", 
                     sec.axis = sec_axis(~ (.-a)/b, name = "Friction Proportion"))+
  coord_cartesian(ylim=ylim.prim) +
  geom_line(aes(y = CMA, colour = "CMA")) + 
  geom_line(aes(y = CMA3, colour = "CMA3")) +
  geom_vline(aes(xintercept = 36), color = "red", linetype = "dashed", size=1) +
  geom_vline(aes(xintercept = 56), color = "green4", linetype = "dashed", size=1) +
  geom_label(aes(x=36, y = 0.9, label="Problematic IL: 36"), color="red") +
  geom_label(aes(x=56, y = 0.9, label="Preferred IL: 56"), color="green4") +
  geom_label(aes(x=46, y = 0.15, label="Histogram of Friction"), color="gray28") +
  theme_bw() + 
  theme(panel.border = element_blank(), axis.line = element_line(color = "black")) +
  scale_x_continuous(name="Friction (SR40)", breaks=seq(24, 66, 4), limits = c(22, 68)) +
  ggtitle("Friciton Investigory Levels for Truck Crashes", 
          subtitle = paste("Sample Size =", nrow(sorted), ", Rolling Window =", n)) +
  theme(plot.title = element_text(hjust = 0.5 , margin = margin(0,0,15,0)), legend.position="bottom") +
  ylab("Smoothed Crash Rates (per 100 million VMT)") + 
  scale_color_manual(name="Smooth Alignment:",
                     breaks=c("CMA", "CMA3"),
                     values = c("brown4", "cyan4"),
                     labels = c("Left","Right")) 
  
ggsave('C:\\Users\\litabook\\repos\\truck\\pics\\il_dense_seg_paper.png', width = 7, height = 6, device='png', dpi=300)
