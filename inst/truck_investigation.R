rm(list=ls())
require(ggplot2)
library(ggpubr)
library(scales)
library(zoo)

data <- read.csv("D:\\repos\\SCRIM\\spf_network\\segs\\data\\step_11_mc.csv")
data$CRASH_RATE = data$MC_TOT/(data$AADT*0.1*365*5)*100000000
data$CURVE = ifelse(data$RADIUS_MIN<2000, 1, 0)

K_row<-nrow(data)
data_ord <- data[order(data$SR40_MIN3PTS),]
data_ord$crash_smooth = rollmean(data_ord$CRASH_RATE, k=K_row/5-1, fill=NA, align='left')
ggscatter(data=data_ord,x='SR40_MIN3PTS', y='crash_smooth', size = 1)+
  labs(x = "SR40 Minimum Moving Average",
       y = "Smoothed Crash Rate per 100M VMT", 
       title = "Investigation Level") +
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(panel.grid.major = element_line(linetype="dashed",color="grey", size=0.25))+ 
  scale_x_continuous(breaks = seq(0, 90, 5))
  #geom_vline(xintercept = 50, linetype="solid",color = alpha("red", 0.4), size=1.5)+
  #geom_hline(yintercept = 167, linetype="solid",color = alpha("red", 0.4), size=1.5) 

ggsave('D:\\repos\\SCRIM\\spf_network\\segs\\investigationary level\\mc_sr40.png', width = 7/1.2, height = 5/1.2, device='png', dpi=300)

table(data$SIGNAL)
table(data$NEAR_SIGNAL)
table(data$CURVE)

data_seg = data[(data$SIGNAL==0) & (data$CURVE==0) & (data$NEAR_SIGNAL==0), ]
K_row<-nrow(data_seg)
data_seg_ord <- data_seg[order(data_seg$SR40_MIN3PTS),]
data_seg_ord$crash_smooth = rollmean(data_seg_ord$CRASH_RATE, k=K_row/3-1, fill=NA, align='right')
ggscatter(data=data_seg_ord,x='SR40_MIN3PTS', y='crash_smooth', size = 1)+
  labs(x = "SR40 Minimum Moving Average",
       y = "Smoothed Crash Rate per 100M VMT", 
       title = "Tagenet Segments") +
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(panel.grid.major = element_line(linetype="dashed",color="grey", size=0.25))+ 
  scale_x_continuous(breaks = seq(0, 90, 5))
#geom_vline(xintercept = 50, linetype="solid",color = alpha("red", 0.4), size=1.5)+
#geom_hline(yintercept = 167, linetype="solid",color = alpha("red", 0.4), size=1.5) 

ggsave('D:\\repos\\SCRIM\\spf_network\\segs\\investigationary level\\mc_sr40.png', width = 7/1.2, height = 5/1.2, device='png', dpi=300)


