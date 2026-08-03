rm(list=ls())
library(ggpubr)

sr_delta = seq(0, 40, 5)
mpd_delta = seq(0, 1.0, 0.1)
iri_delta = seq(0, 50, 5)

sr_std_delta = seq(0,25,5)
mpd_std_delta = seq(0, 0.5, 0.05)

est <- function(b1, b2, x, z, cat) {
  cmf = exp(b1*x + b2*z*x)
  df = data.frame(cat=cat, delta=x, cmf=cmf)
}

# SR40 - min
cmf_1 = est(-0.029, 0, sr_delta, 1.3626, "ALL")
cmf_2 = est(-0.057, 0, sr_delta, 1.3626, "TRUCK")
cmf_combine = rbind(cmf_1, cmf_2)

eq_1 = expression(e^{-0.029 %*% Delta*SR40})
eq_2 = expression(e^{-0.057 %*% Delta*SR40})

ggline(cmf_combine, x="delta", y="cmf", group="cat", color="cat", shape="cat",
       numeric.x.axis = TRUE,
       palette=c("#0072B2", "#E69F00", "#CC79A7", "#32CD32")) +
  geom_point(aes(shape = cat, color=cat), size = 4) +
  geom_vline(xintercept = 10, linetype = "dotted", color = "red", linewidth = 0.8) +
  annotate("segment", x=0, xend=10, y=exp(-0.029*10), yend=exp(-0.029*10),
           linetype = "dashed", color = "#0072B2", linewidth = 0.6) +
  annotate("segment", x=0, xend=10, y=exp(-0.057*10), yend=exp(-0.057*10),
           linetype = "dashed", color = "#E69F00", linewidth = 0.6) +
  annotate("text", x=0, y=exp(-0.029*10)+0.02, label="0.748", color="#0072B2", size=5, hjust=0) +
  annotate("text", x=0, y=exp(-0.057*10)+0.02, label="0.566", color="#E69F00", size=5, hjust=0) +
  labs(x = expression(Delta*SR40),
       y = "CMF", 
       title = "CMF for SR40", 
       subtitle = "Pavement Mix Type: Dense") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = c(0.95, 0.5),
        legend.justification = c(1, 0.5),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        plot.caption = element_text(hjust = 0, size=13),
        panel.grid.major.y = element_line(color = "gray",
                                          size = 0.25,
                                          linetype = 2)) +
  guides(color = guide_legend(title = "Crash Type")) + 
  guides(shape = guide_legend(title = "Crash Type")) +
  annotate("text", label=eq_1, x=20, y=0.92, color="#0072B2", size=6) +
  annotate("text", label=eq_2, x=20, y=0.78, color="#E69F00", size=6)

ggsave('C:\\Users\\litabook\\repos\\truck\\pics\\trb_cmf_dense_sr40.png', width = 7, height = 5, device='png', dpi=300)


# MPD
cmf_1 = est(-0.390, 0, mpd_delta, 80, "ALL")
cmf_2 = est(-0.672, 0, mpd_delta, 80, "TRUCK")
cmf_combine = rbind(cmf_1, cmf_2)

eq_1 = expression(e^{-0.390 %*% Delta*MPD})
eq_3 = expression(e^{-0.672 %*% Delta*MPD})


ggline(cmf_combine, x="delta", y="cmf", group="cat", color="cat", shape="cat",
       numeric.x.axis = TRUE,
       palette=c("#0072B2", "#CC79A7", "#32CD32")) +
  geom_point(aes(shape = cat, color=cat), size = 4) +
  geom_vline(xintercept = 0.2, linetype = "dotted", color = "red", linewidth = 0.8) +
  annotate("segment", x=0, xend=0.2, y=exp(-0.390*0.2), yend=exp(-0.390*0.2),
           linetype = "dashed", color = "#0072B2", linewidth = 0.6) +
  annotate("segment", x=0, xend=0.2, y=exp(-0.672*0.2), yend=exp(-0.672*0.2),
           linetype = "dashed", color = "#CC79A7", linewidth = 0.6) +
  annotate("text", x=0, y=exp(-0.390*0.2)+0.02, label="0.925", color="#0072B2", size=5, hjust=0) +
  annotate("text", x=0, y=exp(-0.672*0.2)+0.02, label="0.874", color="#CC79A7", size=5, hjust=0) +
  labs(x = expression(paste(Delta*MPD, " (mm)")),
       y = "CMF", 
       title = "CMF for MPD",
       subtitle = "Pavement Mix Type: Dense") +
  scale_x_continuous(breaks = seq(0, 1.0, 0.2)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = c(0.02, 0.15),
        legend.justification = c(0, 0),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        panel.grid.major.y = element_line(color = "gray",
                                          size = 0.25,
                                          linetype = 2)) +
  guides(color = guide_legend(title = "Crash Type")) + 
  guides(shape = guide_legend(title = "Crash Type")) +
  annotate("text", label=eq_1, x=0.5, y=0.50, color="#0072B2", size=6) +
  annotate("text", label=eq_3, x=0.5, y=0.45, color="#CC79A7", size=6)

ggsave('C:\\Users\\litabook\\repos\\truck\\pics\\trb_cmf_dense_mpd.png', width = 7, height = 5, device='png', dpi=300)


# MPD STD
cmf_1 = est(0.926, 0, mpd_std_delta, 80, "ALL")
cmf_2 = est(1.646, 0, mpd_std_delta, 80, "TRUCK")
cmf_combine = rbind(cmf_1, cmf_2)

eq_1 = expression(e^{0.926 %*% Delta*STD(MPD)})
eq_2 = expression(e^{1.646 %*% Delta*STD(MPD)})

ggline(cmf_combine, x="delta", y="cmf", group="cat", color="cat", shape="cat",
       numeric.x.axis = TRUE,
       palette=c("#0072B2", "#E69F00", "#CC79A7", "#32CD32")) +
  geom_point(aes(shape = cat, color=cat), size = 4) +
  geom_vline(xintercept = 0.1, linetype = "dotted", color = "red", linewidth = 0.8) +
  annotate("segment", x=0, xend=0.1, y=exp(0.926*0.1), yend=exp(0.926*0.1),
           linetype = "dashed", color = "#0072B2", linewidth = 0.6) +
  annotate("segment", x=0, xend=0.1, y=exp(1.646*0.1), yend=exp(1.646*0.1),
           linetype = "dashed", color = "#E69F00", linewidth = 0.6) +
  annotate("text", x=0, y=exp(0.926*0.1)+0.04, label="1.097", color="#0072B2", size=5, hjust=0) +
  annotate("text", x=0, y=exp(1.646*0.1)+0.04, label="1.179", color="#E69F00", size=5, hjust=0) +
  labs(x = expression(paste(Delta*MPD, " (mm)")),
       y = "CMF", 
       title = "CMF for MPD Variation",
       subtitle = "Pavement Mix Type: Dense") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = c(0.05, 0.9),
        legend.justification = c(0, 1),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        panel.grid.major.y = element_line(color = "gray",
                                          size = 0.25,
                                          linetype = 2)) +
  guides(color = guide_legend(title = "Crash Type")) + 
  guides(shape = guide_legend(title = "Crash Type")) +
  annotate("text", label=eq_1, x=0.25, y=2.18, color="#0072B2", size=6) +
  annotate("text", label=eq_2, x=0.25, y=1.98, color="#E69F00", size=6)

ggsave('C:\\Users\\litabook\\repos\\truck\\pics\\trb_cmf_dense_mpd_std.png', width = 7, height = 5, device='png', dpi=300)


# IRI
cmf_1 = est(0.0027, 0, iri_delta, 80, "ALL")
cmf_2 = est(0.0029, 0, iri_delta, 80, "TRUCK")
cmf_combine = rbind(cmf_1, cmf_2)

eq_1 = expression(e^{0.0027 %*% Delta*IRI})
eq_2 = expression(e^{0.0029 %*% Delta*IRI})

ggline(cmf_combine, x="delta", y="cmf", group="cat", color="cat", shape="cat",
       numeric.x.axis = TRUE,
       palette=c("#0072B2", "#E69F00", "#CC79A7", "#32CD32")) +
  geom_point(aes(shape = cat, color=cat), size = 4) +
  geom_vline(xintercept = 10, linetype = "dotted", color = "red", linewidth = 0.8) +
  annotate("segment", x=0, xend=10, y=exp(0.0027*10), yend=exp(0.0027*10),
           linetype = "dashed", color = "#0072B2", linewidth = 0.6) +
  annotate("segment", x=0, xend=10, y=exp(0.0029*10), yend=exp(0.0029*10),
           linetype = "dashed", color = "#E69F00", linewidth = 0.6) +
  annotate("text", x=0, y=exp(0.0027*10)-0.01, label="1.027", color="#0072B2", size=5, hjust=0) +
  annotate("text", x=0, y=exp(0.0029*10)+0.01, label="1.029", color="#E69F00", size=5, hjust=0) +
  labs(x = expression(paste(Delta*IRI, " (in/mi)")),
       y = "CMF", 
       title = "CMF for IRI",
       subtitle = "Pavement Mix Type: Dense") +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  scale_y_continuous(limits = c(1.0, 1.2), breaks = seq(1.0, 1.2, 0.05)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = c(0.05, 0.9),
        legend.justification = c(0, 1),
        legend.title=element_text(size=13),
        legend.text=element_text(size=13),
        panel.grid.major.y = element_line(color = "gray",
                                          size = 0.25,
                                          linetype = 2)) +
  guides(color = guide_legend(title = "Crash Type")) + 
  guides(shape = guide_legend(title = "Crash Type")) +
  annotate("text", label=eq_1, x=25, y=1.17, color="#0072B2", size=6) +
  annotate("text", label=eq_2, x=25, y=1.15, color="#E69F00", size=6)

ggsave('C:\\Users\\litabook\\repos\\truck\\pics\\trb_cmf_dense_iri.png', width = 7, height = 5, device='png', dpi=300)

