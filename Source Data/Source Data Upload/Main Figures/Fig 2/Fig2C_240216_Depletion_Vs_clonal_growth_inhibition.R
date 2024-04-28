library('here')
data <- read.csv(here("/Users/zisan/Files/CRISPRi-Seq/Rhamnose IC50_clonal/Rhamnose IC Values/IC50 Values_depletion.csv"))
head(data)
library('ggplot2')
#install.packages('ggMarginal')
library("ggpubr")
library(ggExtra)
library(cowplot)
library(dplyr)

# Extract the lowest value for each sgRNA
lowest<-# Extract the lowest value for each sgRNA in terms of OD
  lowest_values <- data %>%
  group_by(sgRNA) %>%
  summarize(lowest_OD = min(OD, na.rm = TRUE))

# Merge with the original data to find corresponding Depletion value
merged_data <- merge(lowest_values, data, by = "sgRNA")

# Select only relevant columns
final_data <- merged_data[, c("sgRNA", "lowest_OD", "Depletion")]
head(final_data)

# Save the lowest values as a new CSV file
write.csv(final_data, "/Users/zisan/Files/CRISPRi-Seq/Rhamnose IC50_clonal/Rhamnose IC Values/lowest_OD_and_Depletion_values.csv", row.names = FALSE)

data2 <- read.csv(here("/Users/zisan/Files/CRISPRi-Seq/Rhamnose IC50_clonal/Rhamnose IC Values/lowest_OD_and_Depletion_values.csv"))


a<-ggscatter(data2, x = "Depletion", y = "lowest_OD2",
          color = "sgRNA", shape = 19, alpha=0.7,size = 3, # Points color, shape and size
          #conf.int = TRUE, 
          add.params = list(color = "darkblue",
                            fill = "lightgray"),
          #cor.coef = TRUE, cor.method = "pearson",
          xlab = "Depletion (Pool)", ylab = "x.expression")+
  geom_smooth(method = "lm", color="blue4") +
  theme(legend.position = "none")
print(a)
b<-a+scale_color_manual(values=c("1828"="grey",
                                 "1831"="grey",
                                 "1832"="grey",
                                 "1834"="grey",
                                 "1835"="grey",
                                 "1836"="grey",
                                 "1837"="grey",
                                 "1838"="grey",
                                 "1404"="grey",
                                 "1405"="grey",
                                 "1406"="grey",
                                 "1407"="grey",
                                 "1408"="grey",
                                 "1412"="grey",
                                 "1413"="grey",
                                 "1419"="grey",
                                 "1421"="grey",
                                 "1425"="grey",
                                 "1839"="grey",
                                 "1841"="grey",
                                 "1842"="grey",
                                 "1843"="grey",
                                 "1844"="grey",
                                 "1845"="grey",
                                 "1846"="grey",
                                 "1847"="grey",
                                 "1848"="grey",
                                 "1849"="grey",
                                 "1850"="grey",
                                 "1851"="grey",
                                 "1852"="grey",
                                 "1853"="grey",
                                 "1854"="grey",
                                 "1856"="grey",
                                 "1438"="grey",
                                 "1439"="grey",
                                 "1443"="grey",
                                 "1445"="grey",
                                 "1456"="grey",
                                 "1466"="grey",
                                 "pgRNA-nontarget"="blue4",
                                 "1338"="grey",
                                 "1339"="grey",
                                 "1340"="grey",
                                 "1341"="grey",
                                 "1342"="grey",
                                 "1343"="grey",
                                 "1350"="grey",
                                 "1351"="grey",
                                 "1352"="grey",
                                 "1353"="grey",
                                 "1354"="grey",
                                 "1355"="grey",
                                 "1362"="grey",
                                 "1363"="grey",
                                 "1872"="grey",
                                 "1875"="grey",
                                 "1876"="grey",
                                 "1879"="grey",
                                 "1883"="grey",
                                 "1886"="grey",
                                 "1891"="grey")) +
  labs(x = "Depletion from Pool (%)", y = "Max Clonal Growth Inhibition (%)")+
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_rect("transparent"),
        panel.background = element_blank() )+
  theme(axis.title.x = element_text(color="black", size=15),
        axis.title.y = element_text(color="black", size=15))+
  theme(axis.text=element_text(size=15), legend.position = "none")+
  theme(legend.text = element_text(size = 12, colour = "black"))+
  theme(legend.title = element_text(size= 12, color = "black"))+
  geom_text(x=17, y=100, label="R = 0.42, p<2.2e-16", color="black", fontface="italic")+
  geom_text(x=70, y=94, label=" -- Control", color="blue4")+
  geom_text(x=75, y=100, label="-- EG Mutants", color="grey")
print(b)
b + ggsave("/Users/zisan/Files/CRISPRi-Seq/Rhamnose IC50_clonal/Rhamnose IC Values/240216_Depletion_Vs_clonal_growth_inhibition_240224.tiff", units="in", width=4.5, height=4.5, dpi=500, compression = 'lzw')

#CoD.lm = lm(mumax ~ CRISPRRater_Score_real, data=data)
#summary(CoD.lm)$r.squared 
#x.expression <- expression("(+)"~Rhamnose~µ[max])
#y.expression <- expression(?[max]~ "/"~ CRISPRRater)
#DepletionVsRank

q<-ggplot(data, aes(x = Concentration, y = OD, group = sgRNA, color = sgRNA)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("1828"="grey",
                             "1831"="grey",
                             "1832"="grey",
                             "1834"="grey",
                             "1835"="grey",
                             "1836"="grey",
                             "1837"="grey",
                             "1838"="grey",
                             "1404"="grey",
                             "1405"="grey",
                             "1406"="grey",
                             "1407"="grey",
                             "1408"="grey",
                             "1412"="grey",
                             "1413"="grey",
                             "1419"="grey",
                             "1421"="grey",
                             "1425"="grey",
                             "1839"="grey",
                             "1841"="grey",
                             "1842"="grey",
                             "1843"="grey",
                             "1844"="grey",
                             "1845"="grey",
                             "1846"="grey",
                             "1847"="grey",
                             "1848"="grey",
                             "1849"="grey",
                             "1850"="grey",
                             "1851"="grey",
                             "1852"="grey",
                             "1853"="grey",
                             "1854"="grey",
                             "1856"="grey",
                             "1438"="grey",
                             "1439"="grey",
                             "1443"="grey",
                             "1445"="grey",
                             "1456"="grey",
                             "1466"="grey",
                             "pgRNA-nontarget"="blue4",
                             "1338"="grey",
                             "1339"="grey",
                             "1340"="grey",
                             "1341"="grey",
                             "1342"="grey",
                             "1343"="grey",
                             "1350"="grey",
                             "1351"="grey",
                             "1352"="grey",
                             "1353"="grey",
                             "1354"="grey",
                             "1355"="grey",
                             "1362"="grey",
                             "1363"="grey",
                             "1872"="grey",
                             "1875"="grey",
                             "1876"="grey",
                             "1879"="grey",
                             "1883"="grey",
                             "1886"="grey",
                             "1891"="grey")) +
  labs(x = "Rhamnose Concentration (%)", y = "Percent Growth")+
  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_rect("transparent"),
        panel.background = element_blank() )+
  theme(axis.title.x = element_text(color="black", size=15),
        axis.title.y = element_text(color="black", size=15))+
  theme(axis.text=element_text(size=15), legend.position = "none")+
  theme(legend.text = element_text(size = 12, colour = "black"))+
  theme(legend.title = element_text(size= 12, color = "black"))+
  ylim(0,100)+
  geom_text(x=0.05, y=6, label=" -- Control", color="blue4")+
  geom_text(x=0.1, y=1, label="-- EG Mutants", color="grey")
print(q)

q+ggsave("/Users/zisan/Files/CRISPRi-Seq/Rhamnose IC50_clonal/Rhamnose IC Values/240216_Rha_dose_response.tiff", units="in", width=4.5, height=4.5, dpi=500, compression = 'lzw')


#q + ggsave("221130_Depletion_Vs_plusrahmnose_mumax.tiff", units="in", width=4.5, height=4.5, dpi=500, compression = 'lzw')

#Connect
# Marginal densities along x axis
xdens <- axis_canvas(a, axis = "x") +
  geom_density(data = data, aes(x = mumax_plusRha, fill = 'grey'),
               alpha = 1, size = 0.2) +
  scale_color_manual(values=c("grey"))+
  scale_fill_manual(values=c("grey"))
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(a, axis = "y", coord_flip = F) +
  geom_density(data = data, aes(y = Dep, fill = "darkred"),
               alpha = 1, size = 0.2) +
  scale_color_manual(values=c("grey"))+
  scale_fill_manual(values=c("grey"))
p1 <- insert_xaxis_grob(a, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)
p2 + ggsave("221204_Depletion_Vs_plusrahmnose_mumax_withdensity.tiff", units="in", width=5, height=5, dpi=500, compression = 'lzw')
