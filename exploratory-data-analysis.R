library(tidyverse)
library("JM")
library("lattice")
library(reshape2)
library(gridExtra)

aids_data <- aids

#Grouping patients by drug type
aids.id %>% 
  group_by(drug) %>% 
  count()

#Number of total deaths
aids %>%
  group_by(patient) %>%
  slice(tail(row_number(), 1)) %>%
  filter(death==1) %>% 
  nrow #188 deaths 

#Death by drug type
aids.id %>% 
  group_by(death, drug) %>% 
  summarise(n = n()) %>%  
  mutate(freq = prop.table(n)*100)

#Frequency of states recorded
ddtta %>%
  group_by(id) %>% 
  group_by(state) %>%
  summarise(n = n()) %>% 
  mutate(Proportion= prop.table(n)*100)

#Graph showing the observed states at various check up times
plot <- ggplot(ddta, aes(x = obstime, fill = factor(state))) +
  geom_histogram(bins = 25, alpha = 0.6, position = "stack") +
  scale_x_continuous(breaks = seq(0, 25, by = 2)) +
  scale_fill_manual(values = c("#808080", "#FEF001", 
                               "#FFCE03", "#FD9A01", "#FD6104", "#F00505")) +
  labs(fill="State", 
       x = "Check up time (months)",
       y = "Count",) +
  theme(text=element_text(size=10)) + theme_minimal()

ggsave("../project_report/figs/plot1.pdf",plot, device=cairo_pdf, width=7, height=4.5, units="in")

#Graph showing examples of different hazards
x <- seq(0, 20, length.out=10000)
exp1 <- rep(1, 10000)
gomp1 <- 0.5*exp(0.1*x)
weib <- 0.5 * 1.5 * x^(1.5-1)
lamda <- 1.5
rho <- 3
loglog <- (lamda * rho * (lamda*x)^(rho-1))/(1+(lamda*x)^(rho))

df. <- data.frame(x, exp1, gomp1, weib, loglog)
colnames(df.) <- c("x" ,"Exp(1), \u03BB = 1", "Gompertz (0.5, 0.1), \u03BB = 0.5 \u03B6=0.1", "Weibull (0.5, 1.5), \u03BB = 0.5 \u03C4=1.5", "Log-logistic (1.5, 3), \u03BB = 1.5 \u03C1 = 3" )
df2. <- melt(data = df., id.vars = "x")
colnames(df2.) <- c("Time", "variable", "Hazard")

intensities <- ggplot(data = df2., aes(x = Time, y = Hazard, colour = variable)) + geom_line(size=1) +
  ylim(0, 4) +
  theme_minimal() +
  theme(legend.position=c(0.2,0.85), legend.title = element_blank(), legend.text = element_text(size=12)) 

ggsave("../project_report/figs/intensity_graph.pdf",intensities, device=cairo_pdf, width=20, height=15, units="cm")