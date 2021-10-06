library(msm)

#Survival curve for patient administered drug 1 and with a previous AIDS diagnosis

#Setting covariates
covariates <- list( list(drug=1, prevOI=1), list(drug=1, prevOI=1), list(drug=1, prevOI=1),
                    list(drug=1, prevOI=1), list(drug=1, prevOI=1) )

#Setting grid for approximation
times <- c(0 ,5, 10, 15)

state1 <- c()
state2 <- c()
state3 <- c()
state4 <- c()

#Survival forecasted for 15 years 
for (i in 1:15) {
  m <- pmatrix.piecewise.msm(final_model, 0, i, times, covariates)
  state1[i] <- 1 - m[1, 5]
  state2[i] <- 1 - m[2, 5]
  state3[i] <- 1 - m[3, 5]
  state4[i] <- 1 - m[4, 5]
  
}


#Plotting survival curve
pdf("../project_report/figs/probability_plot_updated.pdf", width = 7, height = 4.5)
x <- c(1:15)
plot(x, state1, type='l', xlim = c(1, 15), ylim = c(0.3,1), ylab="Fitted Survival Probability",
     xlab = "Time (months)", col="red", lwd=1.8)
lines(state2, col="blue4", lwd=1.8)
lines(state3, col="seagreen", lwd=1.8)

lines(state4, col="orange1", lwd=1.8)
legend(1, 0.6, legend=c("From State 1", "From State 2", "From State 3", "From State 4"),
       col=c("red", "blue4", "seagreen", "orange1"), cex=0.8, lty=1, bty = "n")
dev.off()

#####################
#The string "mean" denotes the means of the covariates in the data 

covariates <- list( "mean", "mean", "mean",
                    "mean", "mean", "mean", "mean")
times <- c(0 ,3, 6, 9, 12, 18 )
for (i in 1:20) {
  m <- pmatrix.piecewise.msm(final_model, 0, i, times, covariates)
  state1[i] <- 1 - m[1,5]
  state2[i] <- 1 - m[2,5]
  state3[i] <- 1 - m[3,5]
  state4[i] <- 1 - m[4,5]
  
}

pdf("../project_report/figs/km_curves.pdf", width = 8, height = 5.5)
par(mfrow = c(2, 2),mar=c(4, 4, 2, 2) )

#Plotting Kaplan Meir curves and overlaying multi-state model survival curve 

#Patients with base state 1
tes1 <- dta_cov %>% 
  dplyr::filter(obstime == 0.00) %>% dplyr::select(-gender) %>%
  dplyr::filter(state == 1) 
#Plotting KM survival curve
plot(survfit(Surv(Time, death) ~ 1, data = tes1), xlim = c(0,20), ylim=c(0.3,1),
     ylab = "Survival Probability", xlab = "Time (months)")
#Plotting multi-state model survival
lines(state1, col="red")
legend(0.5, 0.5, legend=c("From State 1"),
       col=c("red"), cex=0.9, lty=1, bty = "n")

#Patients with base state 2
tes2 <- dta_cov %>% 
  dplyr::filter(obstime == 0.00) %>% dplyr::select(-gender) %>%
  dplyr::filter(state == 2) 
#Plotting KM survival curve
plot(survfit(Surv(Time, death) ~ 1, data = tes2), xlim = c(0,20), ylim=c(0.3,1),
     ylab = "Survival Probability", xlab = "Time (months)")
#Plotting multi-state model survival
lines(state2, col="red")
legend(0.5, 0.5, legend=c("From State 2"),
       col=c("red"), cex=0.9, lty=1, bty = "n")

#Patients with base state 3
tes3 <- dta_cov %>% 
  dplyr::filter(obstime == 0.00) %>% dplyr::select(-gender) %>%
  dplyr::filter(state == 3) 
#Plotting KM survival curve
plot(survfit(Surv(Time, death) ~ 1, data = tes3), xlim = c(0,20), ylim=c(0.3,1),
     ylab = "Survival Probability", xlab = "Time (months)")
#Plotting multi-state model survival
lines(state3, col="red")
legend(0.5, 0.5, legend=c("From State 3"),
       col=c("red"), cex=0.9, lty=1, bty = "n")

#Patients with base state 4
tes4 <- dta_cov %>% 
  dplyr::filter(obstime == 0.00) %>% dplyr::select(-gender) %>%
  dplyr::filter(state == 4) 
#Plotting KM survival curve
plot(survfit(Surv(Time, death) ~ 1, data = tes4), xlim = c(0,20), ylim=c(0.3,1),
     ylab = "Survival Probability", xlab = "Time (months)")
#Plotting multi-state model survival
lines(state4,col="red")
legend(0.2, 0.5, legend=c("From State 4"),
       col=c("red"), cex=0.9, lty=1, bty = "n")

dev.off()
