# Prelim:
library(JM)

#############################
# Data for multi-state model:
dtta <- aids
# Living states by CD4 groups:
CD4rnd <- round(dtta$CD4)
state  <- ifelse(CD4rnd<5,1,NA)
state  <- ifelse(5 <= CD4rnd & CD4rnd < 10,2,state)
state  <- ifelse(10 <= CD4rnd & CD4rnd < 15,3,state)
state  <- ifelse(15 <= CD4rnd ,4,state)
# Orientation:
# A higher CD 4 number indicates a stronger immune system.
# Make higher states less healthy:
state <- 5 - state
dtta$state <- state

# Other vars:
dtta$id    <- dtta$patient
dtta$drug  <- as.numeric(dtta$drug=="ddI")
dtta$prevOI <- as.numeric(dtta$prevOI == "AIDS")
dtta$gender <- as.numeric(dtta$gender == "male")

# Death state:
D <- 5
# Censored state:
censored <- -2
# Select variables:
dtta <- dtta[,c("id","state","obstime","CD4","drug", "death", "Time", "prevOI", "gender")]

# Add death or censoring:
subject <- unique(dtta$id)
N <- length(subject)
for(i in 1:N){
  dtta.i <- dtta[dtta$id==subject[i],]
  death <- dtta.i$death[1] == 1
  last.state  <- ifelse(death,D,censored)
  last.time   <- dtta.i$Time[1]
  last.record <- c(subject[i],last.state,last.time,NA,NA,
                   as.numeric(death),last.time,NA,NA)
  dtta.i <- rbind(dtta.i,last.record)
  if(i ==1){ddtta <- dtta.i}else{ddtta <- rbind(ddtta,dtta.i)}
}
dta_cov <- ddtta

