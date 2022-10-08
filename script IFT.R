################################################################################
################################################################################
#########################                        ###############################
#########################  SIS stochastic model  ###############################
#########################                        ###############################
################################################################################
################################################################################

library(ggplot2)

#Closed-form solution to the mean-field equation

logistic_growth <- function(alpha, beta, n_pop, I_0, time_step){
  K = n_pop * (1 - (beta/alpha))
    return (K/(1 + ((K/I_0) - 1) * exp(-(alpha - beta) * time_step)))
}

#Population size

n_pop = 1000

#Initial number of infectives (minimum 1)

I_0 = 900

#Time step

delta_T = 1

#Parameter which dictates the force of infection

alpha = 0.0005

#Parameter that dictates the recovery rate

beta = 0.001

#Maximum number of iterations per realization

iter_max = 2500

#Number of realizations

realizations = 50

#Vector containing the values of the closed-form solution to the mean-field approximation

total_expectation <- vector("numeric", iter_max)

for(i in 1:iter_max) total_expectation[i] <- ((100 * logistic_growth(alpha, beta, n_pop, I_0, i * delta_T))/n_pop)

#for(i in 1:iter_max) total_expectation[i] <- (100/(alpha * i *delta_T + (n_pop/I_0)))

mean_infectives <- vector("numeric", iter_max)
for(i in 1:iter_max) mean_infectives[i] <- 0

sd_infectives <- vector("numeric", iter_max)
for(i in 1:iter_max) sd_infectives[i] <- 0

infectives = data.frame(matrix(nrow = realizations, ncol = iter_max)) 

for(i in 1:realizations) infectives[i, 1] <- I_0

for(i in 1:realizations){
  for(j in 2:iter_max){
    p <- runif(n = 1, min =  0, max = 1)
    if(p <= beta * delta_T * infectives[i, j - 1]){
      infectives[i, j] <- infectives[i, (j - 1)] - 1
    }
    else{
      if(p <= (beta + alpha * (1 - infectives[i, (j - 1)]/n_pop)) * delta_T * infectives[i, (j - 1)]){
        infectives[i, j] <- infectives[i, (j - 1)] + 1
      }
      else{
        infectives[i, j] <- infectives[i, (j - 1)]
      }
    }
  }
}

#Data visualization
#Computes the mean and standard deviation between all realizations

for(j in 1:iter_max){
  mean_infectives[j] = 100 * mean(infectives[, j])/n_pop
  sd_infectives[j] = 100 * sd(infectives[, j])/n_pop
}

x <- 1:(iter_max/100)

plot_mean_infectives <- vector("numeric", iter_max/100)
for(i in x) plot_mean_infectives[i] <- mean_infectives[i * 100]

plot_sd_infectives <- vector("numeric", iter_max/100)
for(i in x) plot_sd_infectives[i] <- sd_infectives[i * 100]

#Data visualization - scatterplot with few points and error bars

plot(x, plot_mean_infectives,
     ylim=range(c(plot_mean_infectives - plot_sd_infectives, plot_mean_infectives + plot_sd_infectives)),
     pch=19, xlab="Time", ylab="Percentage (%)",
     main = "Epidemic evolution"
)

arrows(x, plot_mean_infectives - plot_sd_infectives, x, plot_mean_infectives + plot_sd_infectives, length = 0.05, angle = 90, code = 3)

#Data visualization - scatterplot with lots of points, no error bars and with fitting curve

z  <- 1:iter_max
y1 <- mean_infectives
y2 <- total_expectation
df <- data.frame(z, y1, y2)

ggplot(df, aes(z)) + geom_line(aes(y = y1), colour="red", size = 0.5) + 
  geom_line(aes(y = y2), colour="green", linetype = "dashed") + xlab("Time") + ylab("Percentage") + 
  ggtitle("Epidemic evolution")

#Computes the fraction of realizations entering the absorbing state (0) in the time
#window of the fixed number of iterations

count = 0

for(i in 1:realizations){
  if(infectives[i, iter_max] == 0) count = count + 1
}

count/realizations
