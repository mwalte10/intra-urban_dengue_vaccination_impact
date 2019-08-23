############################
#CRC info
############################
args = commandArgs(TRUE)
input = as.numeric(args[1])


##########################
#IF DOING VAC COVERAGE
##########################

validation_parms.test <- c(0.1660247, 0.1717892, 0.1826108 , 0.2045885, 0.2716160)
beta_h <- validation_parms.test[input]
beta_l <- beta_h
native <- 1
travel <- 0
vac_h <- 0.8
vac_l <- 0.8

load('pop_1950.RData')
load('birth_1950.RData')
load('death_1950.RData')


hopkins <- c(0.53, 1, 0.115)
hopkins_inverse <- 1 - hopkins


library(deSolve)




############################
#Initial conidtions and parameters
###########################

initial_conditions <- as.data.frame(matrix(NA, nrow = 80*4, ncol = 3))
initial_conditions[,2] <- rep(0:79,4)
for(i in 1:4){
  x <- i - 1
  initial_conditions[x*(80) + (1:80),1] <- rep(i,80)
}
initial_conditions[,3] <- rep(pop,4)

susceptible_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
susceptible_init_h[,3] <- c(initial_conditions[1:80,3] * 0.5, initial_conditions[81:160,3] * 0, 
                            initial_conditions[161:240,3] * 0, initial_conditions[241:320,3] * 0)
susceptible_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
susceptible_init_l[,3] <- c(initial_conditions[1:80,3] * 0.5, initial_conditions[81:160,3] * 0, 
                            initial_conditions[161:240,3] * 0, initial_conditions[241:320,3] * 0)

infected_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
infected_init_h[,3] <- c(rep(1/80, 80), initial_conditions[81:160,3] * 0, 
                         initial_conditions[161:240,3] * 0, initial_conditions[241:320,3] * 0)
infected_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
infected_init_l[,3] <- c(rep(1/80, 80), initial_conditions[81:160,3] * 0, 
                         initial_conditions[161:240,3] * 0, initial_conditions[241:320,3] * 0)

recovered_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
recovered_init_h[,3] <- c(initial_conditions[1:80,3] * 0, initial_conditions[81:160,3] * 0, 
                          initial_conditions[161:240,3] * 0, initial_conditions[241:320,3] * 0)
recovered_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
recovered_init_l[,3] <- c(initial_conditions[1:80,3] * 0, initial_conditions[81:160,3] * 0, 
                          initial_conditions[161:240,3] * 0, initial_conditions[241:320,3] * 0)

susceptible_h <- function(exposure){
  x <- exposure - 1
  susceptible <- c(susceptible_init_h[(x * 80) + 1:80,3])
}
susceptible_total_h <- sum(susceptible_h(1) + susceptible_h(2) + susceptible_h(3) + susceptible_h(4))
susceptible_l <- function(exposure){
  x <- exposure - 1
  susceptible <- c(susceptible_init_l[(x * 80) + 1:80,3])
}
susceptible_total_l <- sum(susceptible_l(1) + susceptible_l(2) + susceptible_l(3) + susceptible_l(4))

infected_h <- function(exposure){
  x <- exposure - 1
  infected <- c(infected_init_h[(x * 80) + 1:80,3])
} 
infected_total_h <- sum(infected_h(1) + infected_h(2) + infected_h(3) + infected_h(4))
infected_l <- function(exposure){
  x <- exposure - 1
  infected <- c(infected_init_l[(x * 80) + 1:80,3])
} 
infected_total_l <- sum(infected_l(1) + infected_l(2) + infected_l(3) + infected_l(4))

recovered_h <- function(exposure){
  x <- exposure - 1
  recovered <- c(recovered_init_h[(x * 80) + 1:80,3])
} 
recovered_total_h <- sum(recovered_h(1) + recovered_h(2) + recovered_h(3) + recovered_h(4))
recovered_l <- function(exposure){
  x <- exposure - 1
  recovered <- c(recovered_init_l[(x * 80) + 1:80,3])
} 
recovered_total_l <- sum(recovered_l(1) + recovered_l(2) + recovered_l(3) + recovered_l(4))

population_h <- sum(susceptible_total_h + infected_total_h + recovered_total_h)
population_l <- sum(susceptible_total_l + infected_total_l + recovered_total_l)



parms.h <- list(beta_h = beta_h,
                beta_l = beta_l,
                gamma = 1/4,
                sigma = 1/(365 * 1.2),
                mu = birth,
                delta = death,
                age_window = rep(1, 80),
                native = native,
                travel = 1 - native,
                vac_h = vac_h,
                vac_l = vac_l,
                sens = 1,
                spec = 0,
                hopkins,
                hopkins_inverse)
parms_null.h <- list(beta_h = beta_h,
                     beta_l = beta_l,
                     gamma = 1/4,
                     sigma = 1/(365 * 1.2),
                     mu = birth,
                     delta = death,
                     age_window = rep(1, 80),
                     native = native,
                     travel = 1 - native,
                     vac_h = 0,
                     vac_l = 0,
                     sens = 1,
                     spec = 0,
                     hopkins,
                     hopkins_inverse)



years = 60
years_vac = 30
times <- seq(from = 0, to = 365 * years, by = 0.1)
times <- times[1:(length(times) - 1)]

############################
#MODEL
############################
model <- function(t, y, parms){
  
  beta_h <- parms[[1]]
  beta_l <- parms[[2]]
  gamma <- parms[[3]]
  sigma <- parms[[4]]
  mu <- parms[[5]][floor(t / 365) + 1]
  delta <- parms[[6]][[floor(t / 365) + 1]]
  age_window <- parms[[7]]
  native <- parms[[8]]
  travel <- parms[[9]]
  t_vac.h <- ifelse((t>(365*(years_vac))), parms[[10]] / 365, 0)
  t_vac.l <- ifelse((t>(365*(years_vac))), parms[[11]] / 365, 0)
  vac_h <- c(rep(0,8), t_vac.h, rep(0,71))
  vac_l <- c(rep(0,8), t_vac.l, rep(0,71))
  sens <- parms[[12]]
  spec <- parms[[13]]
  inf <- parms[[14]]
  ninf <- parms[[15]]
  
  #Gen pop
  S1_h <- y[which(names(y)== 'sh1')]
  I1_h <- y[which(names(y)== 'ih1')]
  R1_h <- y[which(names(y)== 'rh1')]
  S2_h <- y[which(names(y)== 'sh2')]
  I2_h <- y[which(names(y)== 'ih2')]
  R2_h <- y[which(names(y)== 'rh2')]
  S3_h <- y[which(names(y)== 'sh3')]
  I3_h <- y[which(names(y)== 'ih3')]
  R3_h <- y[which(names(y)== 'rh3')]
  S4_h <- y[which(names(y)== 'sh4')]
  I4_h <- y[which(names(y)== 'ih4')]
  R4_h <- y[which(names(y)== 'rh4')]
  
  #Vaccinated
  R1_h.v <- y[which(names(y)== 'rh1.v')]
  S2_h.v <- y[which(names(y)== 'sh2.v')]
  I2_h.v <- y[which(names(y)== 'ih2.v')]
  R2_h.v <- y[which(names(y)== 'rh2.v')]
  S3_h.v <- y[which(names(y)== 'sh3.v')]
  I3_h.v <- y[which(names(y)== 'ih3.v')]
  R3_h.v <- y[which(names(y)== 'rh3.v')]
  S4_h.v <- y[which(names(y)== 'sh4.v')]
  I4_h.v <- y[which(names(y)== 'ih4.v')]
  R4_h.v <- y[which(names(y)== 'rh4.v')]
  
  #Gen pop
  S1_l <- y[which(names(y)== 'sl1')]
  I1_l <- y[which(names(y)== 'il1')]
  R1_l <- y[which(names(y)== 'rl1')]
  S2_l <- y[which(names(y)== 'sl2')]
  I2_l <- y[which(names(y)== 'il2')]
  R2_l <- y[which(names(y)== 'rl2')]
  S3_l <- y[which(names(y)== 'sl3')]
  I3_l <- y[which(names(y)== 'il3')]
  R3_l <- y[which(names(y)== 'rl3')]
  S4_l <- y[which(names(y)== 'sl4')]
  I4_l <- y[which(names(y)== 'il4')]
  R4_l <- y[which(names(y)== 'rl4')]
  
  #Vaccinated
  R1_l.v <- y[which(names(y)== 'rl1.v')]
  S2_l.v <- y[which(names(y)== 'sl2.v')]
  I2_l.v <- y[which(names(y)== 'il2.v')]
  R2_l.v <- y[which(names(y)== 'rl2.v')]
  S3_l.v <- y[which(names(y)== 'sl3.v')]
  I3_l.v <- y[which(names(y)== 'il3.v')]
  R3_l.v <- y[which(names(y)== 'rl3.v')]
  S4_l.v <- y[which(names(y)== 'sl4.v')]
  I4_l.v <- y[which(names(y)== 'il4.v')]
  R4_l.v <- y[which(names(y)== 'rl4.v')]
  
  ###Tracking infected individuals
  
  
  infecteds_h.prim.c <- I1_h * inf[1]
  infecteds_h.sec.c <- cbind(I2_h, I2_h.v) * inf[2]
  infecteds_h.psec.c <- cbind(I3_h, I4_h, I3_h.v, I4_h.v) * inf[3]
  sym_inf_h <- cbind(infecteds_h.prim.c, infecteds_h.sec.c, infecteds_h.psec.c)
  sym_inf_h <- rowSums(sym_inf_h)
  
  infecteds_h.prim <- I1_h * ninf[1]
  infecteds_h.sec <- cbind(I2_h, I2_h.v) * ninf[2]
  infecteds_h.psec <- cbind(I3_h, I4_h, I3_h.v, I4_h.v) * ninf[3]
  inf_h <- cbind(infecteds_h.prim, infecteds_h.sec, infecteds_h.psec)
  inf_h <- rowSums(inf_h)
  
  
  population_h <- cbind(S1_h, I1_h, R1_h,
                        S2_h, I2_h, R2_h,
                        S3_h, I3_h, R3_h,
                        S4_h, I4_h, R4_h,
                        R1_h.v,
                        S2_h.v, I2_h.v, R2_h.v,
                        S3_h.v, I3_h.v, R3_h.v,
                        S4_h.v, I4_h.v, R4_h.v)
  pop_h <- rowSums(population_h)
  birth_pop_h <- sum(pop_h)
  
  infecteds_l.prim.c <- I1_l * inf[1]
  infecteds_l.sec.c <- cbind(I2_l, I2_l.v) * inf[2]
  #infecteds_l.sec.c <- rowSums(infecteds_l.sec.c)
  infecteds_l.psec.c <- cbind(I3_l, I4_l, I3_l.v, I4_l.v) * inf[3]
  #infecteds_l.psec.c <- rowSums(infecteds_l.psec.c)
  sym_inf_l <- cbind(infecteds_l.prim.c, infecteds_l.sec.c, infecteds_l.psec.c)
  sym_inf_l <- rowSums(sym_inf_l)
  
  infecteds_l.prim <- I1_l * ninf[1]
  infecteds_l.sec <- cbind(I2_l, I2_l.v) * ninf[2]
  #infecteds_l.sec <- rowSums(infecteds_l.sec.c)
  infecteds_l.psec <- cbind(I3_l, I4_l, I3_l.v, I4_l.v) * ninf[3]
  #infecteds_l.psec <- rowSums(infecteds_l.psec)
  inf_l <- cbind(infecteds_l.prim, infecteds_l.sec, infecteds_l.psec)
  inf_l <- rowSums(inf_l)
  
  
  population_l <- cbind(S1_l, I1_l, R1_l,
                        S2_l, I2_l, R2_l,
                        S3_l, I3_l, R3_l,
                        S4_l, I4_l, R4_l,
                        R1_l.v,
                        S2_l.v, I2_l.v, R2_l.v,
                        S3_l.v, I3_l.v, R3_l.v,
                        S4_l.v, I4_l.v, R4_l.v)
  pop_l <- rowSums(population_l)
  birth_pop_l <- sum(pop_l)
  # effective_population_h <- population_h + rep(travel,22) * population_l
  # effective_population_h <- sum(effective_population_h)
  # effective_population_l <- population_l + rep(travel,22) * population_h
  # effective_population_l <- sum(effective_population_l)
  
  #first infection
  dS1_h <-  
    mu * c(birth_pop_h, rep(0,79)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_h, -1), 0) -
    native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S1_h * delta -
    S1_h * vac_h * (1 - spec)
  dS1_l <-  
    mu * c(birth_pop_l, rep(0,79)) +
    c(0, 1/365 / head(age_window, -1) * head(S1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S1_l, -1), 0)  -
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) -
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) -
    S1_l * delta - 
    S1_l * vac_l * (1 - spec)
  
  dI1_h <- 
    native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I1_h, -1), 0) -
    I1_h * gamma -
    I1_h * delta
  dI1_l <- 
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I1_l, -1), 0) -
    I1_l * gamma -
    I1_l * delta
  
  dR1_h <-
    I1_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_h, -1), 0) -
    R1_h * sigma -
    R1_h * delta - 
    R1_h * vac_h * sens
  dR1_h.v <- 
    S1_h * vac_h * (1 - spec) +
    c(0, 1/365 / head(age_window, -1) * head(R1_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_h.v, -1), 0) -
    R1_h.v * delta -
    R1_h.v * sigma
  dR1_l <-
    I1_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R1_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_l, -1), 0) -
    R1_l * sigma -
    R1_l * delta - 
    R1_l * vac_l * sens
  dR1_l.v <- 
    S1_l * vac_l * (1 - spec) +
    c(0, 1/365 / head(age_window, -1) * head(R1_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R1_l.v, -1), 0) -
    R1_l.v * delta -
    R1_l.v * sigma
  
  #second infection
  dS2_h <- 
    R1_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h, -1), 0) -
    native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S2_h * delta - 
    S2_h * vac_h * sens
  dS2_h.v <-
    R1_h.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_h.v, -1), 0) -
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S2_h.v * delta
  dS2_l <- 
    R1_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l, -1), 0) -
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S2_l * delta - 
    S2_l * vac_l * sens
  dS2_l.v <-
    R1_l.v * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S2_l.v, -1), 0) -
    native * S2_l.v * beta_l * 2 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S2_l.v * delta
  
  dI2_h <- 
    native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h, -1), 0) -
    I2_h * gamma -
    I2_h * delta
  dI2_h.v <-
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_h.v, -1), 0) -
    I2_h.v * gamma -
    I2_h.v * delta
  dI2_l <- 
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_l, -1), 0) -
    I2_l * gamma -
    I2_l * delta
  dI2_l.v <-
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I2_l.v, -1), 0) -
    I2_l.v * gamma -
    I2_l.v * delta
  
  dR2_h <-  
    I2_h * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_h, -1), 0) -
    R2_h * sigma -
    R2_h * delta -
    R2_h * vac_h * sens
  dR2_h.v <-
    I2_h.v * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_h.v, -1), 0) -
    R2_h.v * sigma -
    R2_h.v * delta +
    R1_h * vac_h * sens +
    S2_h * vac_h * sens
  dR2_l <-  
    I2_l * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_l, -1), 0) -
    R2_l * sigma -
    R2_l * delta -
    R2_l * vac_l * sens
  dR2_l.v <-
    I2_l.v * gamma +  
    c(0, 1/365 / head(age_window, -1) * head(R2_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R2_l.v, -1), 0) -
    R2_l.v * sigma -
    R2_l.v * delta +
    R1_l * vac_l * sens +
    S2_l * vac_l * sens
  
  #third infection
  dS3_h <- 
    R2_h * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h, -1), 0) -
    native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S3_h * delta -
    S3_h * vac_h * sens
  dS3_h.v <-
    R2_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_h.v, -1), 0) -
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S3_h.v * delta 
  dS3_l <- 
    R2_l * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l, -1), 0) -
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S3_l * delta -
    S3_l * vac_l * sens
  dS3_l.v <-
    R2_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S3_l.v, -1), 0) -
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S3_l.v * delta
  
  dI3_h <- 
    native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h, -1), 0) -
    I3_h * gamma -
    I3_h * delta
  dI3_h.v <-
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_h.v, -1), 0) -
    I3_h.v * gamma -
    I3_h.v * delta
  dI3_l <- 
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_l, -1), 0) -
    I3_l * gamma -
    I3_l * delta
  dI3_l.v <-
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I3_l.v, -1), 0) -
    I3_l.v * gamma -
    I3_l.v * delta
  
  dR3_h <- 
    I3_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_h, -1), 0) -
    R3_h * sigma -
    R3_h * delta -
    R3_h * vac_h * sens
  dR3_h.v <-
    I3_h.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_h.v, -1), 0) -
    R3_h.v * sigma + 
    R2_h * vac_h * sens +
    S3_h * vac_h * sens - 
    R3_h.v * delta
  dR3_l <- 
    I3_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_l, -1), 0) -
    R3_l * sigma -
    R3_l * delta - 
    R3_l * vac_l * sens
  dR3_l.v <-
    I3_l.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R3_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R3_l.v, -1), 0) -
    R3_l.v * sigma + 
    R2_l * vac_l * sens +
    S3_l * vac_l * sens -
    R3_l.v * delta
  
  #fourth infection
  dS4_h <- 
    R3_h * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h, -1), 0) -
    native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S4_h * delta - 
    S4_h * vac_h * sens
  dS4_h.v <-
    R3_h.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_h.v, -1), 0) -
    native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) -
    travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) -
    native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) -
    travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) -
    S4_h.v * delta 
  dS4_l <- 
    R3_l * sigma +
    c(0, 1/365 / head(age_window, -1) * head(S4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l, -1), 0) -
    native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S4_l * delta - 
    S4_l * vac_l * sens
  dS4_l.v <-
    R3_l.v * sigma + 
    c(0, 1/365 / head(age_window, -1) * head(S4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(S4_l.v, -1), 0) -
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) -
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) -
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) -
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) -
    S4_l.v * delta 
  
  dI4_h <- 
    native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h, -1), 0) -
    I4_h * gamma -
    I4_h * delta
  dI4_h.v <-
    native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    c(0, 1/365 / head(age_window, -1) * head(I4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_h.v, -1), 0) -
    I4_h.v * gamma -
    I4_h.v * delta
  dI4_l <- 
    native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_l, -1), 0) -
    I4_l * gamma -
    I4_l * delta
  dI4_l.v <-
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    c(0, 1/365 / head(age_window, -1) * head(I4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(I4_l.v, -1), 0) -
    I4_l.v * gamma -
    I4_l.v * delta
  
  dR4_h <- 
    I4_h * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_h, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_h, -1), 0) -
    R4_h * delta 
  dR4_h.v <-
    I4_h.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_h.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_h.v, -1), 0) -
    R4_h.v * delta + 
    R3_h * vac_h * sens +
    S4_h * vac_h * sens
  dR4_l <- 
    I4_l * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_l, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_l, -1), 0) -
    R4_l * delta 
  dR4_l.v <-
    I4_l.v * gamma +
    c(0, 1/365 / head(age_window, -1) * head(R4_l.v, -1)) -
    c(1/365 / head(age_window, -1) * head(R4_l.v, -1), 0) -
    R4_l.v * delta + 
    R3_l * vac_l * sens +
    S4_l * vac_l * sens
  
  I_tot <-
    native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) +
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  I_total <- sum(I_tot)
  
  I_h <- 
    native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
  I_h <- sum(I_h)
  
  I_primary_tot <- 
    native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) 
  I_primary_tot <- sum(I_primary_tot)
  
  I_secondary_tot <-
    native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  I_secondary_tot <- sum(I_secondary_tot)
  
  I_secondary_vac <- 
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)
  I_secondary_vac <- sum(I_secondary_vac)
  
  
  I_post_sec_tot <- 
    native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  I_post_sec_tot <- sum(I_post_sec_tot)
  
  I_post_sec_vac <- 
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  I_post_sec_vac <- sum(I_post_sec_vac)
  
  
  
  I_l <-
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) +
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  I_l <- sum(I_l)
  
  I_l_vac <- 
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  I_l_vac <- sum(I_l_vac)
  
  
  I_l_primary_tot <- 
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)
  I_l_primary_tot <- sum(I_l_primary_tot)
  
  I_l_sec_tot <-
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)
  I_l_sec_tot <- sum(I_l_sec_tot)
  
  I_l_sec_vac <- 
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)
  I_l_sec_vac <- sum(I_l_sec_vac)
  
  I_l_post_sec_tot <- 
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  I_l_post_sec_tot <- sum(I_l_post_sec_tot)
  
  I_l_post_sec_vac <- 
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  I_l_post_sec_vac <- sum(I_l_post_sec_vac)
  
  
  primary <-
    native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
    native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)
  primary <- sum(primary) * inf[1]
  primary_inf <- primary/ inf[1]
  
  
  
  ##primary cases incidence
  {  primary_tot.cases <- 
      native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)
    primary_tot.cases <- sum(primary_tot.cases) * inf[1]
    
    primary_tot.ncases <- 
      native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)
    primary_tot.ncases <- sum(primary_tot.ncases) * ninf[1]
    
    primary_cases.l <-
      native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) 
    primary_cases.l <- sum(primary_cases.l) * inf[1]
    
    primary_ncases.l <- 
      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)
    native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    primary_ncases.l <- sum(primary_ncases.l) * ninf[1]
    
    primary_cases.h <- 
      native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    primary_cases.h <- sum(primary_cases.h) * inf[1]
    
    primary_ncases.h <- 
      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) 
    primary_ncases.h <- sum(primary_ncases.h) * ninf[1]}
  
  
  secondary <-
    native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  secondary <- sum(secondary) * inf[2]
  secondary_inf <- secondary / inf[2]
  
  
  #vaccinated case incidence
  {  secondary_vac_tot.cases <- 
      native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
    secondary_vac_tot.cases <- sum(secondary_vac_tot.cases) * inf[2]
    
    secondary_vac_tot.ncases <- 
      native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)  
    secondary_vac_tot.ncases <- sum(secondary_vac_tot.ncases) * ninf[2]
    
 
    
    secondary_vac_ncases.l <- 
      native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)  
    secondary_vac_ncases.l <- sum(secondary_vac_ncases.l) * ninf[2]
    
    secondary_vac_cases.l <- 
      native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)  
    secondary_vac_cases.l <- sum(secondary_vac_cases.l) * inf[2]

   
     secondary_vac_cases.h <- 
      native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)
    secondary_vac_cases.h <- sum(secondary_vac_cases.h) * inf[2]
    
    secondary_vac_ncases.h <- 
      native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)
    secondary_vac_ncases.h <- sum(secondary_vac_ncases.h) * ninf[2]}
  
  #nonvaccinated case incidence
  {  secondary_tot.cases <-
      native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)
    secondary_tot.cases <- sum(secondary_tot.cases) * inf[2]
    
    secondary_tot.ncases <-
      native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
    secondary_tot.ncases <- sum(secondary_tot.ncases) * ninf[2]
    
    secondary_cases.l <-
      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    secondary_cases.l <- sum(secondary_cases.l) * inf[2]
    
    secondary_ncases.l <-
      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)
    native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    secondary_ncases.l <- sum(secondary_ncases.l) * ninf[2]
    
    secondary_cases.h <-
      native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    secondary_cases.h <- sum(secondary_cases.h) * inf[2]
    
    secondary_ncases.h <-
      native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    secondary_ncases.h <- sum(secondary_ncases.h) * ninf[2]}
  
  
  tertiary <-
    native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
    travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
    native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
    travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
    native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
    native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
    native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
    travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
  tertiary <- sum(tertiary) * inf[3]
  
  
  
  ##post sec vaccinated case incidence
  { postsec_vac_tot.cases <- 
      native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
    postsec_vac_tot.cases <- sum(postsec_vac_tot.cases) * inf[3]
    
    postsec_vac_tot.ncases <- 
      native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
    postsec_vac_tot.ncases <- sum(postsec_vac_tot.ncases) * ninf[3]
    
    postsec_vac_cases.l <- 
      travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h)
    postsec_vac_cases.l <- sum(postsec_vac_cases.l) * inf[3]
    
    postsec_vac_ncases.l <- 
      travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    postsec_vac_ncases.l <- sum(postsec_vac_ncases.l) * ninf[3]
    
    postsec_vac_cases.h <- 
      native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    postsec_vac_cases.h <- sum(postsec_vac_cases.h) * inf[3]
    
    postsec_vac_ncases.h <- 
      native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    postsec_vac_ncases.h <- sum(postsec_vac_ncases.h) * ninf[3]}
  
  ##post sec nonvaccinated case incidence
  {  postsec_tot_cases <- 
      native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
    postsec_tot_cases <- sum(postsec_tot_cases) * inf[3]
    
    postsec_tot_ncases <- 
      native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
    postsec_tot_ncases <- sum(postsec_tot_ncases) * ninf[3]
    
    postsec_cases.l <- 
      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    postsec_cases.l <- sum(postsec_cases.l) * inf[3]
    
    postsec_ncases.l <- 
      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    postsec_ncases.l <- sum(postsec_ncases.l) * ninf[3]
    
    postsec_cases.h <- 
      native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    postsec_cases.h <- sum(postsec_cases.h) * inf[3]
    
    postsec_ncases.h <- 
      native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    postsec_ncases.h <- sum(postsec_ncases.h) * ninf[3]
    
  }
  {  
    quaternary <-
      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) 
    quaternary <- sum(quaternary) * inf[3]
    
    cases <- primary + secondary + tertiary + quaternary
    post_sec <- tertiary + quaternary
    post_sec_inf <- post_sec / inf[3]
    
    primary.h <- 
      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) 
    primary.h <- sum(primary.h) * inf[1]
    
    secondary.h <- 
      native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    
    secondary.h <- sum(secondary.h) * inf[2]
    
    post_sec.h <- 
      native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    
    post_sec.h <- sum(post_sec.h) * inf[3]
    
    cases.h <- primary.h + secondary.h + post_sec.h
    
    primary.l <-
      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l) +
      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
      native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +    
      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) 
    primary.l <- sum(primary.l) * inf[1]
    primary.l.inf <- primary.l / inf[1]
    
    secondary.l <-
      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) 
    secondary.l <- sum(secondary.l) * inf[2]
    secondary.l.inf <- secondary.l / inf[2]
    
    tertiary.l <-
      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    tertiary.l <- sum(tertiary.l) * inf[3]
    
    quaternary.l <-
      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) + 
      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) + 
      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) 
    quaternary.l <- sum(quaternary.l) * inf[3]
    
    cases.l <- primary.l + secondary.l  + tertiary.l + quaternary.l
    post_sec.l <- tertiary.l + quaternary.l
    post_sec.l.inf <- post_sec.l / inf[3]
  }
  
  ##whole population
  {  ##s and h
    {  s <- 
      inf[1] * sum(native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)) +
      inf[2] * sum(native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                     native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
      inf[3] * sum(native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                     native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                     native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                     native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))
    
    h <- 
      ninf[1] * sum(native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)) +
      ninf[2] * sum(native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                      native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
      ninf[3] * sum(native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                      native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))}
    ##s.h and h.h
    {
      s.h <- 
        inf[1] * sum(native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h))  +
        
        inf[2] * sum(native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                       native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) ) +
        
        inf[3] * sum(native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                       native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                       native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                       native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) )
      
      h.h <- 
        ninf[1] * sum(native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                        travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                        travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)) +
        
        ninf[2] * sum(native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                        travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                        travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                        native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                        travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                        travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)) +
        
        ninf[3] * sum(native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                        travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                        travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                        native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                        travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                        travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                        native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                        travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                        travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                        native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                        travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                        travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) )
    }
    
    ##s.l and h.l
    {  s.l <- 
        inf[1] * sum( native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                        travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)) +
        inf[2] * sum(native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                       native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
        inf[3] * sum(native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                       native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                       native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                       native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))
      
      h.l <- 
        ninf[1] * sum(native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                        travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)) +
        ninf[2] * sum(native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                        native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
        ninf[3] * sum(native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                        native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                        native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                        native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))}
  }
  
  ##vaccinated population
  { ##s and h
    {  s.v <- 
      
      inf[2] * sum(native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
      inf[3] * sum(native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                     native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))
    
    h.v <- 
      ninf[2] * sum(native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
      ninf[3] * sum(native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                      native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))}
    ##s.h and h.h
    {
      s.h.v <- 
        
        inf[2] * sum(native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) ) +
        
        inf[3] * sum(native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                       native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                       travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                       native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                       travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) )
      
      h.h.v <- 
        ninf[2] * sum(native * S2_h.v * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                        travel * S2_h.v * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                        native * S2_h.v * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                        travel * S2_h.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)) +
        
        ninf[3] * sum( native * S3_h.v * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                         travel * S3_h.v * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                         native * S3_h.v * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                         travel * S3_h.v * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                         native * S4_h.v * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                         travel * S4_h.v * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                         native * S4_h.v * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                         travel * S4_h.v * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) )
    }
    
    ##s.l and h.l
    {  s.l.v <- 
        
        inf[2] * sum(native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
        inf[3] * sum(native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                       native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                       travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                       native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                       travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))
      
      h.l.v <- 
        ninf[2] * sum(native * S2_l.v * beta_l * 2  * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S2_l.v * beta_h * 2 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S2_l.v * beta_l * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S2_l.v * beta_h * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
        ninf[3] * sum(native * S3_l.v * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S3_l.v * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S3_l.v * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S3_l.v * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                        native * S4_l.v * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                        travel * S4_l.v * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                        native * S4_l.v * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                        travel * S4_l.v * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))}
  }
  
  ##non-vaccinated population
  ##s and h
  {  s.nv <- 
      inf[1] * sum(native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)) +
      inf[2] * sum(native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)) +
      inf[3] * sum(native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                     native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))
    
    h.nv <- 
      ninf[1] * sum(native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)) +
      ninf[2] * sum(native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)) +
      ninf[3] * sum(native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) +
                      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) )}
  ##s.h and h.h
  {
    s.h.nv <- 
      inf[1] * sum(native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h))  +
      
      inf[2] * sum(native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)  ) +
      
      inf[3] * sum(native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                     native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                     travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                     native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                     travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) )
    
    h.h.nv <- 
      ninf[1] * sum(native * S1_h * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S1_h * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S1_h * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S1_h * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)) +
      
      ninf[2] * sum(native * S2_h * 2 * 0.75 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S2_h * 2 * 0.75 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S2_h * 0.75 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S2_h * 0.75 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) ) +
      
      ninf[3] * sum(native * S3_h * 2 * 0.5 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S3_h * 2 * 0.5 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S3_h * 0.5 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S3_h * 0.5 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      native * S4_h * 2 * 0.25 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
                      travel * S4_h * 2 * 0.25 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      native * S4_h * 0.25 * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
                      travel * S4_h * 0.25 * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h)  )
  }
  
  ##s.l and h.l
  {  s.l.nv <- 
      inf[1] * sum( native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)) +
      inf[2] * sum(native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
      inf[3] * sum(native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)  +
                     native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                     travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                     native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                     travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l))
    
    h.l.nv <- 
      ninf[1] * sum(native * S1_l * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
                      travel * S1_l * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S1_l * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
                      travel * S1_l * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l)) +
      ninf[2] * sum(native * S2_l * beta_l * 2 * 0.75 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S2_l * beta_h * 2 *  0.75 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S2_l * beta_l * 0.75 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S2_l * beta_h * 0.75 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) ) +
      ninf[3] * sum(native * S3_l * beta_l * 2 * 0.5 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S3_l * beta_h * 2 *  0.5 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S3_l * beta_l * 0.5 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S3_l * beta_h * 0.5 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l)  +
                      native * S4_l * beta_l * 2 * 0.25 * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h)  / pop_h) +
                      travel * S4_l * beta_h * 2 *  0.25 * (native * (sym_inf_h)  / pop_h + travel * (sym_inf_l) / pop_l) +
                      native * S4_l * beta_l * 0.25 * (native * (inf_l) / pop_l + travel * (inf_h)  / pop_h) +
                      travel * S4_l * beta_h * 0.25 * (native * (inf_h)  / pop_h + travel * (inf_l) / pop_l) )}
  
  
  FOI_h.travel <-
    sum(native * 2 * beta_h * (native * (sym_inf_h) / pop_h + travel * (sym_inf_l) / pop_l) +
          travel * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
          native * beta_h * (native * (inf_h) / pop_h + travel * (inf_l) / pop_l) +
          travel * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h))
  
  FOI_l.travel <-
    sum(native * 2 * beta_l * (native * (sym_inf_l) / pop_l + travel * (sym_inf_h) / pop_h) +
          travel * 2 * beta_h * (native * (sym_inf_h)/ pop_h + travel * (sym_inf_l) / pop_l) +
          native * beta_l * (native * (inf_l) / pop_l + travel * (inf_h) / pop_h) +
          travel * beta_h * (native * (inf_h)/ pop_h + travel * (inf_l) / pop_l))
  zero_h <- sum(S1_h * vac_h * (1 - spec))
  one_h <- sum(sum((R1_h * vac_h * sens) + (S2_h * vac_h * sens)))
  two_h <- sum(sum((R2_h * vac_h * sens) + (S3_h * vac_h * sens)))
  three_h <- sum(sum((R3_h * vac_h * sens) + (S4_h * vac_h * sens)))
  
  
  zero_l <- sum(S1_l * vac_l * (1 - spec))
  one_l <- sum(sum((R1_l * vac_l * sens) + (S2_l * vac_l * sens)))
  two_l <- sum(sum((R2_l * vac_l * sens) + (S3_l * vac_l * sens)))
  three_l <- sum(sum((R3_l * vac_l * sens) + (S4_l * vac_l * sens)))
  
  
  tot_h <- sum(zero_h + one_h + two_h + three_h)
  prop_h <- c(zero_h / tot_h, one_h / tot_h, two_h / tot_h, three_h / tot_h)
  
  tot_l <- sum(zero_l + one_l + two_l + three_l)
  prop_l <- c(zero_l / tot_l, one_l / tot_l, two_l / tot_l, three_l / tot_l)
  
  
  
  
  list(c(
    dS1_h, dI1_h, dR1_h,
    dS2_h, dI2_h, dR2_h,
    dS3_h, dI3_h, dR3_h,
    dS4_h, dI4_h, dR4_h,
    dR1_h.v,
    dS2_h.v, dI2_h.v, dR2_h.v,
    dS3_h.v, dI3_h.v, dR3_h.v,
    dS4_h.v, dI4_h.v, dR4_h.v,
    
    dS1_l, dI1_l, dR1_l,
    dS2_l, dI2_l, dR2_l,
    dS3_l, dI3_l, dR3_l,
    dS4_l, dI4_l, dR4_l,
    dR1_l.v,
    dS2_l.v, dI2_l.v, dR2_l.v,
    dS3_l.v, dI3_l.v, dR3_l.v,
    dS4_l.v, dI4_l.v, dR4_l.v,
    
    I_secondary_vac, I_post_sec_vac,
    I_l_sec_vac, I_l_post_sec_vac,
    
    primary_tot.cases, primary_tot.ncases, 
    primary_cases.l, primary_ncases.l,
    primary_cases.h, primary_ncases.h,
    
    secondary_vac_tot.cases, secondary_vac_tot.ncases,
    secondary_vac_cases.l, secondary_vac_ncases.l,
    secondary_vac_cases.h, secondary_vac_ncases.h,
    secondary_tot.cases, secondary_tot.ncases,
    secondary_cases.l, secondary_ncases.l,
    secondary_cases.h, secondary_ncases.h,
    
    postsec_vac_tot.cases, postsec_vac_tot.ncases,
    postsec_vac_cases.l, postsec_vac_ncases.l,
    postsec_vac_cases.h, postsec_vac_ncases.h,
    postsec_tot_cases, postsec_tot_ncases,
    postsec_cases.l, postsec_ncases.l,
    postsec_cases.h, postsec_ncases.h,
    
    cases.h,
    
    I_total, I_l, cases, cases.l,
    I_primary_tot, I_secondary_tot, I_post_sec_tot,
    I_l_primary_tot, I_l_sec_tot, I_l_post_sec_tot,
    primary, secondary, post_sec,
    primary.l, secondary.l, post_sec.l, 
    primary_inf, secondary_inf, post_sec_inf, 
    primary.l.inf, secondary.l.inf, post_sec.l.inf,
    I_h,
    FOI_h.travel, FOI_l.travel,
    prop_h, prop_l
    # FOI_h, FOI_l,
    # zero_h, one_h, two_h, three_h,
    # zero_l, one_l, two_l, three_l
  ))
  
}

############################
#RUN MODEL
############################
y_init <- c(susceptible_h(1), infected_h(1), recovered_h(1),
            susceptible_h(2), infected_h(2), recovered_h(2),
            susceptible_h(3), infected_h(3), recovered_h(3),
            susceptible_h(4), infected_h(4), recovered_h(4),
            rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80),
            rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80),
            susceptible_l(1), infected_l(1), recovered_l(1),
            susceptible_l(2), infected_l(2), recovered_l(2),
            susceptible_l(3), infected_l(3), recovered_l(3),
            susceptible_l(4), infected_l(4), recovered_l(4),
            rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80),
            rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80), rep(0, 80),
            0,
            rep(0, 30),
            rep(0,29),
            rep(0,8))
names(y_init) <- c(rep('sh1', 80), rep('ih1', 80), rep('rh1', 80),
                   rep('sh2', 80), rep('ih2', 80), rep('rh2', 80),
                   rep('sh3', 80), rep('ih3', 80), rep('rh3', 80),
                   rep('sh4', 80), rep('ih4', 80), rep('rh4', 80),
                   rep('rh1.v', 80),
                   rep('sh2.v', 80), rep('ih2.v', 80), rep('rh2.v', 80),
                   rep('sh3.v', 80), rep('ih3.v', 80), rep('rh3.v', 80),
                   rep('sh4.v', 80), rep('ih4.v', 80), rep('rh4.v', 80),
                   rep('sl1', 80), rep('il1', 80), rep('rl1', 80),
                   rep('sl2', 80), rep('il2', 80), rep('rl2', 80),
                   rep('sl3', 80), rep('il3', 80), rep('rl3', 80),
                   rep('sl4', 80), rep('il4', 80), rep('rl4', 80),
                   rep('rl1.v', 80),
                   rep('sl2.v', 80), rep('il2.v', 80), rep('rl2.v', 80),
                   rep('sl3.v', 80), rep('il3.v', 80), rep('rl3.v', 80),
                   rep('sl4.v', 80), rep('il4.v', 80), rep('rl4.v', 80),
                   'i_sec_vac', 'i_psec_vac',
                   'il_sec_vac', 'il_psec_vac',
                   'prim_tot.cases', 'prim_tot.ncases',
                   'prim_cases.l', 'prim_ncases.l',
                   'prim_cases.h', 'prim_ncases.h',
                   
                   'sec_tot.cases.v', 'sec_tot.ncases.v',
                   'sec_vac.cases.l', 'sec_vac.ncases.l',
                   'sec_vac.cases.h', 'sec_vac.ncases.h',
                   'sec_tot.cases', 'sec_tot.ncases',
                   'sec_cases.l', 'sec_ncases.l',
                   'sec_cases.h', 'sec_ncases.h',
                   
                   'psec_vac_tot.cases', 'psec_vac_tot.ncases',
                   'psec_vac.cases.l', 'psec_vac.ncases.l',
                   'psec_vac.cases.h', 'psec_vac.ncases.h',
                   'psec_tot.cases', 'psec_tot.ncases',
                   'psec_cases.l', 'psec_ncases.l',
                   'psec_cases.h', 'psec_ncases.h',
                   
                   'cases.h',
                   
                   'i_total', 'il', 'cases', 'cases.l',
                   'i_prim_tot', 'i_sec_tot', 'i_post_sec_tot',
                   'il_prim_tot', 'il_sec_tot', 'il_psec_tot',
                   'prim', 'sec', 'psec',
                   'prim.l', 'sec.l', 'psec.l',
                   'prim_inf', 'sec_inf', 'psec_inf',
                   'prim.l.inf', 'sec.l.inf', 'psec.l.inf',
                   'ih',
                   'FOI_h.travel', 'FOI_l.travel',
                   # , 'FOI_h', 'FOI_l',
                   # 'zero_h', 'one_h', 'two_h', 'three_h',
                   # 'zero_l', 'one_l', 'two_l', 'three_l'
                   'prop_h', 'prop_l'
)
years = 60
years_vac = 30
out.h <- ode(times = times, y = y_init, func = model, parms = parms.h)
out_null.h <- ode(times = times, y = y_init, func = model, parms = parms_null.h)


#save(time, file = paste('time_', input, '.RData', sep = ''))

##incidence calcs 

out.h <- out.h[,2:ncol(out.h)]
out_null.h <- out_null.h[,2:ncol(out_null.h)]



{
#   ###Coverage calcs
#   vac_h <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh1.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh2.v')]) +
#                  sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih2.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh2.v')]) +
#                  sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh3.v')]) +
#                  sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh4.v')]))
#   
#   
#   vac_h.1 <- sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh1.v')])
#   vac_h.2 <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh2.v')]) +
#                     sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih2.v')]))
#   vac_h.3 <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh2.v')]) +
#                    sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh3.v')]) +
#                    sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh4.v')]))
#   
#   pop_h.1 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sh1')[1]):(which(colnames(out.h) == 'rh1'))[length(which(colnames(out.h) == 'rh1'))])])
#   pop_h.2 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sh2')[1]):(which(colnames(out.h) == 'rh2'))[length(which(colnames(out.h) == 'rh2'))])])
#   pop_h.3 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sh3')[1]):(which(colnames(out.h) == 'rh4'))[length(which(colnames(out.h) == 'rh4'))])])
#   
#   coverage_h <- c(vac_h.1 / pop_h.1, vac_h.2 / pop_h.2, vac_h.3 / pop_h.3)
#   
#   
#   vac_l.1 <- sum(out.h[nrow(out.h),which(colnames(out.h) == 'rl1.v')])
#   vac_l.2 <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'sl2.v')]) +
#                    sum(out.h[nrow(out.h),which(colnames(out.h) == 'il2.v')]))
#   vac_l.3 <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'rl2.v')]) +
#                    sum(out.h[nrow(out.h),which(colnames(out.h) == 'sl3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'il3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rl3.v')]) +
#                    sum(out.h[nrow(out.h),which(colnames(out.h) == 'sl4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'il4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rl4.v')]))
#   
#   pop_l.1 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sl1')[1]):(which(colnames(out.h) == 'rl1'))[length(which(colnames(out.h) == 'rl1'))])])
#   pop_l.2 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sl2')[1]):(which(colnames(out.h) == 'rl2'))[length(which(colnames(out.h) == 'rl2'))])])
#   pop_l.3 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sl3')[1]):(which(colnames(out.h) == 'rl4'))[length(which(colnames(out.h) == 'rl4'))])])
#   
#   coverage_l <- c(vac_l.1 / pop_l.1, vac_l.2 / pop_l.2, vac_l.3 / pop_l.3)
#   
#   coverage <- c(((vac_h.1 + vac_l.1) / (pop_h.1 + pop_l.1)),
#                 ((vac_h.2 + vac_l.2) / (pop_h.2 + pop_l.2)),
#                 ((vac_h.3 + vac_l.3) / (pop_h.3 + pop_l.3)))
#   
# 
# # 
# cases_averted.func <- function(out_mat, out_mat_null, timepoint_year){
#   indexing <- c((3650 * years_vac + 1):(timepoint_year * 3650 ))
#   # cases <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases')]))
#   # cases.l <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases.l')]))
#   # cases.h <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases.h')]))
#   # cases.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'cases')]))
#   # cases.l.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'cases.l')]))
#   # cases.h.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'cases.h')]))
# 
#   ##vaccinated, whole pop
#   cases_vac<- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.h')]),
#                   diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.l')])) *0.14
# 
# 
#   ##unvaccinated, whole pop
#   cases_uvac <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.h')]),
#                     diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.l')])) * 0.14
#   cases_uvac.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')]),
#                          diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])) * 0.14
# 
#   ##vaccinated, high
#   cases_vac.h <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.h')])) * 0.14
#   
# 
#   ##vaccinated, low
#   cases_vac.l <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.l')])) * 0.14
# 
#   
#   cases_vac.h.null <-sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')])) * coverage_h[2] * 0.14 
#   
#   cases_vac.null <- sum( c(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')]),
#                              diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])) * c(coverage_h[2], coverage_l[2]) * 0.14)
#   cases_vac.l.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])) * coverage_l[2] * 0.14
#   ##unvaccinaetd, high
#   cases_uvac.h <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.h')])) * 0.14
#   cases_uvac.h.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')])) * 0.14
# 
#   ##unvaccinated, low
#   cases_uvac.l <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.l')])) * 0.14
#   cases_uvac.l.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])) * 0.14
# 
#  cases.null <- cases_vac.null + cases_uvac.null
#  cases <- cases_vac + cases_uvac
#  cases.h.null <- cases_vac.h.null + cases_uvac.h.null
#  cases.h <- cases_vac.h + cases_uvac.h
#  cases.l.null <- cases_vac.l.null + cases_uvac.l.null
#  cases.l <- cases_vac.l + cases_uvac.l
#  
#   cases_averted <-  ((cases.null - cases) / cases.null) * 100
#   cases_averted.h <- ((cases.h.null - cases.h) / cases.h.null) * 100
#   cases_averted.l <- ((cases.l.null - cases.l) / cases.l.null) * 100
# 
#   cases_av_vac.h <- ((cases_vac.h.null - cases_vac.h) / cases_vac.h.null) * 100
#   cases_av_vac.l <- ((cases_vac.l.null - cases_vac.l) / cases_vac.l.null) * 100
#   cases_av_vac <- ((cases_vac.null - cases_vac) / cases_vac.null) * 100
# 
#   cases_av_uvac.h <- ((cases_uvac.h.null - cases_uvac.h) / cases_uvac.h.null) * 100
#   cases_av_uvac.l <- ((cases_uvac.l.null - cases_uvac.l) / cases_uvac.l.null) * 100
#   cases_av_uvac <- ((cases_uvac.null - cases_uvac) / cases_uvac.null) * 100
# 
# 
#   output <- c(cases_averted.h, cases_averted, cases_averted.l,
#               cases_av_vac.h, cases_av_vac, cases_av_vac.l,
#               cases_av_uvac.h, cases_av_uvac, cases_av_uvac.l)
#   names(output) <- c('cases_averted.h', 'cases_averted', 'cases_averted.l',
#                      'cases_av_vac.h', 'cases_av_vac', 'cases_av_vac.l',
#                      'cases_av_uvac.h', 'cases_av_uvac', 'cases_av_uvac.l')
# 
#   return(output)
# }
# #
# output.vec.h  <- cases_averted.func(out.h, out_null.h, years)
# 
# save(output.vec.h, file = paste('output_alttrav_test_', input, '.RData', sep = ''))
}


###################
#old calcs
###################
{
  ###Coverage calcs
  vac_h <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh1.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh2.v')]) +
                 sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih2.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh2.v')]) +
                 sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh3.v')]) +
                 sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh4.v')]))
  
  
  vac_h.1 <- sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh1.v')])
  vac_h.2 <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh2.v')]) +
                   sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih2.v')]))
  vac_h.3 <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh2.v')]) +
                   sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh3.v')]) +
                   sum(out.h[nrow(out.h),which(colnames(out.h) == 'sh4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'ih4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rh4.v')]))
  
  pop_h.1 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sh1')[1]):(which(colnames(out.h) == 'rh1'))[length(which(colnames(out.h) == 'rh1'))])])
  pop_h.2 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sh2')[1]):(which(colnames(out.h) == 'rh2'))[length(which(colnames(out.h) == 'rh2'))])])
  pop_h.3 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sh3')[1]):(which(colnames(out.h) == 'rh4'))[length(which(colnames(out.h) == 'rh4'))])])
  
  coverage_h <- c(vac_h.1 / pop_h.1, vac_h.2 / pop_h.2, vac_h.3 / pop_h.3)
  
  
  vac_l.1 <- sum(out.h[nrow(out.h),which(colnames(out.h) == 'rl1.v')])
  vac_l.2 <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'sl2.v')]) +
                   sum(out.h[nrow(out.h),which(colnames(out.h) == 'il2.v')]))
  vac_l.3 <- sum(sum(out.h[nrow(out.h),which(colnames(out.h) == 'rl2.v')]) +
                   sum(out.h[nrow(out.h),which(colnames(out.h) == 'sl3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'il3.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rl3.v')]) +
                   sum(out.h[nrow(out.h),which(colnames(out.h) == 'sl4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'il4.v')]) + sum(out.h[nrow(out.h),which(colnames(out.h) == 'rl4.v')]))
  
  pop_l.1 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sl1')[1]):(which(colnames(out.h) == 'rl1'))[length(which(colnames(out.h) == 'rl1'))])])
  pop_l.2 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sl2')[1]):(which(colnames(out.h) == 'rl2'))[length(which(colnames(out.h) == 'rl2'))])])
  pop_l.3 <- sum(out.h[nrow(out.h),((which(colnames(out.h) == 'sl3')[1]):(which(colnames(out.h) == 'rl4'))[length(which(colnames(out.h) == 'rl4'))])])
  
  coverage_l <- c(vac_l.1 / pop_l.1, vac_l.2 / pop_l.2, vac_l.3 / pop_l.3)
  
  coverage <- c(((vac_h.1 + vac_l.1) / (pop_h.1 + pop_l.1)),
                ((vac_h.2 + vac_l.2) / (pop_h.2 + pop_l.2)),
                ((vac_h.3 + vac_l.3) / (pop_h.3 + pop_l.3)))
  
  
  # 
  cases_averted.func <- function(out_mat, out_mat_null, timepoint_year){
    indexing <- c((3650 * years_vac + 1):(timepoint_year * 3650 ))
    # cases <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases')]))
    # cases.l <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases.l')]))
    # cases.h <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'cases.h')]))
    # cases.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'cases')]))
    # cases.l.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'cases.l')]))
    # cases.h.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'cases.h')]))
    
    ##vaccinated, whole pop
    cases_vac<- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.h')]),
                    diff(out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.h')]),
                    diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.l')]),
                    diff(out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.l')]))
    
    
    ##unvaccinated, whole pop
    cases_uvac <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'prim_cases.h')]),
                      diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.h')]),
                      diff(out_mat[indexing,which(colnames(out_mat) == 'psec_cases.h')]),
                      diff(out_mat[indexing,which(colnames(out_mat) == 'prim_cases.l')]),
                      diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.l')]),
                      diff(out_mat[indexing,which(colnames(out_mat) == 'psec_cases.l')]))
    cases_uvac.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')]),
                           diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')]),
                           diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')]),
                           diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')]),
                           diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')]),
                           diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))
    
    ##vaccinated, high
    cases_vac.h <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.h')]),
                       diff(out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.h')]))
    
    
    ##vaccinated, low
    cases_vac.l <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.l')]),
                       diff(out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.l')]))
    
    
    cases_vac.h.null <-sum( c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')])),
                              sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')])),
                              sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')]))) * coverage_h)
    cases_vac.null <- sum( c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')])),
                             sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')])),
                             sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')])),
                             sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')])),
                             sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])),
                             sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))) * c(coverage_h, coverage_l))
    cases_vac.l.null <- sum( c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')])),
                               sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])),
                               sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))) * coverage_l)
    ##unvaccinaetd, high
    cases_uvac.h <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'prim_cases.h')]),
                        diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.h')]),
                        diff(out_mat[indexing,which(colnames(out_mat) == 'psec_cases.h')]))
    cases_uvac.h.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')]),
                             diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')]),
                             diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')]))
    
    ##unvaccinated, low
    cases_uvac.l <- sum(diff(out_mat[indexing,which(colnames(out_mat) == 'prim_cases.l')]),
                        diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.l')]),
                        diff(out_mat[indexing,which(colnames(out_mat) == 'psec_cases.l')]))
    cases_uvac.l.null <- sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')]),
                             diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')]),
                             diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))
    
    cases.null <- cases_vac.null + cases_uvac.null
    cases <- cases_vac + cases_uvac
    cases.h.null <- cases_vac.h.null + cases_uvac.h.null
    cases.h <- cases_vac.h + cases_uvac.h
    cases.l.null <- cases_vac.l.null + cases_uvac.l.null
    cases.l <- cases_vac.l + cases_uvac.l
    
    cases_averted <-  ((cases.null - cases) / cases.null) * 100
    cases_averted.h <- ((cases.h.null - cases.h) / cases.h.null) * 100
    cases_averted.l <- ((cases.l.null - cases.l) / cases.l.null) * 100
    
    cases_av_vac.h <- ((cases_vac.h.null - cases_vac.h) / cases_vac.h.null) * 100
    cases_av_vac.l <- ((cases_vac.l.null - cases_vac.l) / cases_vac.l.null) * 100
    cases_av_vac <- ((cases_vac.null - cases_vac) / cases_vac.null) * 100
    
    cases_av_uvac.h <- ((cases_uvac.h.null - cases_uvac.h) / cases_uvac.h.null) * 100
    cases_av_uvac.l <- ((cases_uvac.l.null - cases_uvac.l) / cases_uvac.l.null) * 100
    cases_av_uvac <- ((cases_uvac.null - cases_uvac) / cases_uvac.null) * 100
    
    
    output <- c(cases_averted.h, cases_averted, cases_averted.l,
                cases_av_vac.h, cases_av_vac, cases_av_vac.l,
                cases_av_uvac.h, cases_av_uvac, cases_av_uvac.l)
    names(output) <- c('cases_averted.h', 'cases_averted', 'cases_averted.l',
                       'cases_av_vac.h', 'cases_av_vac', 'cases_av_vac.l',
                       'cases_av_uvac.h', 'cases_av_uvac', 'cases_av_uvac.l')
    
    return(output)
  }
  #
  output.vec.h  <- cases_averted.func(out.h, out_null.h, years)
  
  save(output.vec.h, file = paste('output_alttrav_test_', input, '.RData', sep = ''))
}


######check seroprevalence levels
out <- out.h
  nines_h <- c(which(colnames(out) == 'sh1')[10], which(colnames(out) == 'ih1')[10], which(colnames(out) == 'rh1')[10],
               which(colnames(out) == 'sh2')[10], which(colnames(out) == 'ih2')[10], which(colnames(out) == 'rh2')[10],
               which(colnames(out) == 'sh3')[10], which(colnames(out) == 'ih3')[10], which(colnames(out) == 'rh3')[10],
               which(colnames(out) == 'sh4')[10], which(colnames(out) == 'ih4')[10], which(colnames(out) == 'rh4')[10],
               which(colnames(out) == 'rh1.v')[10],
               which(colnames(out) == 'sh2.v')[10], which(colnames(out) == 'ih2.v')[10], which(colnames(out) == 'rh2.v')[10],
               which(colnames(out) == 'sh3.v')[10], which(colnames(out) == 'ih3.v')[10], which(colnames(out) == 'rh3.v')[10],
               which(colnames(out) == 'sh4.v')[10], which(colnames(out) == 'ih4.v')[10], which(colnames(out) == 'rh4.v')[10])
  nines_l <- c(which(colnames(out) == 'sl1')[10], which(colnames(out) == 'il1')[10], which(colnames(out) == 'rl1')[10],
               which(colnames(out) == 'sl2')[10], which(colnames(out) == 'il2')[10], which(colnames(out) == 'rl2')[10],
               which(colnames(out) == 'sl3')[10], which(colnames(out) == 'il3')[10], which(colnames(out) == 'rl3')[10],
               which(colnames(out) == 'sl4')[10], which(colnames(out) == 'il4')[10], which(colnames(out) == 'rl4')[10],
               which(colnames(out) == 'rl1.v')[10],
               which(colnames(out) == 'sl2.v')[10], which(colnames(out) == 'il2.v')[10], which(colnames(out) == 'rl2.v')[10],
               which(colnames(out) == 'sl3.v')[10], which(colnames(out) == 'il3.v')[10], which(colnames(out) == 'rl3.v')[10],
               which(colnames(out) == 'sl4.v')[10], which(colnames(out) == 'il4.v')[10], which(colnames(out) == 'rl4.v')[10])
  nines <- c(nines_h, nines_l)


i <- 365 * 10 * years_vac
no_exposure <- out.h[i, nines[1]] + out.h[i, nines[13]] + out.h[i, nines[14]] + out.h[i, nines[23]] + out.h[i, nines[35]] + out.h[i, nines[36]]
no_exposure_h <- out.h[i, nines_h[1]] + out.h[i, nines_h[13]] + out.h[i, nines_h[14]]
no_exposure_l <- out.h[i, nines_l[1]] + out.h[i, nines_l[13]] + out.h[i, nines_l[14]]

sp9 <- 1 - (no_exposure / sum(out.h[i, nines]))
sp9_h <- 1 - (no_exposure_h / sum(out.h[i, nines_h]))
sp9_l <- 1 - (no_exposure_l / sum(out.h[i, nines_l]))

sp9.vec <- c(sp9, sp9_h, sp9_l)
save(sp9.vec, file = paste('sp9_', input, '.RData', sep = ''))










