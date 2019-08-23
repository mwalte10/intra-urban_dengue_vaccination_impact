###########################
##INITIAL CONDITIONS
###########################
{
load('birth_1950.RData')
load('death_1950.RData')
load('recovered_init_h.RData')
load('recovered_init_l.RData')


##set up population distribution
load('pop_1950.RData')
initial_conditions <- as.data.frame(matrix(NA, nrow = 80*4, ncol = 3))
initial_conditions[,2] <- rep(0:79,4)
for(i in 1:4){
  x <- i - 1
  initial_conditions[x*(80) + (1:80),1] <- rep(i,80)
}
initial_conditions[,3] <- rep(pop,4)

##distributes susceptible individuals for the four serotypes
susceptible_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
susceptible_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
infected_init_h <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))
infected_init_l <- cbind(initial_conditions[,1], initial_conditions[,2], rep(NA, 320))

fill_suscep <- function(prop.h.1, prop.h.2, prop.h.3, prop.h.4,
                        prop.l.1, prop.l.2, prop.l.3, prop.l.4){
  if(sum(prop.h.1 + prop.h.2 + prop.h.3 + prop.h.4 +
         prop.l.1 + prop.l.2 + prop.l.3 + prop.l.4) != 1)stop('The sum of proportions must equal 1')
  h <- c(initial_conditions[1:80,3] * prop.h.1, initial_conditions[81:160,3] * prop.h.2, 
                              initial_conditions[161:240,3] * prop.h.3, initial_conditions[241:320,3] * prop.h.4)
  l <- c(initial_conditions[1:80,3] * prop.l.1, initial_conditions[81:160,3] * prop.l.2, 
                              initial_conditions[161:240,3] * prop.l.3, initial_conditions[241:320,3] * prop.l.4)
  return(list(h,l))
  
}

fill_inf <- function(num.inf.h, num.inf.l){
  h <- c(rep(num.inf.h/80,80), initial_conditions[81:160,3] * 0, 
         initial_conditions[161:240,3] * 0, initial_conditions[241:320,3] * 0)
  l <- c(rep(num.inf.l/80,80), initial_conditions[81:160,3] * 0, 
         initial_conditions[161:240,3] * 0, initial_conditions[241:320,3] * 0)
  return(list(h,l))
  
}






####functions to obtain the correct initial conditions
susceptible_h <- function(exposure){
  x <- exposure - 1
  susceptible <- c(susceptible_init_h[(x * 80) + 1:80,3])
  return(susceptible)
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
}

###########################
##CASES AVERTED
###########################

  cases_averted.func <- function(out_mat, out_mat_null, timepoint_year){
  indexing <- (timepoint_year * 3650)
  # if(cases == 1){ ##vaccinated, whole pop
    cases_vac<- sum((out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.h')]),
                    (out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.h')]),
                    (out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.l')]),
                    (out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.l')]))
    
    cases_uvac_eleg.h.null <- sum((out_mat_null[indexing,which(colnames(out_mat_null) == 'vac_eleg_p.h')]),
                           (out_mat_null[indexing,which(colnames(out_mat_null) == 'vac_eleg_s.h')]),
                           (out_mat_null[indexing,which(colnames(out_mat_null) == 'vac_eleg.ps.h')]))
    cases_uvac_eleg.l.null <-  sum((out_mat_null[indexing,which(colnames(out_mat_null) == 'vac_eleg_p.l')]),
                              (out_mat_null[indexing,which(colnames(out_mat_null) == 'vac_eleg_s.l')]),
                              (out_mat_null[indexing,which(colnames(out_mat_null) == 'vac_eleg.ps.l')]))
    
    cases_uvac_eleg.null <- cases_uvac_eleg.l.null + cases_uvac_eleg.h.null
    
    cases_uvac_eleg.h <- sum((out_mat[indexing,which(colnames(out_mat_null) == 'vac_eleg_p.h')]),
                                  (out_mat[indexing,which(colnames(out_mat_null) == 'vac_eleg_s.h')]),
                                  (out_mat[indexing,which(colnames(out_mat_null) == 'vac_eleg.ps.h')]))
    cases_uvac_eleg.l <-  sum((out_mat[indexing,which(colnames(out_mat_null) == 'vac_eleg_p.l')]),
                                   (out_mat[indexing,which(colnames(out_mat_null) == 'vac_eleg_s.l')]),
                                   (out_mat[indexing,which(colnames(out_mat_null) == 'vac_eleg.ps.l')]))
    
    cases_uvac_eleg <- cases_uvac_eleg.l + cases_uvac_eleg.h
    
    
    ##unvaccinated, whole pop
    cases_uvac <- sum((out_mat[indexing,which(colnames(out_mat) == 'prim_cases.h')]),
                      (out_mat[indexing,which(colnames(out_mat) == 'sec_cases.h')]),
                      (out_mat[indexing,which(colnames(out_mat) == 'psec_cases.h')]),
                      (out_mat[indexing,which(colnames(out_mat) == 'prim_cases.l')]),
                      (out_mat[indexing,which(colnames(out_mat) == 'sec_cases.l')]),
                      (out_mat[indexing,which(colnames(out_mat) == 'psec_cases.l')]))
    cases_uvac.null <- sum((out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')]),
                           (out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')]),
                           (out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')]),
                           (out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')]),
                           (out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')]),
                           (out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))
    
    ##vaccinated, high
    cases_vac.h <- sum((out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.h')]) + 
                       (out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.h')]))
    
    
    
    ##vaccinated, low
    cases_vac.l <- sum((out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.l')]) + 
                       (out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.l')]))
    

    ##unvaccinaetd, high, these should be the same number by they aren't
    cases_uvac.h <- sum((out_mat[indexing,which(colnames(out_mat) == 'prim_cases.h')]),
                        (out_mat[indexing,which(colnames(out_mat) == 'sec_cases.h')]),
                        (out_mat[indexing,which(colnames(out_mat) == 'psec_cases.h')]))
    
    ##psec_cases.h null > not null
    cases_uvac.h.null <- sum((out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')]),
                             (out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')]),
                             (out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')]))
    

    
    ##unvaccinated, low
    cases_uvac.l <- sum((out_mat[indexing,which(colnames(out_mat) == 'prim_cases.l')]),
                        (out_mat[indexing,which(colnames(out_mat) == 'sec_cases.l')]),
                        (out_mat[indexing,which(colnames(out_mat) == 'psec_cases.l')]))
    cases_uvac.l.null <- sum((out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')]),
                             (out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')]),
                             (out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))
    
    cases.h.null <- sum((out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')]),
                        (out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')]),
                        (out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')]),
                        (out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_vac.cases.h')]) + 
                          (out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_vac.cases.h')]))
    
    cases.l.null <- sum((out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')]),
                        (out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')]),
                        (out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]),
                        (out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_vac.cases.l')]) + 
                          (out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_vac.cases.l')]))
    
    cases.null <- cases_uvac.null
    cases <- cases_vac + cases_uvac
    # cases.h.null <- cases_vac.h.null + cases_uvac.h.null
    cases.h <- cases_vac.h + cases_uvac.h
    # cases.l.null <- cases_vac.l.null + cases_uvac.l.null
    cases.l <- cases_vac.l + cases_uvac.l
    
    prop.h <- cases_vac.h / (cases_uvac_eleg.h)
    cases_vac.h.null <- prop.h * cases_uvac_eleg.h.null
    prop.l <- cases_vac.l / (cases_uvac_eleg.l)
    cases_vac.l.null <- prop.l * cases_uvac_eleg.l.null
    
    prop.tot <- (cases_vac.h + cases_vac.l) / (cases_uvac_eleg)
    cases_vac.null <- prop.tot * cases_uvac_eleg.null
    
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
    
    return(output)}
#   if( cases == 0){
#     inf <- hopkins
#     
#     
#     inf_vac<- sum(c(sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.h')])),
#                     sum(diff(out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.h')])),
#                     sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.l')])),
#                     sum(diff(out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.l')]))) / c(inf[2], inf[3], inf[2], inf[3]))
#     
#     
#     ##unvaccinated, whole pop
#     inf_uvac <- sum(c(sum(diff(out_mat[indexing,which(colnames(out_mat) == 'prim_cases.h')])),
#                       sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.h')])),
#                       sum(diff(out_mat[indexing,which(colnames(out_mat) == 'psec_cases.h')])),
#                       sum(diff(out_mat[indexing,which(colnames(out_mat) == 'prim_cases.l')])),
#                       sum( diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.l')])),
#                       sum(diff(out_mat[indexing,which(colnames(out_mat) == 'psec_cases.l')]))) /
#                       c(inf[1], inf[2], inf[3],
#                         inf[1], inf[2], inf[3]))
#     inf_uvac.null <- sum(c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')])),
#                            sum( diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))) /
#                            c(inf[1], inf[2], inf[3],
#                              inf[1], inf[2], inf[3]))
#     
#     ##vaccinated, high
#     inf_vac.h <- sum(c(sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.h')])),
#                        sum(diff(out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.h')]))) / 
#                        c(inf[2], inf[3]))
#     
#     
#     ##vaccinated, low
#     inf_vac.l <- sum(c(sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_vac.cases.l')])),
#                        sum(diff(out_mat[indexing,which(colnames(out_mat) == 'psec_vac.cases.l')]))) / 
#                        c(inf[2], inf[3]))
#     
#     
#     inf_vac.h.null <-sum( c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')])),
#                             sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')])),
#                             sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')]))) / c(inf[1], inf[2], inf[3]) * coverage_h)
#     
#     inf_vac.null <- sum( c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])),
#                            sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))) / c(inf[1], inf[2], inf[3], inf[1], inf[2], inf[3]) * c(coverage_h, coverage_l))
#     inf_vac.l.null <- sum( c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')])),
#                              sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])),
#                              sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))) / c(inf[1], inf[2], inf[3])* coverage_l)
#     ##unvaccinaetd, high
#     inf_uvac.h <- sum(c(sum(diff(out_mat[indexing,which(colnames(out_mat) == 'prim_cases.h')])),
#                         sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.h')])),
#                         sum(diff(out_mat[indexing,which(colnames(out_mat) == 'psec_cases.h')]))) * c(inf[1], inf[2], inf[3]))
#     inf_uvac.h.null <- sum(c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.h')])),
#                              sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.h')])),
#                              sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.h')]))) * c(inf[1], inf[2], inf[3]))
#     
#     ##unvaccinated, low
#     inf_uvac.l <-sum(c(sum(diff(out_mat[indexing,which(colnames(out_mat) == 'prim_cases.l')])),
#                        sum(diff(out_mat[indexing,which(colnames(out_mat) == 'sec_cases.l')])),
#                        sum(diff(out_mat[indexing,which(colnames(out_mat) == 'psec_cases.l')]))) * c(inf[1], inf[2], inf[3]))
#     inf_uvac.l.null <- sum(c(sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'prim_cases.l')])),
#                              sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'sec_cases.l')])),
#                              sum(diff(out_mat_null[indexing,which(colnames(out_mat_null) == 'psec_cases.l')]))) * c(inf[1], inf[2], inf[3]))
#     
#     inf.null <- inf_vac.null + inf_uvac.null
#     inf <- inf_vac + inf_uvac
#     inf.h.null <- inf_vac.h.null + inf_uvac.h.null
#     inf.h <- inf_vac.h + inf_uvac.h
#     inf.l.null <- inf_vac.l.null + inf_uvac.l.null
#     inf.l <- inf_vac.l + inf_uvac.l
#     
#     inf_averted <-  ((inf.null - inf) / inf.null) * 100
#     inf_averted.h <- ((inf.h.null - inf.h) / inf.h.null) * 100
#     inf_averted.l <- ((inf.l.null - inf.l) / inf.l.null) * 100
#     
#     inf_av_vac.h <- ((inf_vac.h.null - inf_vac.h) / inf_vac.h.null) * 100
#     inf_av_vac.l <- ((inf_vac.l.null - inf_vac.l) / inf_vac.l.null) * 100
#     inf_av_vac <- ((inf_vac.null - inf_vac) / inf_vac.null) * 100
#     
#     inf_av_uvac.h <- ((inf_uvac.h.null - inf_uvac.h) / inf_uvac.h.null) * 100
#     inf_av_uvac.l <- ((inf_uvac.l.null - inf_uvac.l) / inf_uvac.l.null) * 100
#     inf_av_uvac <- ((inf_uvac.null - inf_uvac) / inf_uvac.null) * 100
#     
#     
#     output <- c(inf_averted.h, inf_averted, inf_averted.l,
#                 inf_av_vac.h, inf_av_vac, inf_av_vac.l,
#                 inf_av_uvac.h, inf_av_uvac, inf_av_uvac.l)
#     names(output) <- c('inf_averted.h', 'inf_averted', 'inf_averted.l',
#                        'inf_av_vac.h', 'inf_av_vac', 'inf_av_vac.l',
#                        'inf_av_uvac.h', 'inf_av_uvac', 'inf_av_uvac.l')
#     
#     
#     
#     return(output)}
#  
# }
#




###########################
##SEROPREVALENCE
###########################
{
  seroprevalence_fun <- function(out_mat, age){
  out <- out_mat
  index <- age + 1
  nines_h <- c(which(colnames(out) == 'sh1')[index], which(colnames(out) == 'ih1')[index], which(colnames(out) == 'rh1')[index],
               which(colnames(out) == 'sh2')[index], which(colnames(out) == 'ih2')[index], which(colnames(out) == 'rh2')[index],
               which(colnames(out) == 'sh3')[index], which(colnames(out) == 'ih3')[index], which(colnames(out) == 'rh3')[index],
               which(colnames(out) == 'sh4')[index], which(colnames(out) == 'ih4')[index], which(colnames(out) == 'rh4')[index],
               which(colnames(out) == 'rh1.v')[index],
               which(colnames(out) == 'sh2.v')[index], which(colnames(out) == 'ih2.v')[index], which(colnames(out) == 'rh2.v')[index],
               which(colnames(out) == 'sh3.v')[index], which(colnames(out) == 'ih3.v')[index], which(colnames(out) == 'rh3.v')[index],
               which(colnames(out) == 'sh4.v')[index], which(colnames(out) == 'ih4.v')[index], which(colnames(out) == 'rh4.v')[index])
  nines_l <- c(which(colnames(out) == 'sl1')[index], which(colnames(out) == 'il1')[index], which(colnames(out) == 'rl1')[index],
               which(colnames(out) == 'sl2')[index], which(colnames(out) == 'il2')[index], which(colnames(out) == 'rl2')[index],
               which(colnames(out) == 'sl3')[index], which(colnames(out) == 'il3')[index], which(colnames(out) == 'rl3')[index],
               which(colnames(out) == 'sl4')[index], which(colnames(out) == 'il4')[index], which(colnames(out) == 'rl4')[index],
               which(colnames(out) == 'rl1.v')[index],
               which(colnames(out) == 'sl2.v')[index], which(colnames(out) == 'il2.v')[index], which(colnames(out) == 'rl2.v')[index],
               which(colnames(out) == 'sl3.v')[index], which(colnames(out) == 'il3.v')[index], which(colnames(out) == 'rl3.v')[index],
               which(colnames(out) == 'sl4.v')[index], which(colnames(out) == 'il4.v')[index], which(colnames(out) == 'rl4.v')[index])
  nines <- c(nines_h, nines_l)
  
  
  i <- 365 * 10 * years_vac + 1
  no_exposure <- out.h[i, nines[1]] + out.h[i, nines[13]] + out.h[i, nines[14]] + out.h[i, nines[23]] + out.h[i, nines[35]] + out.h[i, nines[36]]
  no_exposure_h <- out.h[i, nines_h[1]] + out.h[i, nines_h[13]] + out.h[i, nines_h[14]]
  no_exposure_l <- out.h[i, nines_l[1]] + out.h[i, nines_l[13]] + out.h[i, nines_l[14]]
  
  sp9 <- 1 - (no_exposure / sum(out.h[i, nines]))
  sp9_h <- 1 - (no_exposure_h / sum(out.h[i, nines_h]))
  sp9_l <- 1 - (no_exposure_l / sum(out.h[i, nines_l]))
  
  sp9.vec <- c(sp9, sp9_h, sp9_l)
  names(sp9.vec) <- c('sp9', 'h', 'l')
  return(sp9.vec)
}


}



###########################
##FOI
###########################
FOI_h.fun <- function(years){
  i <- years* 3650
  vac <- out.h[i,which(colnames(out.h) == 'FOI_h.travel')]
  null <- out_null.h[i,which(colnames(out_null.h) == 'FOI_h.travel')]
  
  return(c(vac,null))
}

FOI_l.fun <- function(years){
  i <- years* 3650
  vac <- out.h[i,which(colnames(out.h) == 'FOI_l.travel')]
  null <- out_null.h[i,which(colnames(out_null.h) == 'FOI_l.travel')]
  
  return(c(vac,null))
}




