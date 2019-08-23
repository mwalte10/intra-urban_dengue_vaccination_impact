#####################
#FIGURE 2
#####################
{
  library(viridis)
  library(RColorBrewer)
  low.full <- plasma(2)[1]
  high.full <- 'darkorange'
  low <- adjustcolor(low.full, 0.7)
  high <- adjustcolor(high.full, 0.7)
  black <- adjustcolor('black', 0.7)
setwd('~/Desktop/MANUSCRIPT/FIGURE_2')
load('birth_1950.RData')
load('death_1950.RData')

#FIGURE 2A
{par(pty = 's', mar = c(1,7,1,1), xpd = FALSE, mfrow = c(1,2), oma = rep(0.2,4))
plot(birth*1000*365, type = 'l', lwd = 2,
     main = '', xaxt= 'n', yaxt = 'n', ylim = c(0,50),
     ylab = '', xlab = '')
axis(1, at = seq(0,90, by = 30), seq(1950, 2040, by = 30), cex.axis = 1.5)
abline(v=60, lwd = 2)
axis(2, at = seq(0,50, by =10), seq(0,50, by =10), cex.axis = 2, las  = 1.5)
mtext(side= 1, line = 4.5, expression('Year'), cex = 2)
mtext(side= 2, line = 4.5, expression('Births per 1,000'), cex = 2)
mtext(at = 0, line = 0.25, 'A', cex = 2)}

#FIGURE 2B
{fun <- function(i,x){i[x]}
virdis.col <- viridis(4)
newborns <- unlist(lapply(death, fun, x =1)) * 1000 * 365
ten <- unlist(lapply(death, fun, x =10)) * 1000 * 365
fourty <- unlist(lapply(death, fun, x = 40)) * 1000 * 365
eighty.plus <- unlist(lapply(death, fun, x = 70)) * 1000 * 365

plot(newborns, ylim = c(0,110), type = 'l',lwd = 2, col = virdis.col[1],
     xaxt= 'n', yaxt = 'n', 
     ylab = '', xlab = '')
lines(ten, lwd = 2, col = virdis.col[2])
lines(fourty, lwd = 2,  col = virdis.col[3])
lines(eighty.plus, lwd = 2, col = virdis.col[4])
axis(1, at = seq(0,90, by = 30), seq(1950, 2040, by = 30), cex.axis = 1.5)
axis(2, at = seq(0,110, by =20), seq(0,110, by =20), cex.axis = 1.5, las  =2)
legend('topright', bty = 'n', legend = c(0, 10, 40, 70), title = 'Age',
       col = virdis.col, lwd = rep(2,4))
mtext(side= 1, line = 4.5, expression('Year'), cex = 2)
mtext(side= 2, line = 4.5, expression('Deaths per 1,000'), cex = 2)
mtext(at = 0, line = 0.25, 'B', cex = 2)
}

}


#####################
#FIGURE 3
#####################
{
setwd('~/Desktop/MANUSCRIPT/FIGURE_3/')
  ##loads in results_h, the prop of cases averted from my model
load('model_val_results.RData')
##loads in case_av_flasche, summary statistics of the prop of cases averted from the models in the flasche paper
load('cases_averted_flasche_models.RData')


#########
#GROUPED BARCHART
#########
par(pty = 's', mar = c(10,10,4,13), xpd = FALSE, mfrow = c(1,1))

means_cohort <- case_av_flasche[[1]]
max_cohort <- case_av_flasche[[2]]
min_cohort <- case_av_flasche[[3]]

par(pty = 's', mar = c(7,7,1,1), xpd = FALSE, mfrow = c(1,1), oma = rep(0.2,4))
i_grouped.bc <- rbind(means_cohort, results_h)
par(xpd = FALSE)
colnames(i_grouped.bc) <- c('10%', '30%', '50%', '70%', '90%')
barCenters <- barplot(i_grouped.bc, beside = TRUE, ylim = c(-100, 27), yaxt = 'n', col = c('grey', 'black'),
        xaxt= 'n')
segments(barCenters, unlist(min_cohort), barCenters, unlist(max_cohort), lwd = 1.5)
arrows(barCenters, unlist(min_cohort), barCenters, unlist(max_cohort), lwd = 1.5,
       angle = 90, code = 3, length = 0.05)
axis(2, labels = paste(seq(-100,25, by = 25),'%',sep = ''), at = seq(-100,25, by = 25), las = 2, cex.axis = 1.5)
axis(1, labels = paste(c(10,30,50,70,90),'%',sep =''), at = c(2,5,8,11,14), las = 1, cex.axis = 1.5)
box()
abline(h = 0)
mtext(side = 1, line = 4.5, expression('Seroprevalence'), cex= 2)
mtext(side = 2, line = 6, expression('Proportion of Cases Averted'), cex = 2)
legend('bottomright', legend = c('Flasche et al.,\nModels', 'Our Model'), 
       fill = c('grey', 'black'), bty = 'n', cex = 1)


}


#####################
#FIGURE 4
#####################
{
  par(pty = 's', mar = c(7,5,1,7.5), xpd = FALSE, mfrow = c(1,1), oma = rep(0.4,4))
  setwd('~/Desktop/MANUSCRIPT/FIGURE_4/')
###loads in the coverage rates for no travel
load('coverage_notravel.RData')
coverage_h <- coverage_list.notravel[[1]]
coverage_l <- coverage_list.notravel[[2]]

vac_h.vec <- seq(0, 1, length.out = 20)
vac_l.vec <- rev(vac_h.vec)



plot(x = vac_h.vec, y = coverage_h, ylim = c(0, 70), type = 'l', yaxt = 'n', xaxt = 'n', col = low,
     ylab = '', xlab = '',
     lwd = 2)
lines(x = vac_h.vec, y = coverage_l, col = high, lwd = 2)
mtext(side = 1, line = 6, expression('Intervention Coverage'), cex = 2)
mtext(side = 2, line = 4.5, expression('Vaccination Coverage'), cex = 2)
axis(2, paste(seq(0, 70, by = 10), '%', sep = ''), at = seq(0, 70, by = 10), las = 2, cex.axis = 1.5)
axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(0,.5, 1), cex.axis = 1.5)
axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(0,.5, 1), cex.axis = 1.5)
legend(x = 1.03, y = 65,  legend = c('High-transmission\nNo mobility\n', 
                                  'Low-transmission\nNo mobility\n'),
       col = (c(high.full,low.full)), lty = c(1,1),bty = 'n', xpd = TRUE, 
       lwd =  rep(2, 4), cex= 1)
}


#####################
#FIGURE 5
#####################
{
  
setwd('~/Desktop/MANUSCRIPT/FIGURE_5/')
  ##LOAD IN THE CHANGE IN FOI UNDER ALL CONDITIONS
load('change_in_foi.RData')
h.notravel <- change_foi[[1]]
l.notravel <- change_foi[[2]]
h.travel <- change_foi[[3]]
l.travel <- change_foi[[4]]


#FIGURE 5A
{par(pty = 's', mar = c(5, 7, 0.5,0),  xpd = FALSE, mfrow = c(1,2), oma = c(0,0,0,2))
plot(h.notravel, col = low, type = 'l', ylim = c(-15, -0), xaxt = 'n', yaxt = 'n',
     ylab = '', xlab = '', lwd = 2)
mtext(at = 0, line = 0.25, 'A', cex = 2)

axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(1,10.5, 20), cex.axis = 1.5)
axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(1,10.5, 20), cex.axis = 1.5)
axis(2, paste((seq(0, -15, by = -5)), '%', sep = ''), at = seq(0, -15, by = -5), las = 2, cex.axis = 1.5)
mtext(side = 2, line = 4.5, expression('Force of Infection'), cex = 2)
mtext(side = 3, line = 2, expression('No Mobility'), cex = 2)
mtext(side = 1, line = 6, expression('Intervention Coverage'), cex = 2)
lines(l.notravel, col = high, lwd = 2)}

#FIGURE 5B
{par(mar = c(5, 6, 0.5,1))
plot(h.travel, col = low, 
     type = 'l', ylim = c(-15, -0), xaxt = 'n', yaxt = 'n',
     ylab = '', xlab = '', lwd = 2)
lines(l.travel, col = high,  lwd = 2)
mtext(at = 0, line = 0.25, 'B', cex = 2)
mtext(side = 3, line = 2, expression('Mobility'), cex = 2)

axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(1,10.5, 20), cex.axis = 1.5)
axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(1,10.5, 20), cex.axis = 1.5)
mtext(side = 1, line = 6, expression('Intervention Coverage'), cex = 2)
axis(2, paste((seq(0, -15, by = -5)), '%', sep = ''), at = seq(0, -15, by = -5), las = 2, cex.axis = 1.5)
par(xpd = TRUE)}

}


#####################
#FIGURE 6A-C
#####################
{
  setwd('~/Desktop/MANUSCRIPT/FIGURE_6/')

  #LOADS IN PROPORTION OF CASES AVERTED
  load("norm_travel_test_list_no_travel.RData")
  high_whole_test <- norm_travel_test_list[[1]]
  whole_test <- norm_travel_test_list[[2]]
  low_whole_test <- norm_travel_test_list[[3]]
  high_vac_test <- norm_travel_test_list[[4]] 
  vac_test <- norm_travel_test_list[[5]] 
  low_vac_test <- norm_travel_test_list[[6]] 
  high_unvac_test <- norm_travel_test_list[[7]]
  unvac_test <- norm_travel_test_list[[8]]
  low_unvac_test <- norm_travel_test_list[[9]]
  
{
  
#FIGURE 6A
{ par(mfcol = c(3,2))
  par(pty = 's', mar = c(2, 0, 2,0),  xpd = FALSE, omi = c(1,0,1,0))
  plot(whole_test, type = 'l', ylim = c(0, 15), yaxt = 'n',lwd = 3,
       xlab = '', xaxt = 'n', ylab = '',
       main = "", cex.main = 2, col = black
  )
  mtext(at = 0, line = 0.25, 'A', cex = 2)
  lines(low_whole_test, lwd = 3,col = high)
  lines(high_whole_test,lwd = 3,  col = low)
  axis(2, labels = paste(seq(0, 15, by = 3), '%', sep = ''), at = seq(0, 15, by = 3), las = 2, cex.axis = 1.5)
  mtext(side= 3, line = 2 ,expression('No Mobility'), cex = 2)
  mtext(side= 2, line = 4.5 ,expression('Cases Averted'), cex = 2)
  abline(v= 10.5)}

#FIGURE 6B 
{  plot(vac_test, type = 'l', ylim = c(0, 15), yaxt = 'n',
       xlab = '', xaxt = 'n', ylab = '',lwd = 3,
       main = "", cex.main = 2, col = black
  )
  lines(c(low_vac_test[1:19], NA), lwd = 3,col = high)
  mtext(at = 0, line = 0.25, 'B', cex = 2)
  lines(c(NA,high_vac_test[1:19]), lwd = 3, col = low)
  axis(2, labels = paste(seq(0, 15, by = 3), '%', sep = ''), at = seq(0, 15, by = 3), las = 2, cex.axis = 1.5)
  abline(v= 10.5)
  mtext(side= 2, line = 4.5 ,expression('Cases Averted'), cex = 2)}
  
#FIGURE 6C
{  plot(high_unvac_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,lwd = 3,
       xlab = '', xaxt = 'n', ylab = '',
       main = "", cex.main = 2
  )
  lines(unvac_test,lwd = 3, col = black)
  lines(low_unvac_test, lwd = 3,col = high)
  mtext(at = 0, line = 0.25, 'C', cex = 2)
  axis(2, labels = paste(seq(0, 15, by = 3), '%', sep = ''), at = seq(0, 15, by = 3), las = 2, cex.axis = 1.5)
  abline(v= 10.5)
  axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(1,10.5, 20), cex.axis = 1.5)
  axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(1,10.5, 20), cex.axis = 1.5)
  mtext(side= 2, line = 4.5 ,expression('Cases Averted'), cex = 2)
  mtext(side= 1, line = 6 ,expression('Intervention'), cex = 2)
  }

}

}


#####################
#FIGURE 6D-F
#####################
{
  setwd('~/Desktop/MANUSCRIPT/FIGURE_6/')
  load("norm_travel_test_list_travel.RData")
  high_whole_test <- norm_travel_test_list[[1]]
  whole_test <- norm_travel_test_list[[2]]
  low_whole_test <- norm_travel_test_list[[3]]
  high_vac_test <- norm_travel_test_list[[4]] 
  vac_test <- norm_travel_test_list[[5]] 
  low_vac_test <- norm_travel_test_list[[6]] 
  high_unvac_test <- norm_travel_test_list[[7]]
  unvac_test <- norm_travel_test_list[[8]]
  low_unvac_test <- norm_travel_test_list[[9]]

  #FIGURE 6D
{plot(high_whole_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,lwd = 3,
     xlab = '', xaxt = 'n', ylab = '',
     main = "", cex.main = 2
)
lines(low_whole_test, lwd = 3,col = high)
lines(whole_test,lwd = 3, col = black)
text(x = 22.5, y = 7.5, expression('Total'), cex = 2, srt = 270, xpd = TRUE)
axis(2, labels = paste(seq(0, 15, by = 3), '%', sep = ''), at = seq(0, 15, by = 3), las = 2, cex.axis = 1.5)
mtext(at = 0, line = 0.25, 'D', cex = 2)
mtext(side= 3, line = 2 ,expression('Mobility'), cex = 2)
abline(v= 10.5)}

  #FIGURE 6E
{plot(high_vac_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,
     xlab = '', xaxt = 'n', ylab = '',lwd = 3,
     main = "", cex.main = 2
)
mtext(at = 0, line = 0.25, 'E', cex = 2)
lines(vac_test, lwd = 3,col = black)
lines(low_vac_test, lwd = 3,col = high)
text(x = 22.5, y = 7.5, expression('Vaccinated'), cex = 2, srt = 270, xpd = TRUE)
axis(2, labels = paste(seq(0, 15, by = 3), '%', sep = ''), at = seq(0, 15, by = 3), las = 2, cex.axis = 1.5)
abline(v= 10.5)
}

  #FIGURE 6F
{plot(high_unvac_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,lwd = 3,
     xlab = '', xaxt = 'n', ylab = '',
     main = "", cex.main = 2
)
lines(unvac_test,lwd = 3, col = black)
lines(low_unvac_test, lwd = 3,col = high)
text(x = 22.5, y = 7.5, expression('Unvaccinated'), cex = 2, srt = 270, xpd = TRUE)
mtext(at = 0, line = 0.25, 'F', cex = 2)
abline(v= 10.5)
axis(2, labels = paste(seq(0, 15, by = 3), '%', sep = ''), at = seq(0, 15, by = 3), las = 2, cex.axis = 1.5)
axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(1,10.5, 20), cex.axis = 1.5)
axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(1,10.5, 20), cex.axis = 1.5)
mtext(side = 1, line = 6, expression('Intervention'), cex = 2)}



}


#####################
#FIGURE 7A-C
#####################
{
  setwd('~/Desktop/MANUSCRIPT/FIGURE_7//')
  #load all cases averted
  load("norm_travel_test_list_skew_h.RData")
  high_whole_test <- norm_travel_test_list[[1]]
  whole_test <- norm_travel_test_list[[2]]
  low_whole_test <- norm_travel_test_list[[3]]
  high_vac_test <- norm_travel_test_list[[4]] 
  vac_test <- norm_travel_test_list[[5]] 
  low_vac_test <- norm_travel_test_list[[6]] 
  high_unvac_test <- norm_travel_test_list[[7]]
  unvac_test <- norm_travel_test_list[[8]]
  low_unvac_test <- norm_travel_test_list[[9]]
  
  {
   #FIGURE 7A
   {    par(mfcol = c(3,2))
    par(pty = 's', mar = c(2, 0, 2,0),  xpd = FALSE, omi = c(1,0,1,0))
    plot(high_whole_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,lwd = 3,
         xlab = '', xaxt = 'n', ylab = '',
         main = "", cex.main = 2
    )
    mtext(at = 0, line = 0.25, 'A', cex = 2)
    mtext(side= 3, line = 6.5 ,expression('Larger'), cex = 2)
    mtext(side= 3, line = 4.5 ,expression('Low-Transmission'), cex = 2)
    mtext(side= 3, line = 2.5 ,expression('Community'), cex = 2)
    abline(v = 10.5)
    lines(whole_test, lwd = 3, col = black)
    lines(low_whole_test, lwd = 3,col = high)
    axis(2, labels = paste(seq(0, 30, by = 3), '%', sep = ''), at = seq(0, 30, by = 3), las = 2, cex.axis = 1.5)
    mtext(side = 2, line = 4.5, expression('Cases Averted'), cex = 2)}
    
   #FIGURE 7B
   {    plot(high_vac_test, type = 'l', ylim = c(00,15), yaxt = 'n', col = low,
         xlab = '', xaxt = 'n', ylab = '',lwd = 3,
         main = "", cex.main = 2
    )
    mtext(at = 0, line = 0.25, 'B', cex = 2)
    lines(vac_test, lwd = 3, col = black)
    lines(low_vac_test,lwd = 3, col = high)
    abline(v = 10.5)
    axis(2, labels = paste(seq(0, 30, by = 3), '%', sep = ''), at = seq(0, 30, by = 3), las = 2, cex.axis = 1.5)
    mtext(side = 2, line = 4.5, expression('Cases Averted'), cex = 2)}
  
     #FIGURE 7C 
   { plot(high_unvac_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,lwd = 3,
         xlab = '', xaxt = 'n', ylab = '',
         main = "", cex.main = 2
    )
    mtext(at = 0, line = 0.25, 'C', cex = 2)
    lines(low_unvac_test,lwd = 3, col = high)
    lines(unvac_test, lwd = 3,col = black)
    abline(v = 10.5)
    axis(2, labels = paste(seq(0, 15, by = 3), '%', sep = ''), at = seq(0, 15, by = 3), las = 2, cex.axis = 1.5)
    mtext(side = 2, line = 4.5, expression('Cases Averted'), cex = 2)
    mtext(side = 1, line = 6, expression('Intervention'), cex = 2)
    axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(1,10.5, 20), cex.axis = 1.5)
    axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(1,10.5, 20), cex.axis = 1.5)}
  }
}


#####################
#FIGURE 7D-F
#####################
{ 

  load("norm_travel_test_list_skew_l.RData")
  high_whole_test <- norm_travel_test_list[[1]]
  whole_test <- norm_travel_test_list[[2]]
  low_whole_test <- norm_travel_test_list[[3]]
  high_vac_test <- norm_travel_test_list[[4]] 
  vac_test <- norm_travel_test_list[[5]] 
  low_vac_test <- norm_travel_test_list[[6]] 
  high_unvac_test <- norm_travel_test_list[[7]]
  unvac_test <- norm_travel_test_list[[8]]
  low_unvac_test <- norm_travel_test_list[[9]]
  {
    
    #FIGURE 7D
 {   plot(high_whole_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,lwd = 3,
         xlab = '', xaxt = 'n', ylab = '',
         main = "", cex.main = 2
    )
    mtext(at = 0, line = 0.25, 'D', cex = 2)
    mtext(side= 3, line = 6.5 ,expression('Larger'), cex = 2)
    mtext(side= 3, line = 4.5 ,expression('High-Transmission'), cex = 2)
    mtext(side= 3, line = 2.5 ,expression('Community'), cex = 2)
    text(x = 22.5, y = 7.5, expression('Total'), cex = 2, srt = 270, xpd = TRUE)
    axis(2, labels = paste(seq(0, 30, by = 3), '%', sep = ''), at = seq(0, 30, by = 3), las = 2, cex.axis = 1.5)
    lines(whole_test, lwd = 3, col = black)
    abline(v= 10.5)
    lines(low_whole_test, lwd = 3,col = high)}
   
    #FIGURE 7E
{    plot(high_vac_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,lwd = 3,
         xlab = '', xaxt = 'n', ylab = '',
         main = "", cex.main = 2
    )
    mtext(at = 0, line = 0.25, 'E', cex = 2)
    
    lines(vac_test, lwd = 3,col = black)
    lines(low_vac_test, lwd = 3,col = high)
    abline(v = 10.5)
    text(x = 22.5, y = 7.5, expression('Vaccinated'), cex = 2, srt = 270, xpd = TRUE)
    axis(2, labels = paste(seq(0, 30, by = 3), '%', sep = ''), at = seq(0, 30, by = 3), las = 2, cex.axis = 1.5)
    }

    #FIGURE 7F
{    plot(high_unvac_test, type = 'l', ylim = c(0, 15), yaxt = 'n', col = low,lwd = 3,
         xlab = '', xaxt = 'n', ylab = '',
         main = "", cex.main = 2
    )
    mtext(at = 0, line = 0.25, 'F', cex = 2)
    
    lines(unvac_test, lwd = 3,col = black)
    lines(low_unvac_test, lwd = 3,col = high)
    abline(v= 10.5)
    text(x = 22.5, y = 7.5, expression('Unvaccinated'), cex = 2, srt = 270, xpd = TRUE)
    axis(2, labels = paste(seq(0, 15, by = 3), '%', sep = ''), at = seq(0, 15, by = 3), las = 2, cex.axis = 1.5)
    mtext(side = 1, line = 6, expression('Intervention'), cex = 2)
    axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(1,10.5, 20), cex.axis = 1.5)
    axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(1,10.5, 20), cex.axis = 1.5)}

  }
}


#####################
#FIGURE 8A-C
#####################
{
	setwd('~/Desktop/MANUSCRIPT/FIGURE_8/')
  load("norm_travel_test_list_spec.RData")
  high_whole_test <- norm_travel_test_list[[1]]
  whole_test <- norm_travel_test_list[[2]]
  low_whole_test <- norm_travel_test_list[[3]]
  high_vac_test <- norm_travel_test_list[[4]] 
  vac_test <- norm_travel_test_list[[5]] 
  low_vac_test <- norm_travel_test_list[[6]] 
  high_unvac_test <- norm_travel_test_list[[7]]
  unvac_test <- norm_travel_test_list[[8]]
  low_unvac_test <- norm_travel_test_list[[9]]
	

colfunc <- colorRampPalette(c('lightskyblue1', low))
par(mfcol = c(3,2))
par(pty = 's', mar = c(2, 0, 2,0),  xpd = FALSE, omi = c(1,0,1,0))

#FIGURE 8A
{
  plot(high_whole_test[1:20], type = 'l', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex.lab = 2,
       ylim = c(-5, 20), col = colfunc(6)[1], lwd= 3,
       main = '', cex.main = 2)
  mtext(at = 0, line = 0.25, 'A', cex = 2)
  mtext(side = 3, line = 2.5, 'Low-Transmission\nCommunity', cex = 2)
  mtext(side = 2, line = 4.5, expression('Cases Averted'), cex = 2)
    for(i in 1:21){
    index <- 20 * i + c(1:20)
    lines(high_whole_test[index], col = colfunc(21)[i + 1], lwd= 3)
  }
  axis(2, paste((seq(-5, 20, by  = 5)), '%', sep = ''), at = seq(-5, 20, by = 5), las = 2, cex.axis = 1.5)
  abline(h = 0)
  abline(v = 10.5)
}

#FIGURE 8B
{
  plot(high_vac_test[1:20], type = 'l', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex.lab = 2,
       ylim = c(-5, 20), col = colfunc(6)[1], lwd = 3,
       main = '', cex.main = 2)
  mtext(at = 0, line = 0.25, 'B', cex = 2)
   mtext(side = 2, line = 4.5, expression('Cases Averted'), cex = 2)
 for(i in 1:21){
    index <- 20 * i + c(1:20)
    lines(high_vac_test[index] , col = colfunc(21)[i + 1], lwd = 3)
  }
  axis(2, paste((seq(-5, 20, by  = 5)), '%', sep =''), at = seq(-5, 20, by  = 5), las = 2, cex.axis = 1.5)
  abline(v = 10.5)
  abline(h = 0)

}

#FIGURE 8C
{    
  plot(high_unvac_test[1:20], type = 'l', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex.lab = 2,
      ylim = c(-5, 20), col = colfunc(6)[1], lwd = 3,
      main = '', cex.main = 2)
  mtext(at = 0, line = 0.25, 'C', cex = 2)
  axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(1,10.5, 20), cex.axis = 1.5)
  axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(1,10.5, 20), cex.axis = 1.5)
  mtext(side = 2, line = 4.5, expression('Cases Averted'), cex = 2)
  for(i in 1:21){
    index <- 20 * i + c(1:20)
    lines(high_unvac_test[index], col = colfunc(21)[i + 1], lwd = 3)
  }
  abline(v = 10.5)
  axis(2, paste((seq(-5, 20, by  =5)), '%', sep = ''), at = seq(-5, 20, by = 5), las = 2, cex.axis = 1.5)
   abline(h = 0)
  mtext(side = 1, line = 6, expression('Intervention'), cex = 2)
  
  } 

}


#####################
#FIGURE 8D-F
#####################
{
  colfunc <- colorRampPalette(c('lightgoldenrod', high))
  
#FIGURE 8D
  {
    plot(low_whole_test[1:20], type = 'l', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex.lab = 1.5,
         ylim = c(-5, 20), col = colfunc(6)[1], lwd= 3,
         main = '', cex.main = 2)
    mtext(side = 3, line = 2.5, 'High-Transmission\nCommunity',cex = 2)
    text(x = 22.5, y = 7.5, expression('Total'), cex = 2, srt = 270, xpd = TRUE)
    axis(2, paste((seq(-5, 20, by  = 5)), '%', sep = ''), at = seq(-5, 20, by = 5), las = 2, cex.axis = 1.5)
    mtext(at = 0, line = 0.25, 'D', cex = 2)
    for(i in 1:21){
      index <- 20 * i + c(1:20)
      lines(low_whole_test[index], col = colfunc(21)[i + 1], lwd= 3)
    }
    abline(h = 0)
    abline(v = 10.5)
    colfunc <- colorRampPalette(c('lightskyblue1', low))
    colfunc <- colorRampPalette(c('lightgoldenrod', high))
    
  }
  
#FIGURE 8E
  {
    plot(low_vac_test[1:20], type = 'l', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex.lab = 2,
         ylim = c(-5, 20), col = colfunc(6)[1], lwd = 3,
         main = '', cex.main = 2)
    text(x = 22.5, y = 7.5, expression('Vaccinated'), cex = 2, srt = 270, xpd = TRUE)
    axis(2, paste((seq(-5, 20, by  = 5)), '%', sep = ''), at = seq(-5, 20, by = 5), las = 2, cex.axis = 1.5)
    mtext(at = 0, line = 0.25, 'E', cex = 2)
    for(i in 1:21){
      index <- 20 * i + c(1:20)
      lines(low_vac_test[index], col = colfunc(21)[i + 1], lwd = 3)
    }
    abline(v = 10.5)
     abline(h = 0)

  }
  
#FIGURE 8F
  { 
    plot(low_unvac_test[1:20], type = 'l', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex.lab = 2,
         ylim = c(-5, 20), col = colfunc(6)[1], lwd = 3,
         main = '', cex.main = 2)
    text(x = 22.5, y = 7.5, expression('Unvaccinated'), cex = 2, srt = 270, xpd = TRUE)
    axis(2, paste((seq(-5, 20, by  = 5)), '%', sep = ''), at = seq(-5, 20, by = 5), las = 2, cex.axis = 1.5)
    mtext(at = 0, line = 0.25, 'F', cex = 2)
    axis(side = 1, line = 0.5, labels= c('100% HT', '50% HT', '0% HT'), at = c(1,10.5, 20), cex.axis = 1.5)
    axis(side = 1, line = 3, labels= c('0% LT', '50% LT', '100% LT'), at = c(1,10.5, 20), cex.axis = 1.5)    # mtext(side = 3, line = 1, expression(bold('Unvaccinated')), cex = 2)
    for(i in 1:21){
      index <- 20 * i + c(1:20)
      lines(low_unvac_test[index], col = colfunc(21)[i + 1], lwd = 3)
    }
      abline(h = 0)
     abline(v = 10.5)
     mtext(side = 1, line = 6, expression('Intervention'), cex = 2)

  } 

}

