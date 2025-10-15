data=read.table('carlson_data_cyanobacteria.txt', sep='\t', header = T)
data_v=read.table('carlson_data_virus_and_percentage_infected.txt', sep='\t', header = T)
data=data[data$Latitude<33,]
data_v=data_v[data_v$Latitude<33,]

pdf('hist_abundances_Syn_Pro_virus_percentage_infected.pdf')
hist(data$Prochlorococcus..cells.mL.[data$Prochlorococcus..cells.mL.!=0]*1000, breaks=20, col=scales::alpha('green', 0.5),
     xlab='Cell concentrataion (ind.L-1)', main='', xlim=c(min(data$Prochlorococcus..cells.mL.[data$Prochlorococcus..cells.mL.!=0]*1000),
                                                               max(data$Prochlorococcus..cells.mL.[data$Prochlorococcus..cells.mL.!=0]*1000)))
box()
abline(v=mean(data$Prochlorococcus..cells.mL.[data$Prochlorococcus..cells.mL.!=0]*1000),col='green' ,lwd=4)
hist(data$Synechococcus..cells.mL.[data$Synechococcus..cells.mL.!=0]*1000, breaks=20, col=scales::alpha('blue', 0.5),
     xlab='Cell concentrataion (ind.L-1)', main='')
abline(v=mean(data$Synechococcus..cells.mL.[data$Synechococcus..cells.mL.!=0]*1000),col='blue' ,lwd=4)
box()


hist(data_v$Tot.V*1000, breaks=20, col=scales::alpha('red', 0.5),
     xlab='Cell concentrataion (ind.L-1)', main='')
abline(v=mean(data_v$Tot.V*1000),col='red' ,lwd=4)
box()

hist(data_v$X.I.Syn.Tot, breaks=20, col=scales::alpha('darkblue', 0.5),
     xlab='% Infected', main='')
abline(v=mean(data_v$X.I.Syn.Tot, na.rm=T),col='darkgreen' ,lwd=4)
box()

hist(data_v$X.I.Pro.Tot, breaks=20, col=scales::alpha('darkgreen', 0.5),
     xlab='% Infected', main='')
abline(v=mean(data_v$X.I.Pro.Tot, na.rm=T),col='darkblue' ,lwd=4)
box()
dev.off()
