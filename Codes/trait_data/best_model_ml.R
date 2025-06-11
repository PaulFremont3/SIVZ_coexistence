source('optimisation_gam.R')
source('optimisation_nn.R')
source('optimisation_nnet.R')
source('optimisation_lm.R')
source('optimisation_rf.R')
library('mgcv')
library('nnet')
library('neuralnet')
library('Metrics')
u=read.table('edwards_2018.txt',sep = '\t', header = T)

to_predict=commandArgs(trailingOnly = T)[1]
model_type=commandArgs(trailingOnly = T)[2]
if (to_predict=='BS'){
	df_full=data.frame(y=log10(u$Burst.size), log10Vhost=log10(u$Host.cell.volume..cubic.microns.), RVirus=u$Virus.particle.size..nm., VirusType=u$Virus.type, HostType=u$Host.taxon, Temperature=u$Temperature..ºC., env=u$Environment)
	lab_ax=paste('log10(',to_predict,')', sep='')
} else if (to_predict=='LP'){
	df_full=data.frame(y=log10(u$Latent.period..hr.), log10Vhost=log10(u$Host.cell.volume..cubic.microns.), RVirus=u$Virus.particle.size..nm., VirusType=u$Virus.type, HostType=u$Host.taxon, Temperature=u$Temperature..ºC., env=u$Environment)
	lab_ax=paste('log10(',to_predict,')', sep='')
}
#df_full=df_full[df_full$env=='Marine',]
df_full$HostType[df_full$HostType %in% c("Prasinophyte"  , "Haptophyte",     "Dinoflagellate" )]="Eukaryote"

is_na_rows=apply(df_full, 1, FUN = function(x){u=sum(is.na(x)); return(ifelse(u==0, T, F))})

df_full=df_full[is_na_rows,]

variables_to_check=c('HostType', 'VirusType')

to_rem=NULL
for (var in variables_to_check){
  u_var=unique(df_full[[var]])
  for (uv in u_var){
    su=sum(df_full[[var]]==uv)
    if (su<4){
      to_rem=c(to_rem,which(df_full[[var]]==uv))
      df_full=df_full[df_full[[var]]!=uv,]
    }
  }
}
to_rem=NULL

df_full$pairs=paste(df_full$VirusType, df_full$HostType, sep='_')
pairs=unique(df_full$pairs)
for (phv in pairs){
  su=sum(df_full$pairs==phv)
  if (su<5){
    df_full=df_full[df_full$pairs!=phv,]
  }
}

if (model_type %in% c('nn', 'nnet')){
	cat_vars=c('HostType', 'VirusType')
	for (cvar in cat_vars){
  	flags = data.frame(Reduce(cbind,lapply(unique(df_full[[cvar]]), function(x){(df_full[[cvar]] == x)*1})
  	))
  	names(flags) = unique(df_full[[cvar]])
  	df_full = cbind(df_full, flags)
	}

	variables=c(2:3, c(9:dim(df_full)[2]))
	predictors = paste(colnames(df_full)[variables], collapse = " + ")
	model_formula=paste(colnames(df_full)[1], '~', predictors)
}


print(dim(df_full))

colors=c('deepskyblue',  'red',  'green')
unique_Htype=unique(df_full$HostType)
df_full$color=colors[match(df_full$HostType,unique_Htype)]

id=1
if (model_type=='gam'){
	gamGrid <-  expand.grid(nsplines = seq(from = -5, to = 3, by = 1))
	variables=2:5
	score_list <- best_models_gam(id)
} else if (model_type=='nn'){
	nnetGrid <-  expand.grid(size = seq(from = 1, to = 10, by = 1),
                         decay = c(seq(1, 10, by= 2) %o% 10^(-5:-2)),
                         mxit = c(200, 500, 1000))
	score_list <- best_models_nn(id, df_full)
} else if (model_type=='lm'){
	gamGrid <-  expand.grid(nsplines = seq(from = 0, to = 0, by = 1))
        variables=2:5
        score_list <- best_models_lm(1)
} else if (model_type=='nnet'){
	nnetGrid <-  expand.grid(hidden = seq(from = 1, to = 4, by = 1),
                         decay = c(10^(-2:-3)),
                         mxit = c(1000000, 5000000), size=seq(from = 1, to = 7, by = 1))
	score_list <- best_models_nnet(id,  model_formula, df_full)
} else if (model_type=='rf'){
	rfGrid <-  expand.grid(mtry = seq(from = 1, to = 7, by = 1),
                       ntree = c(seq(1, 10, by= 2) %o% 10^(2:3)))
	variables=2:5
        score_list <- best_models_rf(id=id, df_full)
}


mod=score_list[[2]]

saveRDS(score_list, paste(to_predict,'_',model_type, '_result.rds', sep=''))

colors=c('deepskyblue',  'red',  'green')
unique_Htype=unique(df_full$HostType)
df_full$color=colors[match(df_full$HostType,unique_Htype)]
#plot(df_full$y, mod$fitted.values, col=df_full$color, pch=19 )

#variables=2:6
#score_list <- best_models_gam(1)

mod=score_list[[2]]


pdf(paste('pred_vs_data_',to_predict,'_',model_type,'.pdf', sep=''))
if (model_type=='rf'){
  fitted=mod$predicted
} else{
  fitted=mod$fitted.values
}
mi=min(c(df_full$y, fitted))
mx=max(c(df_full$y, fitted))
plot(df_full$y, fitted, col=df_full$color, pch=19,xlim=c(mi,mx), ylim=c(mi,mx), xlab=lab_ax, ylab=paste(lab_ax,'predicted', sep=' ') ) 
points(seq(mi, mx, by=(mx-mi)/10), seq(mi, mx, by=(mx-mi)/10), type='l', lty=2) 
co=cor.test(df_full$y, fitted)
rms=rmse(df_full$y, fitted)
title(paste('r=',round(co$estimate,2),',pval=',co$p.value, 'rmse=', rms ))

selec=df_full$color!='green'
x=df_full$y[selec]
y=fitted[selec]

mi=min(c(x, y))
mx=max(c(x, y))
plot(x, y, col=df_full$color[selec], pch=19,xlim=c(mi,mx), ylim=c(mi,mx), xlab=lab_ax, ylab=paste(lab_ax,'predicted', sep=' ') )
points(seq(mi, mx, by=(mx-mi)/10), seq(mi, mx, by=(mx-mi)/10), type='l', lty=2)
co=cor.test(x,y)
rms=rmse(x,y)
title(paste('r=',round(co$estimate,2),',pval=',co$p.value, 'rmse=', rms ))

dev.off()
#print(mean(abs(mod$residuals)))
#print(mean(abs(gam_model$residuals)))



df_full$pairs=paste(df_full$VirusType, df_full$HostType, sep='_')

pdf(paste('model_vs_data_',to_predict,'_host_volume_',model_type,'.pdf', sep=''))
pairs=unique(df_full$pairs)
all_pred_data=rep(list(NULL), length(pairs))
all_data=rep(list(NULL), length(pairs))
#if (to_predict=='BS')
col_grads=list('dsDNA_Eukaryote'=colorRampPalette(colors = c("lightblue",  'blue'))(100), 'dsDNA_Cyanobacteria'=colorRampPalette(colors = c("pink",  'red'))(100), 'ssDNA (partial ds)_Diatom'=colorRampPalette(colors = c("lightgreen",  'green'))(100), 'ssRNA_Diatom'=colorRampPalette(colors = c("mediumpurple1",  'mediumpurple4'))(100) ,
               colorRampPalette(colors = c("grey",  'black'))(100))
#names(col_grads)=pairs
for (phv in pairs){
  vt=strsplit(phv, split = '_')[[1]][1]
  ht=strsplit(phv, split = '_')[[1]][2]
  sel=df_full$HostType==ht & df_full$VirusType==vt
  min_host_size=min(df_full$log10Vhost[sel])
  max_host_size=max(df_full$log10Vhost[sel])
  
  if ( min_host_size==max_host_size){
    min_host_size=min_host_size-0.1
    max_host_size=max_host_size+0.1
  }
  
  
  lin_seq_hostsize=seq(round(min_host_size,1), round(max_host_size,1), by=0.1)
  
  min_virus_size=min(df_full$RVirus[sel])
  max_virus_size=max(df_full$RVirus[sel])
  
  if (min_virus_size==max_virus_size){
    min_virus_size=min_virus_size-2.5
    max_virus_size=max_virus_size+2.5
  }
  
  lin_seq_virus=seq(round(min_virus_size,1), round(max_virus_size,1), by=5)
  
  
  grid=expand.grid(log10Vhost=lin_seq_hostsize, RVirus=lin_seq_virus,HostType=ht, VirusType=vt, stringsAsFactors = F )
  
  if (model_type=='gam'){
    pred_bs=predict.gam(mod, newdata = grid)
  } else if (model_type=='nn'){
    grid=expand.grid(log10Vhost=lin_seq_hostsize, RVirus=lin_seq_virus )
    for (v in colnames(df_full)[(9):(dim(df_full)[2]-1)]){
        if (v==vt){
            grid[[v]]=1
        } else if (v==ht){
            grid[[v]]=1
        } else{
            grid[[v]]=0
        }
    } 

    pred_bs=stats::predict(mod, newdata = grid)
  } else if (model_type=='lm'){
    pred_bs=stats::predict(mod, newdata = grid)
  }  else if (model_type=='rf'){
    print(grid)
    pred_bs=stats::predict(mod, newdata = grid, type='response')
  } 
  #pred_bs=matrix(pred_bs, nrow=length(lin_seq_hostsize), ncol=length(lin_seq_virus))
  
  col_virus_size=(grid$RVirus-min(grid$RVirus))*99/(max(grid$RVirus)-min(grid$RVirus))+1
  col_grad=colorRampPalette(colors = c("yellow", "orange", 'red'))(100)
  
  mi = min(pred_bs, df_full$y[sel])
  mx = max(pred_bs, df_full$y[sel])
  print(mi)
  print(mx)
  plot(grid$log10Vhost, pred_bs,  col=col_grad[col_virus_size], pch=19, main=phv, ylim=c(mi,mx))
  
  all_pred_data[[phv]]=list(pred_bs, grid$log10Vhost, as.integer(col_virus_size))
  
  colo=(df_full$RVirus[sel]-min(grid$RVirus))*99/(max(grid$RVirus)-min(grid$RVirus))+1
  points(df_full$log10Vhost[sel], df_full$y[sel], col=col_grad[colo], pch=8, cex=2)
  all_data[[phv]]=list(df_full$y[sel], df_full$log10Vhost[sel], as.integer(colo))
}
i=1
for (phv in pairs){
  if (i==1){
    co=col_grads[[phv]]
    plot(all_pred_data[[phv]][[2]], all_pred_data[[phv]][[1]], pch=19, ylim=c(min(df_full$y), max(df_full$y)), xlim= c(min(df_full$log10Vhost),max(df_full$log10Vhost))   , col=co[all_pred_data[[phv]][[3]]], ylab=lab_ax, xlab="log10(Vhost)" )
    #points(all_data[[phv]][[2]], all_data[[phv]][[1]], pch=8, cex=2, col=col_grads[[phv]][all_data[[phv]][[3]]])
  } else{
    points(all_pred_data[[phv]][[2]], all_pred_data[[phv]][[1]], pch=19, col=col_grads[[phv]][all_pred_data[[phv]][[3]]])
    #points(all_data[[phv]][[2]], all_data[[phv]][[1]], pch=8, cex=2, col=col_grads[[phv]][all_data[[phv]][[3]]])
  }
  i=i+1 
}
for (phv in pairs){
  points(all_data[[phv]][[2]], all_data[[phv]][[1]], pch=8, cex=2, col=col_grads[[phv]][all_data[[phv]][[3]]])
}
dev.off()


saveRDS(mod,paste(to_predict,'_',model_type, '_model.rds', sep=''))
