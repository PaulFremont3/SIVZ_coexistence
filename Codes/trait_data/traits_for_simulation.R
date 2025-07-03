library('nnet')
library('mgcv')

# load data
u=read.table('edwards_2018.txt',sep = '\t', header = T)

to_predict='BS'
model_type='nn'
if (to_predict=='BS'){
  df_full=data.frame(y=log10(u$Burst.size), log10Vhost=log10(u$Host.cell.volume..cubic.microns.), RVirus=u$Virus.particle.size..nm., VirusType=u$Virus.type, HostType=u$Host.taxon, Temperature=u$Temperature..ºC., env=u$Environment)
  lab_ax=paste('log10(',to_predict,')', sep='')
} else if (to_predict=='LP'){
  df_full=data.frame(y=log10(u$Latent.period..hr.), log10Vhost=log10(u$Host.cell.volume..cubic.microns.), RVirus=u$Virus.particle.size..nm., VirusType=u$Virus.type, HostType=u$Host.taxon, Temperature=u$Temperature..ºC., env=u$Environment)
  lab_ax=paste('log10(',to_predict,')', sep='')
}

# remove under represented data
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

# load models
nn_mod=readRDS('BS_nn_result.rds')
nn_mod_BS=nn_mod[[2]]

gam_mod=readRDS('BS_gam_result.rds')
gam_mod_BS=gam_mod[[2]]

lm_mod=readRDS('BS_lm_result.rds')
lm_mod_BS=lm_mod[[2]]

nn_mod=readRDS('LP_nn_result.rds')
nn_mod_LP=nn_mod[[2]]

lm_mod=readRDS('LP_lm_result.rds')
lm_mod_LP=lm_mod[[2]]

gam_mod=readRDS('LP_gam_result.rds')
gam_mod_LP=gam_mod[[2]]

min_V_diatom=min(u$Host.cell.volume..cubic.microns.[u$Host.taxon=='Diatom'])

# host volumes
Host_volumes=as.numeric(read.table('Vs_5.txt'))
Host_types=NULL
for (i in 1:100){
  if (Host_volumes[i]<min_V_diatom){
    Host_types=c(Host_types,'Eukaryote')
  } else{
    Host_types=c(Host_types,'Diatom')
  }
 
}
Host_types=c(Host_types, rep('Eukaryote', 100), rep('Cyanobacteria', 200))
Host_types_true=c(rep('Diatom', 100,), rep('Eukaryote', 100), rep('Cyanobacteria', 200))

# virus types
Virus_type=c(rep('ssDNA (partial ds)', 100), rep('dsDNA', 300))

# virus quotas
Qn_virus=function(r){
  Na=6.022140857*1e23
  Qn=(1e6/Na)*(16*(r-2.5)^3+36*(7.5*r^2-18.75*r+15.63))
  return(Qn)
}

Qc_virus=function(r){
  Na=6.022140857*1e23
  Qn=(1e6/Na)*(41*(r-2.5)^3+130*(7.5*r^2-18.75*r+15.63))
  return(Qn)
}

# phytoplankton quotas
Qn_algae=function(i, data){
  V=10^data$log10Vhost[i]
  type=data$HostType[i]
  if (type=='Diatom'){
    Qc=10^(-0.541 + 0.811*log10(V))
    Qc_micro=Qc*1e-6
  } else if (type=='Cyanobacteria'){
    dc=470
    Qc=dc*V
    Qc_micro=Qc*1e-9
  } else if (type=='Eukaryote'){
    Qc=10^(-0.665 + 0.939*log10(V))
    Qc_micro=Qc*1e-6
  }
  Qc_micromol=Qc_micro/12
  Qn_micromol=Qc_micromol*16/106
  return(Qn_micromol)
}

# virus radiuses
vr_d=round(mean(df_full$RVirus[df_full$HostType=='Diatom']),0)
vr_euk=round(mean(df_full$RVirus[df_full$HostType=='Eukaryote']),0)
vr_cya=round(mean(df_full$RVirus[df_full$HostType=='Cyanobacteria']),0)
virus_radius=c(rep(vr_d, 100),
               rep(vr_euk, 100), rep(vr_cya, 200) )

# data to predict (the 400 phytoplankton types)
new_data=data.frame(log10Vhost=log10(Host_volumes), RVirus=virus_radius,HostType=Host_types, VirusType=Virus_type, stringsAsFactors = F)
new_data_gam=new_data

Qps=sapply(1:dim(new_data)[1], Qn_algae, data=new_data)

for (v in colnames(df_full)[(9):(dim(df_full)[2])]){
  new_data[[v]]=NA
}


for (i in 1:dim(new_data)[1]){
  vt=new_data$VirusType[i]
  ht=new_data$HostType[i]
  for (v in colnames(df_full)[(9):(dim(df_full)[2])]){
    if (v==vt){
      new_data[[v]][i]=1
    } else if (v==ht){
      new_data[[v]][i]=1
    } else{
      new_data[[v]][i]=0
    }
  }
}
if (model_type=='nn'){
  new_data$HostType=NULL
  new_data$VirusType=NULL
}

# generate predictions
pred_BS=10^stats::predict(nn_mod_BS, newdata=new_data)[,1]
pred_BS=round(pred_BS / 5) * 5
pred_LPs=NULL
for (i in 1:100){
  pred_LP0=10^stats::predict(nn_mod_LP, newdata=new_data)[,1]
  pred_LPs=rbind(pred_LPs, pred_LP0)
}
find_good=function(x){
  v=x[x >= 4] |> table() |> which.max() |> names() |> as.numeric()
  return(v)
}
#print(pred_LPs)
pred_LP=apply(pred_LPs, 2, find_good)


pred_BS_gam=as.numeric(10^predict.gam(gam_mod_BS, new_data_gam))
pred_BS_gam=round(pred_BS_gam / 5) * 5

pred_BS_lm=as.numeric(10^stats::predict(lm_mod_BS, newdata = new_data_gam))
pred_BS_lm=round(pred_BS_lm / 5) * 5


pred_LP_lm=as.numeric(10^stats::predict(lm_mod_LP, newdata = new_data_gam))

pred_LP_gam=as.numeric(10^predict.gam(gam_mod_LP, new_data_gam))

mu_max=as.numeric(read.table('mumax_dutkiewicz_5.txt'))

ward_mumax = function(V){
  mumax=0.1128*V^0.84/(0.1504*V^0.5 + 0.024*V^1.1)
  return(mumax)
}
mu_max_bis= sapply(Host_volumes, ward_mumax)

dut2020_mumax=function(i, Vs){
  V=Vs[i]
  if (i %in% c(1:100)){ # diatom
    mum= 3.9*V^-0.08
  } else if (i %in% c(101:200)){ # eukaryotes
    mum= 1.4*V^-0.08
  } else{
    mum= 0.9*V^0.08
  }
}

mu_max_ter= sapply(1:400, dut2020_mumax, Vs=Host_volumes)

types=c('Diatom', 'Euk', 'Syn', "Proch")
for (i in 0:3){
  ind_min_Qp=which.min(Qps[(i*100+1):(i*100+100)])
  if (i!=0){
    bs=pred_BS[(i*100+1):(i*100+100)][ind_min_Qp]
  } else{
    bs=pred_BS_gam[(i*100+1):(i*100+100)][ind_min_Qp]
  }
  bs_lm=pred_BS_lm[(i*100+1):(i*100+100)][ind_min_Qp]
  bs_gam=pred_BS_gam[(i*100+1):(i*100+100)][ind_min_Qp]

  lp=pred_LP[(i*100+1):(i*100+100)][ind_min_Qp]
  lp_lm=pred_LP_lm[(i*100+1):(i*100+100)][ind_min_Qp]
  lp_gam=pred_LP_gam[(i*100+1):(i*100+100)][ind_min_Qp]
  
  Qp=Qps[(i*100+1):(i*100+100)][ind_min_Qp]
  
  print(types[i+1])
  Hv=Host_volumes[(i*100+1):(i*100+100)][ind_min_Qp]
  print('host volume:')
  print(Hv)
  print('host radius')
  Hr=(3*Hv/(4*pi))^(1/3)
  print(Hr)
  print('burst size:')
  print(bs)
  print(bs_lm)
  print(bs_gam)
  
  
  print('latent period')
  print(lp/24)
  print(lp_lm/24)
  print(lp_gam/24)
  print(' ')
  print('Qp')
  print(Qp)
  print(' ')
  
  print('mu')
  print(max(mu_max[(i*100+1):(i*100+100)]))
  print(' ')
  
  print('mu_ward')
  print(mu_max_bis[(i*100+1):(i*100+100)][ind_min_Qp])
  
  print('mu_dut')
  print(mu_max_ter[(i*100+1):(i*100+100)][ind_min_Qp])
  print(' ')
}

final_preds_BS= c(pred_BS_gam[1:100], pred_BS[101:400])
final_preds_LP= pred_LP

final_preds_BS_lm= pred_BS_lm
final_preds_LP_lm= pred_LP_lm

colos=c('red', 'blue', 'green')
taxos=c('Cyanobacteria', 'Eukaryote', 'Diatom')
new_data_gam$col=colos[match(new_data_gam$HostType, taxos)]


log10Tck <- function(side, type){
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceiling(lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}

pdf('trait_values_size_structured_foodweb.pdf')
plot(Host_volumes, final_preds_BS, log='xy',col=new_data_gam$col, pch=19, axes=FALSE,frame.plot=TRUE , xlab='V', ylab='Model burst size')
axis(1, at=log10Tck('x','major'), tcl= 0.2) # bottom
axis(2, at=log10Tck('y','major'), tcl= 0.2) 
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA)
axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA)


plot(Host_volumes, final_preds_BS_lm, log='xy',col=new_data_gam$col, pch=19, axes=FALSE,frame.plot=TRUE , xlab='V', ylab='Model burst size')
axis(1, at=log10Tck('x','major'), tcl= 0.2) # bottom
axis(2, at=log10Tck('y','major'), tcl= 0.2) 
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA)
axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA)

plot(Host_volumes, pred_LP/24, log='x',col=new_data_gam$col, pch=19,  axes=FALSE,frame.plot=TRUE , xlab='V', ylab='Model latent period')
axis(1, at=log10Tck('x','major'), tcl= 0.2) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA)
axis(2)

plot(Host_volumes, final_preds_LP_lm/24, log='x',col=new_data_gam$col, pch=19 ,axes=FALSE,frame.plot=TRUE  , xlab='V', ylab='Model latent period')
axis(1, at=log10Tck('x','major'), tcl= 0.2) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA)
axis(2)
dev.off()

writeLines(as.character(final_preds_BS), 'model_burst_size_nn-gam.txt', sep = ' ')
writeLines(as.character(final_preds_BS_lm), 'model_burst_size_lm.txt',sep = ' ')

writeLines(as.character(final_preds_LP/24), 'model_latent_period_nn-gam.txt',sep = ' ')
writeLines(as.character(final_preds_LP_lm/24), 'model_latent_period_lm.txt',  sep = ' ')
