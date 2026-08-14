setwd(main_dir)

F_mean=0.2
M=0.1
c=1
VB_func<-function(Linf,k,t0,age)
{
  Lt = Linf*(1-exp(-k*(age-t0))) # po is the length of age 0 fish.
  return(Lt)
}

mat_func<-function(L50,L95,length)
{
  b1=log(0.95/0.05)/(L95-L50)
  bo = -L50*b1
  logit_pt =  bo+b1*length
  matp= exp(logit_pt)/(1+exp(logit_pt))
  return(matp)
}

#---------------------------------------------------------------------------------
# define biological and fishing variables

nyear<-100*12 # number of years
rec.age<-1 # the age of recruitment to fishery
first.year<-2024*12 # the first year of projection

# biological variables
nage<-20 # number of ages

ages<-c(rec.age:20) # list of ages in the simulation
years<-c(first.year:(first.year+nyear-1)) # list of years in the simulation

len_mid<-seq(60,1200,50)
len_border<-seq(44,1216,50)
len_border[1]=-Inf
len_border[length(len_border)]=Inf
nlen<-length(len_mid)
len_lower=len_border[1:nlen]
len_upper=len_border[2:(nlen+1)]


#h=8.04/1E6; j=2.38
#ages<-c(rec.age:20)
#W_at_age = h*(ages^j)
W_at_age=85.4*exp(0.009*(ages*30))/1e6
wgt=W_at_age
#plot(ages,W_at_age)
# natural mortality

vbk=0.065
Linf=903
t0=0.6

LatA=VB_func(Linf,vbk,t0,ages)
#plot(ages,LatA)

# weigth at length
a=4.51E-6/1E6; b=3.31389
W_at_len = a*(len_mid^b)

#plot(len_mid,W_at_len,type="l")

# maturation at length
mat_L50=452; mat_L95=583
mat_L=mat_func(mat_L50,mat_L95,len_mid)
# maturation at age
mat_A50=13.6; mat_A95=15.9###16.9 13.6
mat=mat_func(mat_A50,mat_A95,ages)
mat=c(0,
      0,
      0,
      4.60E-06,
      1.65E-05,
      0.0200214013,
      0.150769447,
      0.202762431,
      0.259866573,
      0.234606434,
      0.214223726,
      0.316887115,
      0.525293092,
      0.557205018,
      0.655742263,
      0.757290931,
      0.85643436,
      0.85900625,
      0.959723555,
      0.999723555
)



cv_L=0.2 # the coefficient of variance in length at age at the beginning year
cv_inc=0.1 # the coefficient of variance in growth increment

std_logR=0.5 # standard deviation of recruitment
std_logN0=1 # standar deviation of initial number at age

# SR relationship (BH model) aS/(1+(bS)^c)
#alpha=2e8/(29941/1e6)
#beta=1/(29941/1e6)


alpha  <- 1e5 
beta=1e-4


# age-length tansition matrix
pla=matrix(NA,nrow=nage,ncol=nlen) # transfer age to length
ml=VB_func(Linf,vbk,t0,ages)
sl=cv_L*ml
for(i in 1:nlen){
  pla[,i]=pnorm(len_border[i+1],ml,sl)-pnorm(len_border[i],ml,sl)
}

# survey variables

std_SN=0.01 # survey measurement error

q_surv_L50=208 # survey catchability at length
q_surv_L95=243
q_surv=mat_func(q_surv_L50,q_surv_L95,len_mid)
#plot(len_mid,q_surv,type="l",xlim=c(0,700))
# fishing variable

sel_A50=6.9# fishing selectivity at length
sel_A95=8.6
sel=mat_func(sel_A50,sel_A95,ages)

#plot(ages,sel,type="l",xlim=c(0,25))
# sel=rep(1,nage)

#-------------------------------------------------------------------------------------------------
# start the simulation
for(iter in 1:100){
  
  
  iter_folder <- paste0("sim_iteration_", iter)
  dir.create(iter_folder, showWarnings = FALSE)
  
  
  set.seed(iter)
  F_yr <- F_mean * exp(arima.sim(list(order=c(1,0,0), ar=0.75),
                                 n = nyear) * 0.2)
  F_at_age <- matrix(NA, nrow=nyear, ncol=nage)
  Z_at_age <- matrix(NA, nrow=nyear, ncol=nage)
  M_at_age <-matrix(M,nrow=nyear,ncol=nage) 
  for(t in 1:nyear){
    F_at_age[t, ] <- F_yr[t] * sel
    Z_at_age[t, ] <- F_at_age[t, ] + M
  }
  
  # 2. Beverton–Holt  --------------------------------------------------
  Rec <- numeric(nyear)
  Rec[1] <- alpha * 1E5 / (1+(beta*1E5)^c)  # 
  dev_logR=arima.sim(list(order=c(1,0,0),ar=0.1),n=nyear)*std_logR
  
  # 3. cohort ---------------------------------------------
  surv_at_age <- matrix(NA, nrow=nyear, ncol=nage)
  NA_imm   <- matrix(0, nrow=nyear, ncol=nage)
  NA_mat   <- matrix(0, nrow=nyear, ncol=nage)
  N_at_age <- matrix(0, nrow=nyear, ncol=nage)
  N_at_len <- matrix(0,    nrow=nyear, ncol=nlen)
  B_at_len <- matrix(0,    nrow=nyear, ncol=nlen)
  B_at_age <- matrix(0,    nrow=nyear, ncol=nage)
  SB_at_len<- matrix(0,    nrow=nyear, ncol=nlen)
  SB_at_age<- matrix(0,    nrow=nyear, ncol=nage)
  SSB      <- numeric(nyear)
  
  
  NA_imm[1,1] <- Rec[1]
  NA_mat[1,1] <- 1e-20
  for(a in 2:nage){
    surv_at_age[1, a] <- NA_imm[1,a-1] * exp(-Z_at_age[1,a-1])
    
    NA_imm[1,a] <- surv_at_age[1, a] * (1 - mat[a])
    NA_mat[1,a] <- surv_at_age[1, a] * mat[a]
  }
  N_at_age[1,]   <- NA_imm[1,] + NA_mat[1,]
  N_at_len[1,]   <- N_at_age[1,] %*% pla
  B_at_age[1,]   <- N_at_age[1,] * W_at_age
  B_at_len[1,]   <- N_at_len[1,] * W_at_len
  SB_at_len[1,]  <- B_at_len[1,] * mat_L
  SB_at_age[1,]  <- B_at_age[1,] * mat
  SSB[1]         <- sum(SB_at_age[1,])
  
  for(t in 2:nyear){
    
    logR=log(alpha * SSB[t-1] / (1+(beta*SSB[t-1])^c))
    Rec[t] <- exp(logR+dev_logR[t])
    
    
    NA_imm[t,1] <- Rec[t]
    NA_mat[t,1] <- 1e-20
    
    
    for(a in 2:nage){
      surv_at_age[t, a] <- NA_imm[t-1,a-1] * exp(-Z_at_age[t-1,a-1])
      NA_imm[t,a] <- surv_at_age[t, a] * (1 - mat[a])
      NA_mat[t,a] <- surv_at_age[t, a] * mat[a]
    }
    
    N_at_age[t,]   <- NA_imm[t,] + NA_mat[t,]
    B_at_age[t,]   <- N_at_age[t,] * W_at_age
    #N_at_len[t,]   <- N_at_age[t,] %*% pla
    #B_at_len[t,]   <- N_at_len[t,] * W_at_len
    SB_at_age[t,]  <- B_at_age[t,] * mat
    SSB[t]         <- sum(SB_at_age[t,])
  }
  
  
  CN_at_age <- matrix(0, nrow=nyear, ncol=nage)
  CN_at_len <- matrix(0, nrow=nyear, ncol=nlen)
  CB_at_len <- matrix(0, nrow=nyear, ncol=nlen)
  CB_at_age <- matrix(0, nrow=nyear, ncol=nage)
  TN        <- numeric(nyear)
  TB        <- numeric(nyear)
  CN        <- numeric(nyear)
  CB        <- numeric(nyear)
  #weight_at_age <- as.numeric( pla %*% W_at_len )
  weight_at_age=W_at_age
  for(t in 1:nyear){
    for(a in 1:nage){
      ca_imm <- NA_imm[t,a]*(1-exp(-Z_at_age[t,a]))*(F_at_age[t,a]/Z_at_age[t,a])
      ca_mat <- NA_mat[t,a]*(1-exp(-F_at_age[t,a]))
      CN_at_age[t,a] <- ca_imm + ca_mat
      CB_at_age[t,a] = CN_at_age[t,a] * weight_at_age[a]
    }
    CN_at_len[t,] <- CN_at_age[t,] %*% pla
    CB_at_len[t,] <- CN_at_len[t,] * W_at_len
    TN[t] <- sum(N_at_age[t,])
    TB[t] <- sum(B_at_age[t,])
    CN[t] <- sum(CN_at_len[t,])
    CB[t] <- sum(CB_at_age[t,])
  }
  
  
  set.seed(iter)
  surv_err   <- rnorm(nyear, 0, std_SN)
  RVN_at_age <- N_at_age * sel * exp(surv_err)*1e-7
  RVB_at_age <- RVN_at_age * W_at_age
  #RVN        <- pmax(rowSums(RVN_at_age), 1e-6)
  RVN=TN*1e-7 *exp(surv_err)
  RVB        <- rowSums(RVB_at_age)
  
  sim.data=list(
    
    q_surv = q_surv,
    len_mid = len_mid,
    nyear=nyear,
    nage=nage,
    nlen=nlen,
    ages=c(1:nage),
    weight=matrix(rep(W_at_len,nyear),nrow=nlen,ncol=nyear),
    mat=matrix(rep(mat,nyear),nrow=nage,ncol=nyear),
    
    sel = sel,
    SSB_at_age = SB_at_age[1151:1200,],
    B_at_age = B_at_age[1151:1200,],
    F_at_age = F_at_age[1151:1200,],
    # M_at_age = M_at_age[(1*12):(100*12),],
    #N_at_len = N_at_len[(81*12):(100*12),],
    N_at_age = N_at_age[1151:1200,],
    # B_at_len = B_at_len[(81*12):(100*12),],
    #SB_at_len = SB_at_len[(81*12):(100*12),],
    #CN_at_len = CN_at_len[(81*12):(100*12),],
    CN_at_age = CN_at_age[1151:1200,],
    CB_at_len = CB_at_len[1151:1200,],
    CB_at_age = CB_at_age[1151:1200,],
    
    Rec = NA_imm[1151:1200,1],
    F_yr=F_yr[1151:1200],
    TN = TN[1151:1200],
    TB = TB[1151:1200],
    CN = CN[1151:1200],
    CB = CB[1151:1200],
    SSB = SSB[1151:1200],
    RVB=RVB[1151:1200],
    RVN=RVN[1151:1200],
    NA_imm    = NA_imm[1151:1200, ],
    NA_mat    = NA_mat[1151:1200, ]
  )
  
  save(sim.data, file=file.path(iter_folder, paste0("sim_rep", iter, ".RData")))
  
  
