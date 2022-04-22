## Hao Sun, sunhao92@jlu.edu.cn

input_geno_h2_seed_return_QTL_PHNO=function(x0,h2,seed){
  
  total_marker_number<-dim(x0)[1]
  sample_size<-dim(x0)[2]
  
  num_QTL=19

set.seed(seed)
marker_a<-sample(1:total_marker_number,num_QTL)  ##additive

set.seed(seed+100)
marker_d<-sample(1:total_marker_number,num_QTL)  ##dominant

effect=1

## A-variance
  AQTL_matrix=(x0[marker_a,])
  AVQTL=apply(AQTL_matrix,1,var)
  Avar_genotype=sum(AVQTL)

## D-variance
  DQTL_matrix=(x0[marker_d,])
  DQTL_matrix[DQTL_matrix==2]=0
  DVQTL=apply(DQTL_matrix,1,var)
  Dvar_genotype=sum(DVQTL)

## g-variance
VAD=Avar_genotype+Dvar_genotype

## e-VARIANCE	
var_erro=VAD*(1-h2)/h2

## e
library(MASS)
mu=rep(0,sample_size)
sigma=diag(var_erro,sample_size,sample_size)
set.seed(seed+200)
e=mvrnorm(1,mu,sigma)

## simulated y
Aeffect_QTL=apply(AQTL_matrix,2,sum)
Deffect_QTL=apply(DQTL_matrix,2,sum)
pre_y <-  as.numeric(Aeffect_QTL) + as.numeric(Deffect_QTL) + as.numeric(e) 

return(list(pre_y,
            marker_a,
			marker_d,
			Avar_genotype,
			Dvar_genotype,
			var_erro))
}



FST<-function(Geno,Y){   
             x=as.matrix(Geno)
             y=as.matrix(Y)

n0_num=as.numeric(table(y)[1])  # 0 
n1_num=as.numeric(table(y)[2])  # 1
N=nrow(y)

nc=N-n0_num^2/N-n1_num^2/N
n_mean=N/2

id_n0=which(y==0)	
id_n1=which(y==1)
f0=as.matrix(apply(x[,id_n0],1,sum)/( 2*n0_num))   ##等位基因频率
f1=as.matrix(apply(x[,id_n1],1,sum)/( 2*n1_num))

f_mean=(f0*n0_num+f1*n1_num)/N	

s2=(n0_num*(f0-f_mean)^2+n1_num*(f1-f_mean)^2)/n_mean

h0= (apply(x[,id_n0],1,function(x){length(which(x==1))}))/n0_num      ##观测杂合度
h1= (apply(x[,id_n1],1,function(x){length(which(x==1))}))/n1_num
h_mean=(h0*n0_num+h1*n1_num)/N

c=0.5*h_mean

a=n_mean/nc*(s2-(f_mean*(1-f_mean)-0.5*s2-0.25*h_mean)/(n_mean-1))

b=n_mean/(n_mean-1)*(f_mean*(1-f_mean)-0.5*s2-h_mean*(2*n_mean-1)/(4*n_mean))

				
fst<-a/(a+b+c)

return(fst)
}







for(seeduse in 1:10){
gen<-read.csv("sim.gen.txt",header=TRUE)

gen[gen==1]=2
gen[gen==0]=1
gen[gen==-1]=0

x0=gen
dim(x0)

sim=input_geno_h2_seed_return_QTL_PHNO(x0,0.5,seeduse)
pre_y=sim[[1]]

# sim[[4]]
# sim[[5]]
# sim[[6]]

marker_a=sim[[2]]
marker_d=sim[[3]]

## population 0 and 1
pre_y1=pre_y
mu=mean(pre_y1)

for(ii in 1:length(pre_y1)){
if(pre_y1[ii]<=mu){
                    pre_y1[ii]=0
					}else{pre_y1[ii]=1}
}

#####################################################################
## LMM

## gcta VARIANCE
system(paste("vcftools --vcf sim.gen.vcf --plink --out temp"))
system(paste("plink --file temp --make-bed --out tempbed")) 
system(paste("gcta64 --bfile tempbed  --make-grm --out A "))  #tmp_",phenotype_name,".txt --pheno tmp_",phenotype_name,".pheno --out tmp_",phenotype_name,sep=""))
system(paste("gcta64 --bfile tempbed  --make-grm-d --out D "))

system(paste(" rm temp* "))

kin=data.frame(rbind("A","D.d"))
write.table(as.matrix(kin),"2grm",quote=FALSE,row.names=FALSE,col.names=FALSE)

phe.txt=data.frame(c(1:1000),c(1:1000),pre_y1)

write.table(phe.txt,"phe.txt",row.names=FALSE,col.names=FALSE)
system(paste("gcta64 --reml --reml-alg 2 --reml-no-constrain  --reml-maxit 2000 --mgrm  2grm --pheno phe.txt --out add_dom" ))
variance=read.table("add_dom.hsq",fill=TRUE)
variance


system(paste(" rm phe* "))
system(paste(" rm A* "))
system(paste(" rm D* "))
system(paste(" rm add* "))

vadd=as.numeric(variance[2,2])
vdom=as.numeric(variance[3,2])
verr=as.numeric(variance[4,2])

##############################
## GET gcta grm
sample_size<-dim(x0)[2]
marker_number<-dim(x0)[1]
A_matrix=x0

f<-as.matrix(apply(A_matrix,1,sum))/(sample_size*2)
A_matrix_scale<-(A_matrix-2*f)/sqrt(2*f*(1-f))
dim(A_matrix_scale)
A_matrix_scale<-as.matrix(A_matrix_scale)
A_cor=t(A_matrix_scale)%*%(A_matrix_scale)/marker_number  ###Additive grm

d_matrix<-matrix(0,marker_number,sample_size)
for(i in 1:marker_number){
  x<-as.matrix(A_matrix[i,])
  x[x==1]=2*f[i,1]
  x[x==2]=4*f[i,1]-2
  d_matrix[i,]<-x
}
d_matrix_scale<-(data.frame(d_matrix)-2*f*f)/(2*f*(1-f))
d_matrix_scale<-as.matrix(d_matrix_scale)
D_cor=t(d_matrix_scale)%*%(d_matrix_scale)/marker_number  ###Dominant grm

##############################
## get p-value ADDO 

V=A_cor*vadd+D_cor*vdom+diag(sample_size)*verr
ni = solve(eigen(V)$vectors %*% diag(x=sqrt(eigen(V)$values)) %*% solve(eigen(V)$vectors)) 
		
yt=ni%*%pre_y1
u=as.matrix(rep(1,sample_size))


result_matrix1=matrix(0,marker_number,6)
for(i in 1:marker_number){
  xa<-as.matrix(A_matrix_scale[i,])
  xd<-as.matrix(d_matrix_scale[i,])

  xuad<-ni %*% cbind(u,xa,xd)

result=summary(lm(yt~-1+xuad))$coefficients

beta_a=result[2,1]
sd_a=result[2,2]
p_a=result[2,4]

beta_d=result[3,1]
sd_d=result[3,2]
p_d=result[3,4] 

if(i%%100 ==0){print(i)}

result_matrix1[i,]=c(beta_a,sd_a,p_a,beta_d,sd_d,p_d)

}

dim(result_matrix1)
pa_lmm=result_matrix1[,3]
pd_lmm=result_matrix1[,6]


############################
## fst

sta=FST(x0,pre_y1)


power_test<-function(have,p_sta){   ##have effect， statistic

SNP_number=length(p_sta) ##total num
all_SNP_rank<-rank(-p_sta)/SNP_number ##rank


right_0.001=0
right_0.005=0
right_0.01=0
right_0.05=0


for(i in 1:length(have)){
     aim<-have[i]
     score<-all_SNP_rank[aim]    #the rank of the QTL 
	 
	 if(score <=0.001){right_0.001=right_0.001+1}
     if(score <=0.005){right_0.005=right_0.005+1}
	 if(score <=0.01){right_0.01=right_0.01+1}
	 if(score <=0.05){right_0.05=right_0.05+1}
	 }
	 
right<-c(right_0.001,right_0.005,right_0.01,right_0.05)

return(right)			  
}
   

fsta=power_test(marker_a,sta)
fstd=power_test(marker_d,sta)
lmma=power_test(marker_a,-pa_lmm)
lmmd=power_test(marker_d,-pd_lmm)

# fsta
# fstd
# lmma
# lmmd

out=rbind(
fsta,
fstd,
lmma,
lmmd)

write.table(cbind(pre_y1,pre_y),paste("result_pre_y_",seeduse,".txt",sep=""),col.names=FALSE)
write.table(cbind(marker_a,marker_d),paste("result_marker_ad_",seeduse,".txt",sep=""),row.names=FALSE,col.names=FALSE)
write.table(variance,paste("result_variance_ad_",seeduse,".txt",sep=""),row.names=FALSE,col.names=FALSE)
write.table((out),paste("result_power_",seeduse,".txt",sep=""),col.names=FALSE)
}

