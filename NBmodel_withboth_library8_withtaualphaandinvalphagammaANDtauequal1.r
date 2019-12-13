#August 26,2015 Gainesville FL
#it considers a Negative Binomial model 
#See model_extension_withboth6.pdf section for model specification<- i stii need to update this pdf document
#This model/program was started August 2015 while visiting Lauren
#The correction bias for counts in Line1,q1, can be different from
#the correction bias for counts in Line2
#The model has the following structure 
#x, y and z are negative binomial with means
#The mean structure is
#                         Tester                      Line 1                       both
#Line 1(mated)            q11*beta_i(1/alpha)         q12*betai*alpha              tau[(1-q11)alpha+(1-q12)alpha]betai
#Line 2(virgin)           q21*beta_i(1/delta)gamma    q22*betai*delta*gamma        tau[(1-q21)/delta+(1-q12)delta]betai*gamma
#It is NOT required.
#q11 is the prob of assigning given read comes from tester
#q22 is prob of assignemnt given read comes from line
#source("C:/Users/Luis Leon-Novelo/Dropbox/projects/RNAextended/NBmodel_WITHBOTH5_gammaalphadeltainboth.R")

#To induce the negative binomial distribution 
#we introduce a set latent rvs nu
#See NBmodel_testing.R where data the main function of this library is called
#and the functions in this library are tested

# CALCULATION OF THE MLES for gamma  GIVEN THE SAMPLE
gam.mles.data = function(x){
  n = length(x)
  xb = mean(x)
  xd = mean(log(x))
  s = log(xb)-xd
  a0 = (3.0-s+sqrt((s-3.0)^2+24.0*s))/12.0/s
  l = 1
  repeat{
    ans = (log(a0)-digamma(a0)-s)
    a1 = a0-ans/(1.0/a0-trigamma(a0))
    if(abs(ans) <= 1.0e-7 | l >= 30){break}
    a0 = a1
    l = l+1}
  ah = a1; bh = xb/a1
  return(c(ah,bh))
}
  

ralpha=function(alpha,x1,y1,z1,beta1,nu,tau,q1,sigma_alpha,sd_MH=0.05){
  #log transform alpha to improve the MH mixing
  

  out=list(alpha=alpha,accept=0)
  logalpha=log(alpha)
  candidate=rnorm(1,logalpha,sd=sd_MH)
    expcan=exp(candidate)   
    logaccept=
      (1/alpha-1/expcan)*(q1[1]*sum(beta1*nu$x1)+ (1-q1[1])*sum(beta1*nu$z1)*tau)+
      (alpha-expcan)*    (q1[2]*sum(beta1*nu$y1)+ (1-q1[2])*sum(beta1*nu$z1)*tau)+
      (sum(y1)-sum(x1))*(candidate-logalpha)+
      sum(z1)  *log(((1-q1[1])/expcan+(1-q1[2])*expcan    )/
                    ((1-q1[1])/alpha+(1-q1[2])*alpha      ))+
      #likelihood  
      #-(a_alpha+1)*log(candidate/alpha)-b_alpha*(1/candidate-1/alpha)#inverse gamma prior for alpha
      #(alpha-expcan)*b_alpha+(a_alpha)*(candidate-logalpha)  #gamma prior alpha
      (logalpha^2-candidate^2)/(2*sigma_alpha^2)#lognormal prior for alpha
    #print(paste("logaccept",logaccept))
    if(log(runif(1))<logaccept) out=list(alpha=expcan,accept=1);
  
  return(out)
}



rdelta=function(delta,x2,y2,z2,beta2,gamma_,nu,tau,q2,sigma_delta,sd_MH=0.05){
  #log transform alpha to improve the MH mixing
  out=list(delta=delta,accept=0)
  logdelta=log(delta)
  candidate=rnorm(1,logdelta,sd=sd_MH)
  expcan=exp(candidate)
  logaccept=
    (1/delta-1/expcan)*gamma_*(q2[1]*sum(beta2*nu$x2)+ (1-q2[1])*sum(beta2*nu$z2)*tau)+
    (delta-expcan)*gamma_*    (q2[2]*sum(beta2*nu$y2)+ (1-q2[2])*sum(beta2*nu$z2)*tau)+
    (sum(y2)-sum(x2))*(candidate-logdelta)+
    sum(z2)  *log(((1-q2[1])/expcan+(1-q2[2])*expcan      )/((1-q2[1])/delta+(1-q2[2])*delta      ))+   #likelihood
      #-(a_delta+1)*log(candidate/delta)-b_delta*(1/candidate-1/delta)#inverse gamma prior for delta
      #(delta-expcan)*b_delta+a_delta*(candidate-logdelta)  #gamma prior alpha and Jacobian
    (logdelta^2-candidate^2)/(2*sigma_delta^2)#lognormal prior for alpha
  
    #print(paste("logaccept",logaccept))
    
    if(log(runif(1))<logaccept) out=list(delta=expcan,accept=1);
     #print(paste("delta=",out$delta))
  return(out)
}




rphi_MH=function(phi,nu,I1,I2,a_phi,b_phi,sd_MH=0.05){
	out=list(phi=phi,accept=0)
	candidate=rnorm(1,phi,sd=sd_MH)
	if(candidate>0){
		
		logaccept=	
		(1/candidate-1/phi)*sum(log(unlist(nu)))-
		(1/candidate-1/phi)*sum(   (unlist(nu)))+
		3*(I1+I2)*(	(1/candidate)*log(1/candidate)-(1/phi)*log(1/phi)-
					lgamma(1/candidate)+lgamma(1/phi))+
		#(a_phi-1)*log(candidate/phi)-b_phi*(candidate-phi)	#gamma prior for phi
		-(a_phi+1)*log(candidate/phi)-b_phi*(1/candidate-1/phi)#Inverse gamma prior
	if(log(runif(1))<logaccept) out=list(phi=candidate,accept=1);
	}
	return(out)
	}


rphi_MH_new=function(phi,nu,I1,I2,a_phi,b_phi,sd_MH){
  #Tried to use the transformation u=-log phi
  #to improve mixong but the mixing actually gets worse
  #print(phi)
  
  out=list(phi=phi,accept=0)
  u=log(1/phi)#This is the transformation we are using
  candidate=rnorm(1,u,sd=sd_MH)
  expu=1/phi
  expcan=exp(candidate)
    
  logaccept=	
      (candidate*expcan-u*expu)*                  3*(I1+I2)+
      (lgamma(expu)-lgamma(expcan))* 3*(I1+I2)+
       (expu-expcan)*                 sum(   (unlist(nu)))+
       (expcan-expu)*                 sum(log(unlist(nu)))+
      #(a_phi-1)*log(candidate/phi)-b_phi*(candidate-phi)	#gamma prior for phi
       (candidate-u)*(a_phi+1)+b_phi*(expu-expcan)+#Inverse gamma prior
        (u-candidate)#Jacobian
    
    if(log(runif(1))<logaccept) out=list(phi=1/expcan,accept=1);
  
  return(out)
}


#############################################################
#############################################################
gibbs_NB_niters_withboth=function(
n_iter=10,x1,y1,z1,x2,y2,z2,
all,q1,q2=q1,a,b,sigma,I1,I2,alpha_sd_MH=0.5,delta_sd_MH,phi_sd_MH=0.5){
acceptaceratealpha=acceptaceratedelta=acceptaceratephi=0;
y..=sum(c(y1,y2));z..=sum(c(z1,z2));

for(ii in 1:n_iter){	
all$b_beta=rgamma(1,a$b_beta+(I1+I2)*a$beta,rate=b$b_beta+sum(all$beta1+all$beta2))
  
  
  
#  (all$gamma<-rgamma(1,a$gamma+sum(x2)+sum(y2)+sum(z2),
#rate=with(all,b$gamma+
#				q2[1]*		                            sum(beta2*nu$x2)/delta+
#				q2[2]*                                sum(beta2*nu$y2)*delta+
#        ((1-q2[1])/delta+(1-q2[2])*      delta)*    sum(beta2*nu$z2)*tau )))

###updating betas and nus
for(i in 1:I1){
all$beta1[i]<-rgamma(1,a$beta+x1[i]+y1[i]+z1[i],
                     rate=with(all,b_beta+
                                 q1[1]*nu$x1[i]/alpha+
                                 q1[2]*nu$y1[i]*alpha+
                                 ((1-q1[1])/alpha+(1-q1[2])*alpha)      *nu$z1[i]*tau   ))

all$nu$x1[i]=rgamma(1,1/all$phi+x1[i],rate=with(all,1/phi+q1[1]*(1/alpha)*                  beta1[i]))
all$nu$y1[i]=rgamma(1,1/all$phi+y1[i],rate=with(all,1/phi+q1[2]*alpha*	                    beta1[i]))
all$nu$z1[i]=rgamma(1,1/all$phi+z1[i],rate=with(all,1/phi+((1-q1[1])/alpha+(1-q1[2])*alpha)*beta1[i]*tau))
}

for(i in 1:I2){
all$beta2[i]<-rgamma(1,a$beta+x2[i]+y2[i]+z2[i],
                     rate=with(all,b_beta+gamma*(
                                 q2[1]*nu$x2[i]/delta+
                                 q2[2]*nu$y2[i]*delta+
                                 ((1-q2[1])/delta+(1-q2[2])*delta)*nu$z2[i]*tau )))

(all$nu$x2[i]<-rgamma(1,1/all$phi+x2[i],rate=with(all,1/phi+q2[1]*    (1/delta)*beta2[i]*        gamma)))
(all$nu$y2[i]<-rgamma(1,1/all$phi+y2[i],rate=with(all,1/phi+q2[2]*        delta*beta2[i]* 	     gamma)))
(all$nu$z2[i]<-rgamma(1,1/all$phi+z2[i],rate=with(all,1/phi+((1-q2[1])/delta+(1-q2[2])*delta)*beta2[i]*gamma*tau )))
}



tem=with(all,rdelta(delta=delta,x2=x2,y2=y2,z2=z2,beta2=beta2,
                    gamma_=gamma,nu=nu,tau=tau,q2=q2,sigma_delta=sigma$delta,sd_MH=delta_sd_MH))
all$delta=tem$delta
acceptaceratedelta=acceptaceratedelta+tem$accept

(tem=with(all,ralpha(alpha=alpha,x1=x1,y1=y1,z1=z1,beta1=beta1,
                     nu=nu,tau=tau,q1=q1,sigma_alpha=sigma$alpha,sd_MH=alpha_sd_MH)))
all$alpha=tem$alpha
acceptaceratealpha=acceptaceratealpha+tem$accept

#all$tau=rgamma(1,a$tau+z..,rate=b$tau+with(all,
#                                      ((1-q1[1])/alpha+(1-q1[2])*alpha)*sum(beta1*nu$z1)+     
#                                      ((1-q2[1])/delta+(1-q2[2])*delta)*sum(beta2*nu$z2)*gamma
#                                      ))

tem=rphi_MH(all$phi,all$nu,I1,I2,a$phi,b$phi,sd_MH=phi_sd_MH)
all$phi=tem$phi
acceptaceratephi=acceptaceratephi+tem$accept
}
return(list(all=all,acceptaceratephi=acceptaceratephi/n_iter,acceptaceratealpha=acceptaceratealpha/n_iter,
            acceptaceratedelta=acceptaceratedelta/n_iter));
}

##############################################
#
#xs=c(1,4,0,2,0,3)
#ys=c(0,0,0,3,5,2)
#zs=c(0,12,34,21,32,12)

## fusion_id        qsim_both        qsim_line      qsim_tester       M_num_reps 
#      "r101"      "F10001_SI
#xs=c( 2,  	7, 10, 		15,  	6, 53)
#ys=c(5, 	12, 3, 		16, 	36, 55)
#zs=c(118, 	159,205, 	278, 	508, 640)
#q_=c(0.08333333, 0.08333333, 0.83333333)
#############################################
gibbs_NegBin_withboth=function(nsim=100,nburnin=100,nlag=10,
xs,		#Vector of reads assigned to "tester" allele
ys,		#Vector of reads assigned to "line" allele
zs,		#Vector of reads assigned to both alleles
is=rep(1:3,2),		#index for biorep
js=c(rep(1,3),rep(2,3)),	#index for Line1=1, vigin=Line2
sigma_alpha=sigma_alpha,
a_beta=1,
#b_beta=0.01, 
a_b_beta=0.01,
b_b_beta=0.01,
a_gamma=0.01,
b_gamma=0.01, 
sigma_delta=1,   #standard deviation for the normal prior distribution for log delta
a_tau=1,
b_tau=1,
a_phi=2.01,
b_phi=0.05,				#phi is apriori small	
alpha_sd_MH=0.1,   #Standard deviations for the nurmal proposal to sample via MH
delta_sd_MH=0.1,
phi_sd_MH=0.1,			#Initial SD for the normal proposal to sample from phi via MH
q1=c(0.2,0.2,0.6),		#Bias correction hyperparameter in Line 1
q2=q1,					#Bias correction hyperparameter in Line 2
plots=FALSE					#If true plots are depicted
){
I1=length(unique(is[which(js==1)]))
I2=length(unique(is[which(js==2)]))

#tem=xs
#xs=ys
#ys=tem
#Initiating parameters		
x1=xs[which(js==1)]; 
y1=ys[which(js==1)]; 
z1=zs[which(js==1)]; 
x2=xs[which(js==2)]; 
y2=ys[which(js==2)]; 
z2=zs[which(js==2)]; 

#Checking that I have not a vector of all zeros (this messes with the convergence, 

if(sum(x1+y1+z1)==0) z1[1]=q1[3];
if(sum(x2+y2+z2)==0) z2[1]=q2[3];

if(sum(z1)==0)	z1[which.max(x1+y1)[1]]= q1[3]
if(sum(z2)==0)	z2[which.max(x2+y2)[1]]= q2[3]
if(sum(x1)==0)	x1[which.max(z1)[1]]=    q1[1]
if(sum(x2)==0)	x2[which.max(z2)[1]]=    q2[1]
if(sum(y1)==0)	y1[which.max(z1)[1]]=    q1[2]
if(sum(y2)==0)	y2[which.max(z2)[1]]=    q2[2]
#x1[which(x1==0)]=q1[1]
#y1[which(y1==0)]=q1[2]
#z1[which(z1==0)]=q1[3]*2

#x2[which(x2==0)]=q2[1]
#y2[which(y2==0)]=q2[2]
#z2[which(z2==0)]=q2[3]*2


xs=c(x1,x2);ys=c(y1,y2);zs=c(z1,z2)
y..=sum(ys);z..=sum(zs);
########
#print("data after correcting for zero counts")
#print("counts")
#print(paste(paste(x1,collapse=","),"      ",paste(y1,collapse=","),"     ",paste(z1,collapse=",")))
#print(paste(paste(x2,collapse=","),"      ",paste(y2,collapse=","),"     ",paste(z2,collapse=",")))
#print("q1 correction bias")
#print(q1)
#Checking model assumption
print("For line 1 z., x. q1/(1-q1)+y.q2/(1-q2) and the log of ratio")
tem1=sum(z1)
tem2=sum(x1)*(1-q1[1])/q1[1]+sum(y1)*(1-q1[2])/q1[2]
print(c(tem1,tem2,log(tem1/tem2)         ))
tau=tem1/tem2
print("For line 2")
tem1=sum(z2)
tem2=sum(x2)*(1-q2[1])/q2[1]+sum(y2)*(1-q2[2])/q2[2]
print(c(tem1,tem2,log(tem1/tem2)         ))
tau=mean(c(tau,tem1/tem2))



print(paste("point estimate of tau",tau))
print(paste("but we are making tau equal ",tau<-1))
###Initiating parameter values
delta_est=1
alpha_est=1 

ftem=function(xx){max(c(xx,0.1))}
beta1=sapply(x1+y1+z1,ftem)/            sum(c(q1[1],q1[2]*alpha_est,((1-q1[1])+(1-q1[2])*alpha_est)*tau))
gammatimesbeta2=sapply((x2+y2+z2),ftem)/sum(c(q2[1],q2[2]*delta_est,((1-q2[1])+(1-q2[2])*delta_est)*tau))

gamma_est=1#mean(gammatimesbeta2)/mean(beta1)
beta2=gammatimesbeta2/gamma_est
#mean(beta2)/mean(beta1)


alpha_est=  sqrt((q1[1]*sum(y1))/(q1[2]*max(sum(x1),.5)))

delta_est=  sqrt((q2[1]*sum(y2))/(q2[2]*max(sum(x2),.5)))




#print(paste("alpha_est=",alpha_est))
#print(paste("alpha_est=",gamma_est))
#print(paste("delta_est=",delta_est))
#print("betas")
#print(beta1)
#print(beta2)

tem=gam.mles.data(c(beta1,beta2))
tem[1]=min(max(tem[1],10^(-3)),10^5)  
tem[2]=min(max(tem[2],10^(-3)),10^5)  
a_beta=tem[1]
a_b_beta=2*tem[2]^(-1)#tem[2]^(-3/2)
b_b_beta=2#tem[2]^(-1/2)
#print(paste("mleS FOR BETAI",round(tem[1],3),round(tem[2],3)))
#a_b_beta/b_b_beta^2



#print(paste("setting a_beta equal to the MLE when assuming alpha=delta=gamma=1: 
#            a_beta=",a_beta,
#            "\n The MLE for b_b_beta is set to 1/MLE b_beta=",b_b_beta))

a=  list(beta=a_beta,b_beta=a_b_beta,gamma=a_gamma,tau=a_tau,phi=a_phi)
b=	list(b_beta=b_b_beta,gamma=b_gamma,tau=b_tau,phi=b_phi)
sigma=list(alpha=sigma_alpha,delta=sigma_delta)


nu=list(x1=rep(1,I1),y1=rep(1,I1),z1=rep(1,I1),
		x2=rep(1,I2),y2=rep(1,I2),z2=rep(1,I2))


		#mean((x2./I2)/max(.5,(x1./I1)),(y2./I2)/max(.1,(y1./I1)),(z1./I1)/
		#max(.1,(z2./I2))        )
		#mean(x2./(q_[1]*I2),y2./(alpha*q_[2]*I2),z1./(((alpha+1)/2)*q_[3]*I1))/mu

#print(paste("alphax=",alphax," alphay=",alphay," gamma=",gamma_))
#print("beta1");print(beta1);
#print("beta2");print(beta2);

#print("Expected counts assuming initial values")
#print(
#rbind(	c(q_[1]*sum(beta1)*alphax,		 q_[2]*sum(beta1)*		alphay,			q_[3]*sum(beta1)),
#		c(q_[1]*sum(beta2)*gamma_*alphax,q_[2]*sum(beta2)*gamma_*alphay*delta,	q_[3]*sum(beta2)))
#)

#Printing counts
totals=rbind(c(sum(x1),sum(y1)),c(sum(x2),sum(y2)))
rownames(totals)=c("row1","row2")
colnames(totals)=c("x","y")
#totals/q_[c(1,2)]
phi=b_phi/(a_phi+1)#Mode of an inverse gamma



all=list(beta1=beta1,beta2=beta2,alpha=alpha_est,gamma=gamma_est,delta=delta_est,nu=nu,b_beta=(I1+I2)/sum(beta1+beta2),tau=tau,phi=phi)


#print("Sample Proportion")
sample.props<-addmargins(prop.table(totals))
totalswboth=cbind(totals,c(sum(z1),sum(z2)))
#colnames(totalswboth)=c("x","y","z")
#print("Total counts")
#print(totalswboth)
#print("Total counts divided by q")
#print(cbind(totalswboth[,1]/q_[1],totalswboth[,2]/q_[2],totalswboth[,3]/q_[3]))


####Adjusting the variance of the normal proposal for  the Metropolis Hasting to sample both alpha and phi
#algorith to sample alpha in the Gibbs sampler
acceptaceratealpha=acceptaceratedelta=	acceptaceratephi=	1;
maxiter=100
counter=0
while(
    (acceptaceratealpha<0.25 | acceptaceratealpha>0.75 |
     acceptaceratedelta<0.25 | acceptaceratedelta>0.75 |   
		 acceptaceratephi<0.25 	| acceptaceratephi	>0.75) & counter<maxiter){
  
#alpha_sd_MH=delta_sd_MH=phi_sd_MH=0.1
  
tem=gibbs_NB_niters_withboth(n_iter=100,x1,y1,z1,x2,y2,z2,all,q1,q2,a,b,sigma,I1,I2,
                             alpha_sd_MH=alpha_sd_MH,delta_sd_MH=delta_sd_MH,phi_sd_MH=phi_sd_MH)
all=tem$all
acceptaceratephi=		  tem$acceptaceratephi
acceptaceratealpha=		tem$acceptaceratealpha
acceptaceratedelta=  	tem$acceptaceratedelta
if(acceptaceratephi<0.25){phi_sd_MH=phi_sd_MH/2}
if(acceptaceratephi>0.75){phi_sd_MH=2*phi_sd_MH}

if(acceptaceratealpha<0.25){alpha_sd_MH=alpha_sd_MH/2}
if(acceptaceratealpha>0.75){alpha_sd_MH=2*alpha_sd_MH}

if(acceptaceratedelta<0.25){delta_sd_MH=delta_sd_MH/2}
if(acceptaceratedelta>0.75){delta_sd_MH=2*delta_sd_MH}


counter=counter+1

if(plots==TRUE){
  library("MCMCpack")
print(paste("Iteration for MH alpha and phi", counter))
print(paste("acceptace rate alpha", acceptaceratealpha))
print(paste("next SD for MH for alpha ", alpha_sd_MH))
print(paste("acceptace rate delta", acceptaceratedelta))
print(paste("next SD for MH for delta ", delta_sd_MH))


print(paste("acceptace rate phi", acceptaceratephi))
print(paste("next SD for MH for phi ", phi_sd_MH))

print("alpha delta phi")
with(tem,print(c(acceptaceratealpha,acceptaceratedelta,acceptaceratephi)))
}

}

tem=gibbs_NB_niters_withboth(n_iter=nburnin,x1,y1,z1,x2,y2,z2,all,q1,q2,a,b,sigma,I1,I2,
                             alpha_sd_MH=alpha_sd_MH,delta_sd_MH=delta_sd_MH,phi_sd_MH=phi_sd_MH)
all=tem$all

if(plots==TRUE){
print("Acceptance rates in burn in alpha delta phi")
with(tem,print(c(acceptaceratealpha,acceptaceratedelta,acceptaceratephi)))
}


tau_chain=b_beta_chain=phi_chain=nux11_chain=nuy11_chain=delta_chain=alpha_chain=gamma_chain=rep(NA,nsim)
theta_chain=		matrix(NA,nrow=nsim,ncol=2);colnames(theta_chain)=c("11","21");
beta1_chain=     matrix(NA,nrow=nsim,ncol=I1)  
beta2_chain=     matrix(NA,nrow=nsim,ncol=I2)  


for(ii in 1:nsim){
	all=gibbs_NB_niters_withboth(n_iter=nlag,x1,y1,z1,x2,y2,z2,all,q1,q2,a,b,sigma,I1,I2,
	                             alpha_sd_MH=alpha_sd_MH,delta_sd_MH=delta_sd_MH,phi_sd_MH=phi_sd_MH)$all
	
	delta_chain[ii]=	all$delta	
	beta1_chain[ii,]=	all$beta1
	beta2_chain[ii,]=	all$beta2
	alpha_chain[ii]=	all$alpha
	gamma_chain[ii]=	all$gamma
	nux11_chain[ii]=	all$nu$x1[1]
	nuy11_chain[ii]=	all$nu$y1[1]
  tau_chain[ii]=    all$tau
	phi_chain[ii]=		all$phi
	b_beta_chain[ii]= all$b_beta
  
	theta_chain[ii,"11"]=1/all$alpha   /(1/all$alpha+all$alpha);
	theta_chain[ii,"21"]=(1/all$delta) /(1/all$delta+all$delta)
	}
	
#print(paste(" gamma=",mean.gamma<-mean(gamma_chain)))	
#print(paste(" alphax=",mean.alphax<-mean(alphax_chain)," alphay=",mean.alphay<-mean(alphay_chain)))
#print(paste("interaction =",mean.delta<-mean(delta_chain)))	


#print(paste("posterior means beta1=",beta1.mean," beta2=",beta2.mean))	

#print(paste("correction bias q=",q_))
#print("Expected values")
#expected=
#rbind(beta1.mean*c(q_[1]*mean.alphax		   ,q_[2]*mean.alphay,			q_[3]),
#	  beta2.mean*c(q_[1]*mean.gamma*mean.alphax,q_[2]*mean.gamma*mean.alphay*mean.delta,q_[3]))

#print(expected)
#expected.countswocorrection=
#rbind(beta1.mean*c(mean.alphax,				mean.alphay,					1),
#	  beta2.mean*c(mean.alphax*mean.gamma,	mean.gamma*mean.alphay*mean.delta,1))
	  
##checking convergence	
if(plots==TRUE){
	par(ask=TRUE)
  oldpar<-par()
	plot(mcmc(cbind(gamma_chain)))
	plot(mcmc(cbind(alpha_chain,delta_chain)))
  

  
  
	plot(mcmc(cbind(beta1_chain[,1],beta2_chain[,2],b_beta_chain)),main="beta1[1],beta2[1]")
	plot(mcmc(theta_chain[,c(1,2)]))
	
  plot(mcmc(cbind(nux11_chain,nuy11_chain,phi_chain)))	
	plot(mcmc(tau_chain),main="Tau")
  
  ##Plots for presentation
  par(mfrow=c(2,1),mar=c(3, 5.5, 0.2, 0.2),oma = c(1, 1, 0.2, 0.2))
	plot(alpha_chain,main="",type="l",xlab="",ylab=expression(alpha),cex.lab=2);abline(h=1,col="Blue",lwd=3,xlab="");
	abline(h=mean(alpha_chain),col="brown",lwd=3)
	abline(h=quantile(alpha_chain,c(0.025)),col="brown",lwd=3,lty=2)
	abline(h=quantile(alpha_chain,c(0.975)),col="brown",lwd=3,lty=2)
	
	plot(delta_chain,main="",type="l",ylab=expression(delta),xlab="",cex.lab=2);abline(h=1,col="Blue",lwd=3);
	abline(h=mean(delta_chain),col="brown",lwd=3)
	abline(h=quantile(delta_chain,c(0.025)),col="brown",lwd=3,lty=2)
	abline(h=quantile(delta_chain,c(0.975)),col="brown",lwd=3,lty=2)
	
	par(mfrow=c(3,1))
  plot(beta1_chain[,1],main="",type="l",ylab=expression(beta[1]),xlab="",cex.lab=2)
	abline(h=mean(beta1_chain[,1]),col="brown",lwd=3)
	abline(h=quantile(beta1_chain[,1],c(0.025)),col="brown",lwd=3,lty=2)
	abline(h=quantile(beta1_chain[,1],c(0.975)),col="brown",lwd=3,lty=2)
	
	plot(beta2_chain[,1],main="",type="l",ylab=expression(beta[4]),xlab="",cex.lab=2)
	abline(h=mean(beta2_chain[,1]),col="brown",lwd=3)
	abline(h=quantile(beta2_chain[,1],c(0.025)),col="brown",lwd=3,lty=2)
	abline(h=quantile(beta2_chain[,1],c(0.975)),col="brown",lwd=3,lty=2)
	
	plot(b_beta_chain,main="",type="l",ylab=expression(b[beta]),xlab="",cex.lab=2)
	abline(h=mean(b_beta_chain),col="brown",lwd=3)
  
  
  plot(tau_chain,main="",type="l",ylab=expression(tau),xlab="",cex.lab=2);abline(h=1,col="Blue",lwd=3);
	abline(h=mean(tau_chain),col="brown",lwd=3);abline(h=1,col="Blue",lwd=3);
	abline(h=quantile(tau_chain,c(0.025)),col="brown",lwd=3,lty=2)
	abline(h=quantile(tau_chain,c(0.975)),col="brown",lwd=3,lty=2)
	
	
  plot(nux11_chain,main="",type="l",ylab=expression(nu[x][11]),xlab="",cex.lab=2);
	abline(h=mean(nux11_chain),col="brown",lwd=3)
	
  plot(phi_chain,main="",type="l",ylab=expression(phi),xlab="",cex.lab=2);
	abline(h=mean(phi_chain),col="brown",lwd=3)
  
  par(ask=FALSE)
par(oldpar)
#Computing expected counts with CI
i1=3
print(paste("Observed and Expected counts for biorep ",i1))
tem2=cbind(q1[1]/alpha_chain,q1[2]*alpha_chain,((1-q1[1])/alpha_chain+(1-q1[2])*alpha_chain)*tau_chain)
tem=beta1_chain[,i1]*tem2
temmat=round(
cbind( c(x1[i1],y1[i1],z1[i1]),
      apply(tem,2,FUN=mean),
      apply(tem,2,FUN=function(xx){quantile(xx,c(0.025))}),
      apply(tem,2,FUN=function(xx){quantile(xx,c(0.975))})     
      ),2)
colnames(temmat)=c("obs","mean","q025","q975")
rownames(temmat)=c("xi1","yi1","zi1")
print(temmat)


print("CI for total counts in row 1")
tem=rowSums(beta1_chain)      *tem2

temmat=round(
  cbind( c(sum(x1),sum(y1),sum(z1)),
         apply(tem,2,FUN=mean),
         apply(tem,2,FUN=function(xx){quantile(xx,c(0.025))}),
         apply(tem,2,FUN=function(xx){quantile(xx,c(0.9755))})     
  ),2)
colnames(temmat)=c("obs","mean","q025","q975")
rownames(temmat)=c("xi1","yi1","zi1")
print(temmat)



#Computing expected counts with CI
i2=3
print(paste("Observed and Expected counts for biorep in row 2",i2))
tem2=cbind(q2[1]/delta_chain,q2[2]*delta_chain,((1-q2[1])/delta_chain+(1-q2[2])*delta_chain)*tau_chain)
tem=beta2_chain[,i2]*tem2
temmat=round(
  cbind( c(x1[i2],y1[i2],z1[i2]),
         apply(tem,2,FUN=mean),
         apply(tem,2,FUN=function(xx){quantile(xx,c(0.025))}),
         apply(tem,2,FUN=function(xx){quantile(xx,c(0.9755))})     
  ),2)
colnames(temmat)=c("obs","mean","q025","q975")
rownames(temmat)=c("xi2","yi2","zi2")
print(temmat)


print("CI for total counts in row 2")
tem=rowSums(beta2_chain)      *tem2

temmat=round(
  cbind( c(sum(x2),sum(y2),sum(z2)),
         apply(tem,2,FUN=mean),
         apply(tem,2,FUN=function(xx){quantile(xx,c(0.025))}),
         apply(tem,2,FUN=function(xx){quantile(xx,c(0.9755))})     
  ),2)
colnames(temmat)=c("obs","mean","q025","q975")
rownames(temmat)=c("xi1","yi1","zi1")
print(temmat)

}#end of if(plots==1)


me=   	apply(theta_chain,MARGIN=2,mean)
q025=		apply(theta_chain,MARGIN=2,function(xx){quantile(xx,0.025)})
q975=		apply(theta_chain,MARGIN=2,function(xx){quantile(xx,0.975)})
q50=    apply(theta_chain,MARGIN=2,median)

#print("estimated props")

tem=rbind(me[c("11","12")],me[c("21","22")]);

rownames(tem)=c("row1","row2")
colnames(tem)=c("x","y")
#print(tem);
#print(round(
estimated.props<-addmargins(tem)
#,3))



theta_q025=rbind(q025[c("11","12")],q025[c("21","22")])
theta_q975=rbind(q975[c("11","12")],q975[c("21","22")])
theta_q50=rbind(q50  [c("11","12")],q975[c("21","22")])

summaryt=matrix(NA,nrow=2,ncol=12 )
colnames(summaryt)=
  c("genotype i","x_i.","y_i.","z_i.","x_i./(x_i.+y_i.)","median","mean","q_025","q_975",
    "AIbayesian-pval","flagdifonehalf","binomial-pval")


for(i in 1:2){
  tem2=theta_chain[,i]
  
  summaryt[i,"genotype i"]=i
  summaryt[i,c("x_i.","y_i.","z_i.")]=totalswboth[i,]
  summaryt[i,"x_i./(x_i.+y_i.)"    ]=totals[i,1]/(sum(totals[i,c(1,2)]))  
  
  summaryt[i,"mean"]=                round( mean(tem2),4)
  summaryt[i,c("q_025","median","q_975")]=round(	quantile(tem2,c(0.025,0.5,0.975)),4)
  
  pval=round(	2*min(c(mean(tem2<0.5),mean(tem2>0.5))),4)
  
  
  summaryt[i,"AIbayesian-pval"]=pval
  summaryt[i,"flagdifonehalf"]=ifelse(pval<0.05,1,0)	
  
  
  
  op <- options(warn = (-1)) # suppress warnings
  summaryt[i,"binomial-pval"]=
  ifelse(floor(nn<-sum(totals[i,c(1,2)]))!=nn,NA,binom.test(x=totals[i,1],n=nn,p=0.5)$p.value)
  options(op) # reset the default value, if you want
}


Bayes.pval_independence=2*min(c(mean(alpha_chain/delta_chain<1),mean(alpha_chain/delta_chain>1)))
alphagreater1greaterdelta=mean(alpha_chain>1 & 1>delta_chain)
deltagreater1greateralpha=mean(alpha_chain<1 & 1<delta_chain)

op <- options(warn = (-1)) # suppress warnings
independenceChisqpvalue=chisq.test(x=totals)$p.value
options(op) # reset the default value, if you want

out=list(
tests=summaryt,
counts=totalswboth,
independencetest=list(Chisqpvalue=independenceChisqpvalue,bayesianpvalue=Bayes.pval_independence),
moretests=list(alphagreater1greaterdelta=alphagreater1greaterdelta,deltagreater1greateralpha=deltagreater1greateralpha),
modelparameters_postmean=list(alpha=mean(alpha_chain),gamma=mean(gamma_chain),delta=mean(delta_chain),tau=mean(tau_chain))
)
return(out)
	}
