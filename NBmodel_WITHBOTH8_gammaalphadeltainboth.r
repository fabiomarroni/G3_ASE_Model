#Extension of the Poisson-Gamma (PG) model with both.
#Instead of assuming a poisson sampling model 

#it considers a Negative Binomial model 
#See model_extension_withboth9.pdf section for model specification<- i stii need to update this pdf document
#This model/program was started September 2015 
#The correction bias for counts in Line1,q1, can be different from
#the correction bias for counts in Line2
#The mean structure is
#                         Tester           Line 1                       both
#Line 1(mated)            q11*(1/alpha)beta_i       q12*betai*alpha              tau[(1-q11)/alpha+(1-q12)alpha]betai
#Line 2(virgin)           q21*(1/delta)beta_igamma  q22*betai*delta*gamma        tau[(1-q11)/delta+(1-q12)*delta]betai*gamma
#It is required that q11+q12+q13=1  (notice that q13 is in the model implicitely)


# Get command line args with file names for processing
args <- commandArgs(TRUE)
        # c("./data/w52_fusion_example.csv","./data/test_out.csv")#   commandArgs(TRUE)
        #c("./data/split_10.csv","./data/split10_out.csv")#   commandArgs(TRUE)
         #c("./data/NB_simulated_data.csv","./data/NB_simulated_data_out.csv")
         #c("./data/input_rows_for_luis.csv","./data/input_rows_for_luis_out.csv")
         c("./fabioRprograms/test_M_noai_V_ai.csv","./fabioRprograms/out_test_M_noai_V_ai.csv")     
   # Import NB Functions
#setwd("C:/Users/Luis Leon-Novelo/Dropbox/projects/RNAextended/")


#setwd("C:/Users/Luis Leon-Novelo/Documents/luis/projects/RNAextended/")


#source("./rprograms/NBmodel_withboth_library8_withtaualphaandinvalpha.R")
#source("./rprograms/NBmodel_withboth_library8_withtaualphaandinvalphagammaequal1.R")
 source("NBmodel_withboth_library8_withtaualphaandinvalphagammaANDtauequal1.r")

# Prepare output File
fileout = args[2]
fornames1=c("TesterinMated","TesterVirgin")# 
fornames2=c("sampleprop","mean","q025","q975","AI_Bayes_pval")
fornames=paste(paste(rep(fornames1,each=length(fornames2)),fornames2,sep="_"),collapse=",")
fornames3=paste("AI_Mated_decision","AI_Virgin_decision","AI_diffinVirginandMated",sep=",")
firstheaders=paste(c("line","fusion_id","M_num_reps","V_num_reps","counts_M_tester","counts_M_line","counts_M_both",
                     "counts_V_tester","counts_V_line","counts_V_both",
                     "qsim_tester","qsim_line","chisqpvalue","independence_Bayesian_pvalue",
                     "prob_alphagreater1greaterdelta","prob_deltagreater1greateralpha"
                     ),collapse=",")

headers_out=paste(firstheaders,fornames,fornames3,sep=",")
headers_out=paste(headers_out,"alpha_postmean,delta_postmean",sep=",")
cat(headers_out,file=fileout,append=FALSE,sep="\n")

# Make Connection to input
con = file(args[1],"r")
newline=readLines(con,n=1) #moving the pointer away from the headers
headers_in=strsplit(newline,split=",")[[1]]

# Begin Processing
print(paste("saving results in ",fileout))
activeflag1="M_flag_analyze"#The flag that I will consider to select the exons
activeflag2="V_flag_analyze"
xactiveflag1=which(headers_in == activeflag1)
xactiveflag2=which(headers_in == activeflag2)

while(length(newline) != 0 ){
	flaganalyze=0
	while(flaganalyze==0){
	newline<-readLines(con,n=1);if(length(newline) == 0){break};
	mydata<-as.vector(strsplit(newline,split=",")[[1]])
	names(mydata)=headers_in[1:length(mydata)]
	flaganalyze=as.numeric(mydata[xactiveflag1])*as.numeric(mydata[xactiveflag2]);
	}
  if(length(newline) == 0){break}
	print("----------------------")
	print(paste("line:",mydata["line"]))
	print(paste("fusion_id",mydata["fusion_id"]))
	
	I1=as.numeric(mydata["M_num_reps"])		#Mates 	are line1
	I2=as.numeric(mydata["V_num_reps"])		#Virgin are line2
	xT=c(paste("M_tester_total_",1:I1,sep=""),paste("V_tester_total_",1:I2,sep=""));
	xL=c(paste("M_Line_total_",1:I1,sep=""),paste("V_Line_total_",  1:I2,sep=""));
	xB=c(paste("M_Both_total_",1:I1,sep=""),paste("V_Both_total_",  1:I2,sep=""));
	xs	<-	as.numeric(mydata[xT])#Vector of reads assigned to "tester" allele
	ys	<-	as.numeric(mydata[xL])	#Vector of reads assigned to "line" allele
	zs=	as.numeric(mydata[xB])	#Vector of reads assigned to both alleles
	q_=	c(as.numeric(mydata["qsim_tester"]),as.numeric(mydata["qsim_line"]),as.numeric(mydata["qsim_both"]))#Bias correction hyperparameter expected percentage of counts aligning to line when no AI
	      #REMEBER THAT SOME DATASETS REQUIRE Q to be multiplied by 2!!!!
  
  if(sum(xs+ys)>=1 && (q_[1]-q_[1]^2)!=0 && (q_[2]-q_[2]^2)!=0 && (q_[3]-q_[3]^2)!=0){
	  res<-
	    gibbs_NegBin_withboth(nsim=10^3,nburnin=10^3,nlag=10,
	                          xs=xs,    #Vector of reads assigned to "tester" allele I SWAPPED THE LABELS PUT THEM BACK
	                          ys=ys,  	#Vector of reads assigned to "line" allele
	                          zs=zs,		#Vector of reads assigned to both alleles
	                          is=c(1:I1,1:I2),		#index for biorep
	                          js=c(rep(1,I1),rep(2,I1)),	#index for Line1=1, vigin=Line2 
	                          sigma_alpha=0.25,
	                          a_beta=NA,#1,#0.01,<-------------Setting it to the MLE under the null model  (alpha=gamma=delta=1)
	                          #b_beta=0.01, 
	                          a_b_beta=NA,#2.1,#This parameters only if NBmodel_withboth_library5_gammaalphadeltainboth 
	                          b_b_beta=NA,#1,#where a priori for b_beta is used, this prior is vague
	                          a_gamma=0.01,
	                          b_gamma=0.01, 
	                          sigma_delta=0.25,					
	                          a_tau=0.1,
	                          b_tau=0.1,
	                          a_phi=2.01,#1,
	                          b_phi=0.05,#100,				#phi is apriori small, if inverse gamma is used prior mean b_phi/(a_phi-1)	       
	                          alpha_sd_MH=0.1,
	                          delta_sd_MH=0.1,    #
	                          phi_sd_MH=.05,			#Initial SD for the normal proposal to sample from phi via MH
	                          q1=2*q_	,  #2*q_[c(2,1,3)] when swapping xs and ys			#Bias correction hyperparameter for alleles in Line 1
	                          q2=2*q_,	 #2*q_[c(2,1,3)]			                        #Bias correction for alleles in Line 2
	                          plots=FALSE			#If true plots are depicted
	                          
	    )

	#binom.test(x=sum(res$counts[1,]), n=sum(res$counts), p=0.5, conf.level=0.95)	
	#Sorting the results to be saved  
  
out=
  paste(
	paste(
  mydata["line"],mydata["fusion_id"],I1,I2,paste(c(res$counts[1,],res$counts[2,]),collapse=","),
	q_[1],q_[2],
  round(res$independencetest$Chisqpvalue,4),round(res$independencetest$bayesianpvalue,4),
	paste((res$moretests)$alphagreater1greaterdelta,(res$moretests)$deltagreater1greateralpha,sep=","),
  paste(round(as.vector(t(res$tests[1,c(5,7,8,9,10)])),4),collapse=","),#length(vec1)
	paste(round(as.vector(res$tests[2,c(5,7,8,9,10)]),4),collapse=","),#lengt(vec2)
	paste(round(as.vector(c(res$tests[1,11],res$tests[2,11])),4),collapse=","),#lengt(vec2)
  ifelse(res$independencetest$bayesianpvalue<0.05,1,0),
  paste(round(unlist(res$modelparameters_postmean[c("alpha","delta")]),4),collapse=","),
	sep=","))


	}else{
		#If all the counts are zero I fill it with NaN
	out=paste(mydata["line"],mydata["fusion_id"],I1,I2,                         
		sum(xs[1:I1]),#"counts_M_tester"
		sum(ys[1:I1]),#"counts_M_line"
		sum(zs[1:I1]),#"counts_M_both"
		sum(xs[(I1+1):(I1+I2)]),#"counts_V_tester"
		sum(ys[(I1+1):(I1+I2)]),#"counts_V_line"
		sum(zs[(I1+1):(I1+I2)]),#"counts_V_both"
		q_[1],q_[2],	paste(rep(NA,19),collapse=","),sep=",")
		}
	#If one wants to see the results in a very long vector
	#(outmatrix=cbind(unlist(strsplit(headers_out, ",")),unlist(strsplit(out, ","))))
#sum(as.numeric(outmatrix[52:55,2]))
	cat(out,file=fileout,append=TRUE,sep="\n")
}
close(con)
