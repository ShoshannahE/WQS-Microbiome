************************************************************;
** read in data in form of relative abundance;
proc import 
	file='C:\...\HMP_Pre_processed_Stool_1stvisit_10%Detect_RELATIVE_Abundances.csv'
     	out=one
	 dbms=csv replace;

data one; 
    set one; 
    dummy=1;  * variable for merging means;

*** 20 OTUs selected based on literature related to smoking;
*** set 2 as strongly related, 8 as moderately related, and 10 as weakly related;
****
 strong:  x425, x450
 medium:  x545, x616, x845, x58, x63, x133, x141, x235 
 low:     x378, x584, x737, x6, x127, x202, x239, x262, x359, x391 

** order macros to link selected OTUs with numbered variables;
%macro selOTUs; OTU_97_17 OTU_97_173 OTU_97_189  OTU_97_201 OTU_97_244 OTU_97_1074 OTU_97_1084 OTU_97_1200
                OTU_97_1214 OTU_97_1367 OTU_97_1634 OTU_97_1960 OTU_97_2220 OTU_97_10052 OTU_97_11868
				OTU_97_13013 OTU_97_13765 OTU_97_14182 OTU_97_15976 OTU_97_16504 %mend;
%macro selected; x425 x450 x545 x616 x845 x58 x63 x133 x141 x235 x378 x584 x737 x6 x127 x202 x239
					x262 x359 x391 %mend;
proc means data=one;
   var %selotus %selected; run;
proc means data=one noprint;
   var %selOTUs;
   id dummy;
   output out=max max=max1 max2 max3 max4 max5 max6 max7 max8 max9 max10
                      max11 max12 max13 max14 max15 max16 max17 max18 max19 max20;
		run;


data simul;
   merge one (drop= drop psn  run_center hmp_body_side hmb_body_subsite srs_sample visitno) 
	max; by dummy;
   b0 = -5;
   beta_s = 8; 	** beta for strong association;
   beta_m = 4; 	** beta for moderate association;
   beta_w = 2; *	** beta for weak association;
   bsex = -1;
   if sex='Male' then female=0;
   if sex='Female' then female=1;

   array otu %otuname;
   array xname %xnamefull;
   do over otu;
      xname = log10(otu+1);
	  end;
*case 1;   
sum1 = sum(log2(x425/log10(max1+1)+1), log2(x450/log10(max2+1)+1));
sum2 = sum(log2(x545/log10(max3+1)+1), log2(x616/log10(max4+1)+1), log2(x845/log10(max5+1)+1),
                          log2(x58/log10(max6+1)+1),  log2(x63/log10(max7+1)+1), 
                          log2(x133/log10(max8+1)+1), log2(x141/log10(max9+1)+1), log2(x235/log10(max10+1)+1));
sum3 = sum(log2(x378/log10(max11+1)+1), log2(x584/log10(max12+1)+1), log2(x737/log10(max13+1)+1), 		  		       log2(x6/log10(max14+1)+1), log2(x127/log10(max15+1)+1), 
                       log2(x202/log10(max16+1)+1), log2(x239/log10(max17+1)+1), log2(x262/log10(max18+1)+1),   
                       log2(x359/log10(max19+1)+1),  log2(x391/log10(max20+1)+1)) ;

logit = min(b0 + beta_s*sum1 + beta_m*sum2 + beta_w*sum3 + bsex*female,
            3.5);

   mu = 1/(1+exp(-logit));

  do sim_exp= 1 to 10;  ** run multiple times to evaluate similarity across simulated sets;
       seed=100597 * sim_exp;
	   call ranbin(seed, 1, mu ,sim_smk);
       rv0 = (ranuni(100597*sim_exp)<=0.5);
     output;
	 end;
   keep sim_exp sim_smk female logit rv0 x1--x868 max1--max20 rsid sum1 sum2 sum3;
   run;
proc sort data=simul; by sim_exp;
proc univariate data=simul;
   var logit ;*sum1 sum2 sum3;
   histogram; run; quit;
proc freq data=simul;
   table sim_smk*sim_exp; run;
proc means data=simul;
 	by sim_exp;
 	var logit; run;

%macro xvar; x1--x868 %mend;
%macro qvar; q1 - q868 %mend;
%macro truth; x425 x450 x545 x616 x845 x58 x63 x133 x141 x235 x359  x391 x378 x584 x737 x6 x127 x202 x239 x262  %mend;
%macro qtruth;q425 q450 q545 q616 q845 q58 q63 q133 q141 q235 q359  q391 q378 q584 q737 q6 q127 q202 q239 q262  %mend;

** to rank  keeping zeros as 0 and then tertile the remaining values – and then combine;
data forrank;
   set simul;
   array xvar %xvar;
   do over xvar;
     if xvar=0 then xvar=.;
   end; 
proc rank data=forrank  out=ranks groups=3;
   var %xvar;
   ranks  %qvar;
run;
proc freq data=ranks;
   table q1 q2 q3 q425;
run;
data ranks;
   set ranks;
   array qvar %qvar;
   do over qvar; 
      qvar = qvar+1;
      if qvar=. then qvar=0;
	  end;
proc freq data=ranks;
   where sim_exp=1;
   table %qtruth;
run;

proc univariate data=simul noprint;
   var %truth;
   histogram;
   run;
proc genmod data=ranks descending;
*  where sim_exp=1;
   by sim_exp;
   model sim_smk = %qtruth q1--q5 female/dist=bin;
   ods output parameterestimates=output;
   run;
proc print data=output; by sim_exp; run;

proc sort data=simul; by sim_exp rsid;
proc means data=simul n sum mean min max;
   by sim_exp;
   var sim_smk logit rv0 ;*%selected;
   run;

data forR;
   set ranks; drop %xvar max1--max20 sum1 sum2 sum3; run;
proc contents data=forR varnum; run;

data forSHANNAH_070821;
   set forR;
   drop q1--q868;

proc export
	data= forSHANNAH_070821
	outfile='C:\Users\gennic01\Desktop\research\shoshannah eggers\hmp study\forShannah070821.csv'
	 dbms=csv replace;
run; 

