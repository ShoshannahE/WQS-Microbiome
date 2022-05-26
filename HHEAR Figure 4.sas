********************************
**graphics for the weights in michels HMP**
********************************;


**import data;

	*keep the r output first column that numbers the rows--this is what the var1 and the seed vars become; 
proc import datafile="C:\Users\bixbym01\OneDrive - The Mount Sinai Hospital\HMP Microbiome\Michels Microbiome\Results\Michels_BMI_Negcont_WtsRSRH30_051922.csv"
dbms=csv out=weights replace; run;

*import the OTU taxonomic information; 
proc import datafile='C:\Users\bixbym01\OneDrive - The Mount Sinai Hospital\HMP Microbiome\Michels Microbiome\Final Processed Data\micheals_HMP_tax_table.csv'
out=tax
dbms=csv replace;
run;
data tax;
set tax;
rename var1=asv;
run;/*
data weights (DROP=VAR1);
set weights;
seed=input(VAR1,3.);
run;*/
data weights ;
set weights;
seed=VAR1;
drop VAR1;
run;
data vars;
set weights;
drop seed;
run;
%global slist;
proc sql;
      select name into :slist separated by ' '
      from dictionary.columns
      where memname=upcase("vars");
quit;
%put &slist;

proc transpose data=weights out=test2; 
var &slist; 
by seed; run; 

data test2; 
set test2;
rename col1=average_weight _name_=ASV;
run;

proc sort data=test2; by asv; run;


	proc means data=test2 uclm lclm noprint;
	by asv;
	var average_weight; output out=test3 mean=weight_rep; run;

	
data genera;
set tax;
keep asv phylum genus; 
run;
proc sort data=genera; by asv;
proc sort data=test3; by asv; 
data test3;
merge genera test3; 
by asv;
if weight_rep=. then delete;
if phylum="Firmicutes" and genus="NA" then genus="Firmicutes Unclassified";
if phylum="Actinobacter" and genus="NA" then genus="Actinobacter Unclassified";
run;


proc sql;
    create table want as select
    phylum,
	genus,
    sum(weight_rep) as sumweight_rep
    from test3
    group by genus;
quit;
proc sort data=want nodupkey out=want2; by genus; run; *48 unique genera; 


data want3;
set want2; 
threshold=(1/48);
if sumweight_rep > threshold then flag=1; *0.021;
else flag=0;
run;

**need to do something about the asvs--not sure that all the NAs are the same......
;proc freq data=want3; where flag=1; table genus*phylum; run;

proc freq data=tax; where genus="NA"; table phylum asv; run; *there are 100 asvs; 
**now do this for the sum within genera PER seed? ;
proc transpose data=weights out=perseed; 
var &slist; 
by seed; run;

**sum within genus;

data perseed; 
set perseed;
rename col1=average_weight _name_=asv;

run;
proc sort data=perseed; by asv;
proc sort data=tax; by asv;
data perseed2;
merge perseed tax; by asv;
if average_weight=. then delete;
keep seed asv average_weight phylum genus;
if phylum="Firmicutes" and genus="NA" then genus="Firmicutes Unclassified";
if phylum="Actinobacter" and genus="NA" then genus="Actinobacter Unclassified";
run;


proc sql;
    create table perseed3 as select
    phylum,
	genus, 
	seed,
    sum(average_weight) as sumaverage_weight
    from perseed2
    group by genus, seed;
quit;


proc sort data=perseed3 nodupkey out=perseed4; by genus seed; run; *89 genera here;
**just keep only the values that were above the trehsold;
proc sort data=perseed4; by genus; 
proc sort data=want3; by genus;   
data merged; 
merge perseed4 want3; by genus;
run;
data seedlevel;
set merged;
keep phylum genus seed sumaverage_weight flag; 
run;
data overall;
set merged;
keep phylum genus sumweight_rep flag; 
run;
proc sort data=overall nodupkey out=overall2; by genus; run;

data overall3;
set overall2;
where flag=1;
run;


data final; 
set overall2 seedlevel;
where flag eq 1;
run;
proc sort data=final nodupkey out=names; by genus; run;

proc sql;
alter table work.final
  modify TextVar char(36) format=$36.;
quit;
proc contents data=final; run;
data final;
length genus $36.;
format genus $36.;
set final;
if genus="Firmicutes Unclas" then genus="Firmicutes Unclassified";
if genus="Actinobacter Uncl" then genus="Actinobacter Unclassified";
if phylum="Bacteroidete" and genus="NA" then genus="Bacteroidete Unclassified";
run;
proc print data=names; var phylum genus; run;

*to create the final summed weights with the full taxonomic naming informaiton;

proc sort data=final; by phylum genus; 
proc sort data=tax; by phylum genus;

data tax_final;
merge final tax; by phylum genus; 
if flag ne 1 then delete;
run;


ods graphics / reset;
	ods path work.testtemp(update) sasuser.templat(update) sashelp.tmplmst(read);
proc template;
	define statgraph weights;
	dynamic TITLE CHEMICAL_DISTRIBUTION WEIGHT_DISTRIBUTION CHEMICAL  WEIGHTREP /*PERCENT_BAD*/ CHEMTYPE CHEMTYPEDIST;
	begingraph / pad=5 attrpriority=color;
	*entrytitle /*"% Weights from OTUs were 90% of mean weights were > threshold over 100 repeated holdouts" */
			"BMI: Sum of OTU mean weights within genera across 30 Repeated Holdouts >threshold (sorted by phylum)";
	entryfootnote textattrs=(size=8 family="Times New Roman") halign=left "Notes: Data points indicate sum of weights within genera for each of the 30 holdouts. Box plots show 25th, 50th, and 75th percentiles, and whiskers show 10th and 90th percentiles of sum of weights within genera for the 30 holdouts. Closed diamonds show sum of mean weights within genera for the 30 holdouts.";
		layout overlay / xaxisopts=(display=(label tickvalues) label="Genus" labelattrs=(Size=11 weight=bold family="Times New Roman") tickvalueattrs=(Size=8 family="Times New Roman") discreteopts=(sortorder=data))
 						 yaxisopts=(label="			Weight" labelattrs=(Size=11 weight=normal family="Times New Roman") 
									linearopts=(viewmin=0 viewmax=0.60 includeranges=(0-.6) tickvaluelist=(0 .1 .2 .3 .4 .5 .6 )) tickvalueattrs=(Size=10 family="Times New Roman"));
						 /*y2axisopts=(label="		# Reps Above Chemical of Concern Threshold" labelattrs=(Size=11 weight=normal family="Times New Roman") 
									linearopts=(viewmin=0 viewmax=30 includeranges=(0-30) tickvaluelist=(0 5 10 15 20 25 30) tickdisplaylist=('0' '5' '10' '15' '20' '25' '30')) tickvalueattrs=(Size=10 family="Times New Roman"));*/
		/*These are plotted underneath only in order to have correct legend;
					*Plot Mean Weight of 100 Reps;*/
					scatterplot y=WEIGHTREP 
								x=CHEMICAL / 
									yaxis=y 
									name="valid"
									group=CHEMTYPE
									groupdisplay=cluster
									clusterwidth=0.85
									markerattrs=(color=black size=6pt symbol=diamondfilled) datalabelattrs=(size=4)  ; 
			/**Bar Chart for %Times Bad Actors over 100 Reps;
			barchartparm y=PERCENT_BAD 
						 x=CHEMICAL / 
						 	yaxis=y2 
							name="percent"
							group=CHEMTYPE
							groupdisplay=cluster
							display=(fill) 
							fillattrs=(transparency=0.8)
							baselineattrs=(thickness=0);*/
			*Box Plot for Distribution of Weights over 100 Reps;
			boxplot y=WEIGHT_DISTRIBUTION 
					x=CHEMICAL_DISTRIBUTION / 
							yaxis=y 
							name="box"
							group=phylum
							groupdisplay=cluster
							display=(fill median) 
							datatransparency=0.3 
							boxwidth=0.9
							clusterwidth=0.85
							fillattrs=(transparency=0.4)
							outlineattrs=(thickness=0.5) 
							medianattrs=(thickness=1.5)
							whiskerattrs=(thickness=1.5) 
							whiskerpercentile=10 
							capscale=0.5; 
			*Plot Weights for Each of 100 Reps;
			scatterplot y=WEIGHT_DISTRIBUTION 
						x=CHEMICAL_DISTRIBUTION / 
							yaxis=y 
							name="30 weights"
							group=CHEMTYPEDIST
							datatransparency=0.75
									groupdisplay=cluster
									clusterwidth=0.85
							markerattrs=(size=2pt symbol=circlefilled) datalabelattrs=(size=4)  
							jitter=auto;
			*Plot Mean Weight of 100 Reps;
			scatterplot y=WEIGHTREP 
						x=CHEMICAL / 
							yaxis=y 
							group=CHEMTYPE
								groupdisplay=cluster
									clusterwidth=0.85
							markerattrs=(size=4pt symbol=diamondfilled);
			*Bad actor threshold;
			referenceline y=0.011 /  
							yaxis=y 
							lineattrs=(color=black)
							curvelabel="Threshold (0.02)" 
							curvelabelattrs=(Size=8 color=black family="Times New Roman") 
							curvelabelposition=max;
			*Symbol legend for mean weights;
			/*discretelegend "full" "valid" / 
							title="Mean Weight:" 
							titleattrs=(weight=bold size=9 family="Times New Roman")
							valueattrs=(size=9 family="Times New Roman")
							autoalign=(topright) 
							location=inside
							across=1;*/
			*Color legend for chemical groups;
			discretelegend "box" / 
							title="Phylum:" 
							titleattrs=(weight=normal size=9 family="Times New Roman")
							valueattrs=(size=9 family="Times New Roman");
		endlayout;
	endgraph;
	end;
	define style noheaderborder;
    	parent = styles.default;
		*Define colors for chemical groups;
		class Graphdata1 / color=DEPPK contrastcolor=DEPPK;
		class Graphdata2 / color=BIPB contrastcolor=BIPB;
		class Graphdata3 / color=blue contrastcolor=blue;
		class Graphdata4 / color=VIBG contrastcolor=VIBG;
		class Graphdata5 / color=green contrastcolor=green;
		*Eliminate excessive borders by making them white or turning off;
    	class graphborderlines / contrastcolor=white;	
        class graphbackground / color=white ;
		style graphwalls from graphwalls / frameborder=off;
	end;
	run;
 
	ods html image_dpi=400;
ods graphics on / width=225mm attrpriority=none border=off;
ods html style=noheaderborder;
proc sgrender data=final template=weights;
	dynamic /*CHEMTYPEDIST="Typedist"*/ CHEMICAL_DISTRIBUTION="genus" WEIGHT_DISTRIBUTION="sumaverage_weight"  
			/*CHEMTYPE="type"*/ CHEMICAL="genus" WEIGHTREP="sumweight_rep" /*PERCENT_BAD="times_bad"*/ ;
	run;
	ods graphics off;
