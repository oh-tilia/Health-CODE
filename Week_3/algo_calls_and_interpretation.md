**Call details from C5.0 algorithm**
Call:
C5.0.formula(formula = cell_line ~ ., data = train)

C5.0 [Release 2.07 GPL Edition]  	Thu Jun 18 13:15:30 2026
-------------------------------

Class specified by attribute `outcome'

Read 42 cases (499 attributes) from undefined.data

Number of inner nodes:    4
Number of terminal nodes: 5


Evaluation on training data (42 cases):

	    Decision Tree   
	  ----------------  
	  Size      Errors  

	     4    1( 2.4%)   <<

	   (a)   (b)   (c)   (d)   (e)    <-classified as
	  ----  ----  ----  ----  ----
	    12                            (a): class HS578T
	                       1          (b): class MCF10A
	                 3                (c): class MCF7
	                      13          (d): class MDAMB231
	                            13    (e): class SUM159

	Attribute usage:

	100.00%	C14orf132
	 66.67%	MYH14
	 59.52%	MMP17

	 **Decision tree:**

	C14orf132 <= 1.680699: MDAMB231 (14/1)
	C14orf132 > 1.680699:
	:...MYH14 > 0.3644396: MCF7 (3)
		MYH14 <= 0.3644396:
		:...MMP17 <= 2.448398: HS578T (12)
			MMP17 > 2.448398: SUM159 (13)

	**Decision tree from party_model:**

	Fitted party:
	[1] root
	|   [2] C14orf132 <= 1.37143
	|   |   [3] state in Cancerous: MDAMB231 (n = 12, err = 0.0%)
	|   |   [4] state in Healthy: MCF10A (n = 3, err = 0.0%)
	|   [5] C14orf132 > 1.37143
	|   |   [6] IGFBP3 <= 6.2422
	|   |   |   [7] A2ML1 <= 0.05261: SUM159 (n = 12, err = 0.0%)
	|   |   |   [8] A2ML1 > 0.05261: MCF7 (n = 2, err = 0.0%)
	|   |   [9] IGFBP3 > 6.2422: HS578T (n = 13, err = 7.7%)

Time: 0.1 secs



**Call details from Ranger**

Call:
 ranger::ranger(x = maybe_data_frame(x), y = y, mtry = min_cols(~floor(sqrt(ncol(train) -      1)), x), num.trees = ~700, importance = ~"permutation", num.threads = 1,      verbose = FALSE, seed = sample.int(10^5, 1), probability = TRUE) 

Type:                             Probability estimation 

Number of trees:                  700 

Sample size:                      43 

Number of independent variables:  498 

Mtry:                             22 

Target node size:                 10 

Variable importance mode:         permutation 

Splitrule:                        gini 

OOB prediction error (Brier s.):  0.01271114 
