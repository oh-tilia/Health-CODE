Call details from C5.0 algorithm
Call:
C5.0.formula(formula = cell_line ~ ., data = train)

C5.0 [Release 2.07 GPL Edition]  	Thu Jun 18 13:15:30 2026
-------------------------------

Class specified by attribute `outcome'

Read 42 cases (499 attributes) from undefined.data

Decision tree:

C14orf132 <= 1.680699: MDAMB231 (14/1)
C14orf132 > 1.680699:
:...MYH14 > 0.3644396: MCF7 (3)
    MYH14 <= 0.3644396:
    :...MMP17 <= 2.448398: HS578T (12)
        MMP17 > 2.448398: SUM159 (13)


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

Time: 0.1 secs
