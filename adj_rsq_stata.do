
*======================================
*
*     regression for EB and IE model: calibrate parameters and calculate adjusted R-squared
*		  
*======================================

  global PP   "/Users/HaoLI/Dropbox/Network/Network2019/QJ"  // please set your own root directory
  global D    "$PP/data"      //data
  cd "$D"
  import delimited indegree_ccd.csv, clear

** network where isolated nodes are removed
* total number of nodes = 83663
* total number of isolated nodes = 0
* number of edges = 275764
* average outdegree = 3.296129
* secondtier = ({p}+(mu-1)*{r})

* network with isolated nodes
* total number of nodes = 129003
* total number of isolated nodes = 45340, Pc = 1-45340/129003 = 0.6485
* number of edges = 351713
* average outdegree = 2.276

*-> EB model using network where isolated nodes are removed
generate lgh_kpmul_kploop_rmisonodes_3296 = log10(gh_kpmul_kploop_rmisonodes_3296)
nl (lgh_kpmul_kploop_rmisonodes_3296 = (3.296129/(15*{eta=0.1}))*log10(((3.296129*{q=0.1}/{eta})/(indegree_8 + 3.296129*{q}/{eta}))))
est store reg_EB

*-> IE model using network with isolated nodes
generate l10gh_kpmul_kploop_2726 = log10(gh_kpmul_kploop_2726) 
nl (l10gh_kpmul_kploop_2726 = (351713/(15*({p}+(2.726-1)*{r=0.0001})*129003))*log10(((0.6485*2.726*{q=0.01}/({p=0.01}+(2.726-1)*0.01) + (1-0.6485)*({q}-(2.726-1)*0.01)/({p}+(2.726-1)*{r})))/(indegree_1 + (0.6485*2.726*{q}/({p}+(2.726-1)*{r}) + (1-0.6485)*({q}-(2.726-1)*0.01)/({p}+(2.726-1)*{r})))))


*-> JR model using gh_kpmul_kploop_rmisonodes_3296, 
import delimited indegree_ccd.csv, clear
generate lgh_kpmul_kploop_rmisonodes_3296 = log10(gh_kpmul_kploop_rmisonodes_3296)
nl (lgh_kpmul_kploop_rmisonodes_3296 = (1+{r=0.1})*log10(({r}*3.296)/(indegree_8+{r}*3.296)))

*-> JR model using gh_kpmul_kploop_2726
import delimited indegree_ccd.csv, clear
generate l10gh_kpmul_kploop_2726 = log10(gh_kpmul_kploop_2726) 
nl (l10gh_kpmul_kploop_2726 = (1+{r=0.1})*log10(({r}*2.726)/(indegree_1+{r}*2.726)))


esttab reg_EB reg_IE using indegreeCDF_IE_EB.tex, replace
