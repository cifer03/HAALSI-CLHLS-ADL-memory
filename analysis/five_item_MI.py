"""MI of binary adjustment variables ONLY, for known mood/memory and outcome.
Two explicit scenarios for CL: ordinary missing with unknown responses excluded,
and broader MAR sensitivity also treating 'don't know' as missing covariates.
No cognitive inability/proxy score, outcome or death imputation.
"""
from pathlib import Path
import json
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.special import expit
from scipy.stats import t as student_t
R=Path(__file__).resolve().parent
import pickle
with open(R/'five_item_frames.pkl','rb') as stream: frames=pickle.load(stream)
ha4=frames['HAALSI'];cl=frames['CLHLS']
fixed=['mood_z','memory_low_z','agec','agec2','male'];covs=['school_any','hypertension','diabetes','cvd','smoke','drink','partner'];cols=fixed+covs;keys=['mood_z','memory_low_z'];rng=np.random.default_rng(20261011)
out={'models':[],'imputations':40,'iterations_per_chain':20,'seed':20261011,'limitations':['Known memory, mood and functional endpoint only. No cognitive test noncompletion, proxy score, death or unknown endpoint imputation.','Binary logistic FCS with asymptotic normal parameter draws, includes observed functional outcome to preserve associations.','Conditional MAR for covariates is an assumption; unknown disease may be related to access to diagnosis or poor health. Broad MAR is a stronger sensitivity, not correction proof.','Does not recover all baseline people or eliminate informative baseline score noncompletion, survival or loss.','Imputation predictors restricted to analysis variables and known outcome, no rich auxiliary clinical or socioeconomic variables. Preliminary sensitivity, not final publication MI specification.','Fixed five-item baseline-reference within-cohort score SD, robust completed-data SE, Rubin pooling with finite-sample Barnard-Rubin df.','No model or threshold selection, no project-wide multiplicity correction; conditional ORs not causal effects.']}
def fit(p):
 X=sm.add_constant(p[cols],has_constant='add').astype(float);assert np.linalg.matrix_rank(X)==len(X.columns)
 f=sm.GLM(p.end,X,family=sm.families.Binomial()).fit(cov_type='HC0',maxiter=100);assert f.converged;return f
def simple_effects(f):
 return {k:{'OR':float(np.exp(f.params[k])),'CI':np.exp([f.params[k]-1.959964*f.bse[k],f.params[k]+1.959964*f.bse[k]]).tolist()} for k in keys}
for name,q in [('HAALSI W2-W4',ha4.copy()),('CLHLS2018-2021',cl.copy())]:
 risk=q[q.start.eq(0)].copy();s=risk[['mood','memory']].dropna();risk['mood_z']=(risk.mood-s.mood.mean())/s.mood.std();risk['memory_low_z']=(s.memory.mean()-risk.memory)/s.memory.std()
 known=risk.dropna(subset=['end']+fixed).copy();cc=known.dropna(subset=covs)
 scenarios=[('all_covariate_missing_MAR',known)]
 if name.startswith('CL'):
  # Retain raw genuine missing; exclude explicit don't-know ONLY if that derived predictor is missing.
  flag=(known.school_any.isna()&(known.f1.eq(88)|(~known.f1.between(0,30)&~known.f1.isin([88,99])&known.f1.notna())))
  flag|=known.hypertension.isna()&known.g15a1.eq(8);flag|=known.diabetes.isna()&known.g15b1.eq(8)
  flag|=known.cvd.isna()&(known[['g15c1','g15d1']].eq(8).any(axis=1))
  # Lifestyle unknown coding remains excluded conservatively if not a NaN/ordinary 9 missing.
  for d,rawcol,valid in [('smoke','d71',[1,2]),('drink','d81',[1,2]),('partner','f41',[1,2,3,4,5])]:
   flag|=known[d].isna()&known[rawcol].notna()&~known[rawcol].isin(valid+[9])
  scenarios=[('ordinary_missing_excluding_explicit_unknown',known[~flag]),('all_covariate_missing_MAR',known)]
 for label,pool in scenarios:
  p=pool.reset_index(drop=True);missing=p[covs].isna();missingcols=[k for k in covs if missing[k].any()]
  assert not p[fixed+['end']].isna().any().any()
  for k in covs:assert set(p[k].dropna().unique())<={0.,1.}
  pp=p.dropna(subset=covs);fc=fit(pp)
  r={'cohort':name,'scenario':label,'n':len(p),'events':int(p.end.eq(1).sum()),'complete_case_n':len(pp),'complete_case_events':int(pp.end.eq(1).sum()),'recovered_by_covariate_imputation_n':len(p)-len(pp),'excluded_from_known_scores_outcomes_n':len(known)-len(p),'missing_by_variable':{k:int(missing[k].sum()) for k in covs},'complete_case_effects':simple_effects(fc),'imputed_variables':missingcols}
  estimates=[];variances=[];rates=[];traces=[]
  for m in range(40):
   complete=p.copy()
   for k in missingcols:
    mask=missing[k];complete.loc[mask,k]=rng.binomial(1,float(p[k].mean()),int(mask.sum()))
   chain=[]
   for iteration in range(20):
    for k in missingcols:
     mask=missing[k].to_numpy();others=[v for v in cols if v!=k]+['end'];Z=sm.add_constant(complete[others],has_constant='add').astype(float).to_numpy()
     f=sm.GLM(p.loc[~mask,k],Z[~mask],family=sm.families.Binomial()).fit(maxiter=100);assert f.converged
     params=rng.multivariate_normal(np.asarray(f.params),np.asarray(f.cov_params()),check_valid='raise')
     complete.loc[mask,k]=rng.binomial(1,expit(Z[mask]@params))
    chain.append({k:float(complete.loc[missing[k],k].mean()) for k in missingcols})
   final=fit(complete);estimates.append([float(final.params[k]) for k in keys]);variances.append([float(final.bse[k]**2) for k in keys]);rates.append(chain[-1]);traces.append(chain)
   if (m+1)%10==0:print(name,label,'MI',m+1,flush=True)
  E=np.array(estimates);U=np.array(variances).mean(axis=0);B=E.var(axis=0,ddof=1);T=U+(1+1/40)*B;mean=E.mean(axis=0);lam=(1+1/40)*B/T
  olddf=np.where(lam>0,39/np.maximum(lam,1e-20)**2,np.inf);dfcomplete=len(p)-(len(cols)+1);obsdf=((dfcomplete+1)/(dfcomplete+3))*dfcomplete*(1-lam);df=1/(1/olddf+1/obsdf);critical=student_t.ppf(.975,df)
  r['pooled_effects']={k:{'OR':float(np.exp(mean[j])),'CI':np.exp([mean[j]-critical[j]*np.sqrt(T[j]),mean[j]+critical[j]*np.sqrt(T[j])]).tolist(),'Rubin_Barnard_df':float(df[j]),'between_variance':float(B[j]),'within_variance':float(U[j]),'relative_increase_variance':float((1+1/40)*B[j]/U[j]),'MC_error_logOR':float(np.sqrt(B[j]/40)),'pooled_SE_logOR':float(np.sqrt(T[j]))} for j,k in enumerate(keys)}
  r['imputed_binary_prevalence_range']={k:[min(v[k] for v in rates),max(v[k] for v in rates)] for k in missingcols};r['chain_prevalence_traces']=traces
  out['models'].append(r);print(json.dumps({k:v for k,v in r.items() if k!='chain_prevalence_traces'},ensure_ascii=False),flush=True)
(R/'five_item_MI_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
