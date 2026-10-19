"""Conservative score-validity check and independent fresh-seed FCS MI.
Keeps all raw-unknown HA recall items unknown; no participant-level output.
"""
from independent_audit import *
old=json.loads((O/'independent_audit_results.json').read_text());out={'score_rule':'Conservative alternative: require all ten HA delayed word indicators valid 0/1; invalid/unknown derived zeros not accepted. No assertion that publisher score is erroneous.','cohorts':[]}
for name,q in frames.items():
 q=q.copy()
 if name=='HAALSI':
  word=q[[f'w2cn012s{i}' for i in range(1,11)]];mask=~word.isin([0,1]).all(axis=1);q.loc[mask,'memory']=np.nan;q=std(q);frames[name]=q
 main=model(q[q.start.eq(0)]);r={'cohort':name,'flow':flow(q),'main':main,'sensitivities':{}}
 if name=='HAALSI':
  r['sensitivities']={'age65_89':model(q[q.start.eq(0)&q.agec.between(-1,1.4)]),'all_six_independent_at_baseline_five_item_endpoint':model(q[q.six_start.eq(0)]),'cohort_specific_all_six':model(q[q.six_start.eq(0)],'six_end'),'five_help_or_cannot':model(q[q.broad_start.eq(0)],'broad_end'),'exclude_followup_proxy':model(q[q.start.eq(0)&q.w4rproxy.eq(0)])}
 out['cohorts'].append(r)
 print(name,main,flush=True)
(O/'score_quality_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
# Independent MI uses its own NumPy logistic solver, independently reconstructed data,
# fresh seed 20261017. Pooling draws retained for verification (aggregate coefficients only).
rng=np.random.default_rng(20261017);fixed=C[:5];covs=C[5:];miout={'seed':20261017,'m':40,'iterations':20,'models':[]}
for name,q in frames.items():
 pool=q[q.start.eq(0)].dropna(subset=fixed+['end']).copy();scenarios=[('all_covariate_missing_MAR',pool)]
 if name=='CLHLS':
  flag=pool.school_any.isna()&(pool.f1.eq(88)|(~pool.f1.between(0,30)&~pool.f1.isin([88,99])&pool.f1.notna()))
  flag|=pool.hypertension.isna()&pool.g15a1.eq(8);flag|=pool.diabetes.isna()&pool.g15b1.eq(8);flag|=pool.cvd.isna()&pool[['g15c1','g15d1']].eq(8).any(axis=1)
  for target,col,valid in [('smoke','d71',[1,2]),('drink','d81',[1,2]),('partner','f41',[1,2,3,4,5])]:flag|=pool[target].isna()&pool[col].notna()&~pool[col].isin(valid+[9])
  scenarios=[('ordinary_missing_excluding_explicit_unknown',pool[~flag]),('all_covariate_missing_MAR',pool)]
 for label,p in scenarios:
  values=p[C+['end']].to_numpy(float);missing=np.isnan(values);indices=[i for i in range(5,12) if missing[:,i].any()];E=[];U=[];trace=[];M=40
  for m in range(M):
   a=values.copy()
   for j in indices:a[missing[:,j],j]=rng.binomial(1,np.nanmean(values[:,j]),missing[:,j].sum())
   chain=[]
   for it in range(20):
    for j in indices:
     mm=missing[:,j];other=[i for i in range(13) if i!=j];Z=np.column_stack([np.ones(len(a)),a[:,other]]);b,_,cov,_=logistic(Z[~mm],a[~mm,j]);draw=rng.multivariate_normal(b,cov,check_valid='raise');a[mm,j]=rng.binomial(1,expit(Z[mm]@draw))
    chain.append({C[j]:float(a[missing[:,j],j].mean()) for j in indices})
   b,V,_,_=logistic(np.column_stack([np.ones(len(a)),a[:,:12]]),a[:,12]);E.append(b[[2,1]].tolist());U.append(np.diag(V)[[2,1]].tolist());trace.append(chain)
   if (m+1)%10==0:print(name,label,m+1,flush=True)
  e=np.array(E);u=np.array(U);B=e.var(axis=0,ddof=1);ubar=u.mean(axis=0);T=ubar+(1+1/M)*B;lam=(1+1/M)*B/T;dfc=len(p)-13;df=1/(1/((M-1)/lam**2)+1/((dfc+1)/(dfc+3)*dfc*(1-lam)));mean=e.mean(axis=0);se=np.sqrt(T)
  r={'cohort':name,'scenario':label,'n':len(p),'events':int(p.end.sum()),'effects':{},'coefficient_draws_order':K,'coefficient_draws':E,'HC0_variance_draws':U,'prevalence_traces':trace}
  for j,k in enumerate(K):r['effects'][k]={'OR':float(np.exp(mean[j])),'CI':np.exp(mean[j]+np.array([-1,1])*t.ppf(.975,df[j])*se[j]).tolist(),'p':float(2*t.sf(abs(mean[j]/se[j]),df[j])),'df':float(df[j]),'SE':float(se[j]),'MC_error_logOR':float(np.sqrt(B[j]/M))}
  miout['models'].append(r);(O/'independent_MI_results.json').write_text(json.dumps(miout,ensure_ascii=False,indent=2));print(r['effects'],flush=True)
print('INDEPENDENT MI COMPLETE',flush=True)
