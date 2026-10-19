from independent_audit import *
rng=np.random.default_rng(20261018);fixed=C[:5];covs=C[5:];miout={'seed':20261018,'m':40,'iterations':20,'models':[]}
for name,q in [('HAALSI',frames['HAALSI'])]:
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
  miout['models'].append(r);(O/'original_score_HA_MI_results.json').write_text(json.dumps(miout,ensure_ascii=False,indent=2));print(r['effects'],flush=True)
print('INDEPENDENT MI COMPLETE',flush=True)
