"""Independent raw-data reconstruction and numerical audit, 2026-10-17.
Does not execute/import earlier analysis scripts or use saved frames to construct data.
Saved frames/results are opened only AFTER raw reconstruction for comparison.
No individual records are written to deliverables. Requires numpy pandas scipy.
"""
from pathlib import Path
import zipfile,json,pickle,hashlib,os
import numpy as np,pandas as pd
from scipy.special import expit
from scipy.stats import norm,t
ROOT=Path(__file__).resolve().parent.parent; D=Path(os.environ.get('HAALSI_CLHLS_INPUTS', ROOT/'inputs')); O=ROOT/'audit'
C=['mood_z','memory_low_z','agec','agec2','male','school_any','hypertension','diabetes','cvd','smoke','drink','partner'];K=['memory_low_z','mood_z']
def readzip(z,name=None,columns=None):
 with zipfile.ZipFile(D/z) as f:
  name=name or next(n for n in f.namelist() if n.endswith('.dta'))
  with f.open(name) as x:return pd.read_stata(x,columns=columns,convert_categoricals=False)
def union(x):
 a=np.asarray(x,float);r=np.full(len(a),np.nan);r[np.all(a==0,axis=1)]=0;r[np.any(a==1,axis=1)]=1
 return pd.Series(r,index=x.index)
def logistic(X,y):
 X=np.asarray(X,float);y=np.asarray(y,float);b=np.zeros(X.shape[1]);ll=-np.inf
 assert np.linalg.matrix_rank(X)==X.shape[1]
 for it in range(100):
  p=expit(X@b);H=X.T@((p*(1-p))[:,None]*X);step=np.linalg.solve(H,X.T@(y-p));rate=1.
  old=np.sum(y*(X@b)-np.logaddexp(0,X@b))
  while True:
   bn=b+rate*step;new=np.sum(y*(X@bn)-np.logaddexp(0,X@bn))
   if new>=old-1e-10:break
   rate/=2
   if rate<1e-9:raise RuntimeError('line search')
  b=bn
  if np.max(np.abs(rate*step))<1e-10:break
 else:raise RuntimeError('nonconvergence')
 p=expit(X@b);bread=np.linalg.inv(X.T@((p*(1-p))[:,None]*X));scores=X*(y-p)[:,None];V=bread@(scores.T@scores)@bread
 assert np.max(np.abs(X.T@(y-p)))<1e-6
 return b,V,bread,it+1

def model(q,y='end'):
 p=q.dropna(subset=C+[y]);X=np.column_stack([np.ones(len(p)),p[C].to_numpy()]);b,V,_,it=logistic(X,p[y]);se=np.sqrt(np.diag(V));out={'n':len(p),'events':int(p[y].sum()),'iterations':it,'effects':{}}
 for k in K:
  j=C.index(k)+1;out['effects'][k]={'OR':float(np.exp(b[j])),'CI':np.exp(b[j]+np.array([-1,1])*1.959964*se[j]).tolist(),'p':float(2*norm.sf(abs(b[j]/se[j]))),'beta':float(b[j]),'SE':float(se[j])}
 return out

def std(q,start='start'):
 s=q.loc[q[start].eq(0),['mood','memory']].dropna();q=q.copy();q['mood_z']=(q.mood-s.mood.mean())/s.mood.std();q['memory_low_z']=(s.memory.mean()-q.memory)/s.memory.std();return q

def flow(q):
 a=q[q.start.eq(0)];b=a.dropna(subset=['mood','memory','agec','male']);c=b.dropna(subset=['end']);d=c.dropna(subset=C)
 return dict(age_eligible=len(q),baseline_positive=int(q.start.eq(1).sum()),baseline_unknown=int(q.start.isna().sum()),baseline_negative=len(a),scores_age_sex_known=len(b),scored_followup_known=len(c),scored_followup_unknown=len(b)-len(c),complete_case=len(d),excluded_covariates=len(c)-len(d),events=int(d.end.sum()))
H=['prim_key','w1c_bd_educ4','w2status','w3status','w2c_rage_calc','w2c_rsex','w2rproxy','w2c_cn_delrec','w2c_cd_cesd']+[f'w2cd{i}' for i in range(201,221)]+['w2'+k for k in ['cm002','cm007_males','cm007_females','cm025','cm029','cm038','cm041','cm059','cm061','cm068','cm069','c_bd_mar']]+[f'w2pf{i:03}' for i in [8,12,18,22,14,13,5,9,15,19,1,4]]+[f'w2cn012s{i}' for i in range(1,11)]+['w2cn012s88']
h=pd.read_stata(D/'HAALSI W1-3 Longitudinal Dataset.dta',columns=H,convert_categoricals=False)
w=readzip('dataverse_files.zip','HAALSI Wave 4 Public Dataset.dta',['prim_key','w4status','w4rproxy']+[f'w4pf{i:03}' for i in [8,12,18,22,14,13,5,9,15,19,1,4]])
assert h.prim_key.is_unique and w.prim_key.is_unique
h=h[h.w2status.eq(1)&h.w2c_rage_calc.between(65,120)].merge(w,on='prim_key',how='left',validate='one_to_one')
h['male']=h.w2c_rsex.map({1:1.,2:0.});h['school_any']=h.w1c_bd_educ4.map({1:0.,2:1.,3:1.,4:1.});h['agec']=(h.w2c_rage_calc-75)/10;h['agec2']=h.agec**2
h['hypertension']=h.w2cm002.map({1:1.,2:0.});h['diabetes']=h.w2cm007_males.where(h.male.eq(1),h.w2cm007_females).map({1:1.,2:0.}).where(h.male.notna());h['cvd']=union(h[['w2cm025','w2cm029','w2cm038','w2cm041']].replace({1:1.,2:0.}).where(lambda x:x.isin([0.,1.])))
for target,ever,current in [('smoke','w2cm059','w2cm061'),('drink','w2cm068','w2cm069')]:
 h[target]=h[current].map({1:1.,2:0.});h.loc[h[ever].eq(2),target]=0.
h['partner']=h.w2c_bd_mar.map({0:0.,1:0.,2:0.,3:1.})
items=h[[f'w2cd{i}' for i in range(201,221)]].where(lambda x:x.isin([1,2,3,4]));v=items.to_numpy();v[:,[3,7,11,15]]=5-v[:,[3,7,11,15]];h['mood']=pd.DataFrame(v-1,index=h.index).sum(axis=1,min_count=20)
h['memory']=h.w2c_cn_delrec.where(h.w2c_cn_delrec.between(0,10));cesd=int((h.mood.dropna()-h.w2c_cd_cesd).abs().gt(1e-8).sum())
for wave in [2,4]:
 x=h[[f'w{wave}pf{i:03}' for i in [8,12,18,22,14]]].apply(lambda s:s.map({1:1.,2:0.}));x.loc[h[f'w{wave}pf013'].eq(2)&h[f'w{wave}pf014'].isna(),f'w{wave}pf014']=0
 h[f'four{wave}']=union(x.iloc[:,:4]).where(h[f'w{wave}status'].eq(1));h[f'five{wave}']=union(x).where(h[f'w{wave}status'].eq(1))
 sixth=h[f'w{wave}pf004'].map({1:1.,2:0.});h[f'six{wave}']=union(pd.concat([h[f'five{wave}'],sixth],axis=1)).where(h[f'w{wave}status'].eq(1))
 for j,dp in enumerate([5,9,15,19,13]):x.loc[h[f'w{wave}pf{dp:03}'].eq(3),x.columns[j]]=1
 h[f'broad{wave}']=union(x).where(h[f'w{wave}status'].eq(1))
bad=h.w3status.eq(0)&h.w4status.eq(1)
h.loc[bad,['four4','five4','six4','broad4']]=np.nan
h['start']=h.five2;h['end']=h.five4;h['start4']=h.four2;h['end4']=h.four4;h['six_start']=h.six2;h['six_end']=h.six4;h['broad_start']=h.broad2;h['broad_end']=h.broad4
c=readzip('clhls_2018_cross_sectional_dataset_15874.zip',columns=['id','trueage','a1','f1','g15a1','g15b1','g15c1','g15d1','d71','d81','f41','yearin']+[f'b3{i}' for i in range(1,10)]+['c41a','c41b','c41c']+[f'e{i}' for i in range(1,7)])
d=readzip('dataverse_files (1).zip','CLHLS2021_released.dta',columns=['ID','sex','age','year']+[f'e{i}' for i in range(1,7)])
assert c.id.is_unique and d.ID.is_unique
c=c[c.trueage.between(65,120)].merge(d.rename(columns={k:k+'_follow' for k in d if k!='ID'}),left_on='id',right_on='ID',how='left',validate='one_to_one')
c['male']=c.a1.map({1:1.,2:0.});c['school_any']=c.f1.where(c.f1.between(0,30)).gt(0).astype(float).where(c.f1.between(0,30));c['agec']=(c.trueage-75)/10;c['agec2']=c.agec**2
for k,col in [('hypertension','g15a1'),('diabetes','g15b1'),('smoke','d71'),('drink','d81')]:c[k]=c[col].map({1:1.,2:0.})
c['cvd']=union(c[['g15c1','g15d1']].apply(lambda s:s.map({1:1.,2:0.})));c['partner']=c.f41.map({1:1.,2:0.,3:0.,4:0.,5:0.})
a=c[[f'b3{i}' for i in range(1,10)]].where(lambda x:x.isin([1,2,3,4,5])).to_numpy();a=5-a;a[:,[4,6]]=4-a[:,[4,6]];c['mood']=pd.DataFrame(a,index=c.index).sum(axis=1,min_count=9)
c['memory']=c[['c41a','c41b','c41c']].where(lambda x:x.isin([0,1])).sum(axis=1,min_count=3)
for suffix,target in [('', 'start'),('_follow','end')]:
 x=c[[f'e{i}{suffix}' for i in [1,2,3,4,6,5]]].apply(lambda s:s.map({1:0.,2:1.,3:1.}));c[target]=union(x.iloc[:,:5]);c[target+'4']=union(x.iloc[:,:4]);c['six_'+target]=union(x)
c['sex_flag']=c.a1.isin([1,2])&c.sex_follow.isin([1,2])&c.a1.ne(c.sex_follow);dy=c.year_follow-c.yearin;da=c.age_follow-c.trueage;valid=c.yearin.between(2017,2019)&c.year_follow.between(2021,2022);c['age_flag']=da.lt(0)|(valid&((da-dy).abs().gt(2)|dy.le(0)));c['union_flag']=c.sex_flag|c.age_flag
frames={'HAALSI':std(h),'CLHLS':std(c)}
previous=json.loads((ROOT/'analysis/five_item_results.json').read_text());oldframes=pickle.load(open(ROOT/'analysis/five_item_frames.pkl','rb')) if (ROOT/'analysis/five_item_frames.pkl').exists() else {};out={'method':'Independent raw reconstruction; NumPy Newton logistic MLE with direct HC0 sandwich; no old analysis code executed','models':[],'raw_checks':{'HA_CESD_mismatches':cesd,'HA_status_reversal':int(bad.sum())}}
comparisons=[]
for name,q in frames.items():
 ref=next(x for x in previous['models'] if x['cohort']==name);key='prim_key' if name=='HAALSI' else 'id';old=oldframes[name].set_index(key) if name in oldframes else None;new=q.set_index(key);assert old is None or set(old.index)==set(new.index)
 errors={}
 for col in C+['mood','memory','start','end','start4','end4','six_start','six_end']+(['broad_start','broad_end'] if name=='HAALSI' else ['sex_flag','age_flag','union_flag']):
  if old is None:continue
  a=new[col].sort_index();b=old[col].sort_index();assert np.array_equal(a.isna(),b.isna()),(name,col,'missing mismatch');delta=np.max(np.abs(a.dropna().astype(float)-b.dropna().astype(float))) if a.notna().any() else 0.;errors[col]=float(delta);assert delta<1e-10,(name,col,delta)
 main=model(q[q.start.eq(0)]);f=flow(q);assert f==ref['five_flow'],(f,ref['five_flow']);four=q.copy();four['start']=four.start4;four['end']=four.end4;four=std(four);models=[('five_main',main,ref['five_main']),('original_four',model(four[four.start.eq(0)]),ref['original_four'])]
 fixed=q.copy();fixed['mood_z']=four.mood_z;fixed['memory_low_z']=four.memory_low_z
 sensitivities={'five_with_original_four_SD':model(fixed[fixed.start.eq(0)]),'age65_89':model(q[q.start.eq(0)&q.agec.between(-1,1.4)]),'all_six_independent_at_baseline_five_item_endpoint':model(q[q.six_start.eq(0)]),'cohort_specific_all_six':model(q[q.six_start.eq(0)],'six_end')}
 if name=='HAALSI':sensitivities['five_help_or_cannot']=model(q[q.broad_start.eq(0)],'broad_end')
 else:
  for flag in ['sex_flag','age_flag','union_flag']:sensitivities['exclude_'+flag]=model(q[q.start.eq(0)&~q[flag]])
 for s in ref['sensitivities']:models.append((s['analysis'],sensitivities[s['analysis']],s))
 maxdiff=0.
 for label,r,s in models:
  assert (r['n'],r['events'])==(s['n'],s['events'])
  for k in K:
   for field in ['OR','CI','p']:
    delta=float(np.max(np.abs(np.asarray(r['effects'][k][field])-np.asarray(s['effects'][k][field]))));maxdiff=max(maxdiff,delta);assert delta<1e-7,(name,label,k,field,delta)
 p=q[q.start.eq(0)].dropna(subset=C+['end']);extras={'baseline_summary':{v:{'mean':float(p[v].mean()),'SD':float(p[v].std())} for v in C+['mood','memory']}}
 if name=='HAALSI':
  extras['W2proxy_main_counts']={str(k):int(v) for k,v in p.w2rproxy.value_counts(dropna=False).items()};extras['W4proxy_main_counts']={str(k):int(v) for k,v in p.w4rproxy.value_counts(dropna=False).items()}
  words=q[[f'w2cn012s{i}' for i in range(1,11)]];rawsum=words.where(words.isin([0,1])).sum(axis=1,min_count=10);ok=rawsum.notna()&q.memory.notna();extras['delayed_words_both_complete']=int(ok.sum());extras['delayed_words_vs_derived_mismatch']=int(rawsum[ok].ne(q.memory[ok]).sum());extras['derived_score_valid_but_raw_unknown_n']=int((q.memory.notna()&~words.isin([0,1]).all(axis=1)).sum());extras['word_code_values']={str(k):int(v) for k,v in pd.Series(words.to_numpy().ravel()).value_counts(dropna=False).items()}
  extras['eating_screens']={str(w):{'screen_no_helpmissing':int((q[f'w{w}pf013'].eq(2)&q[f'w{w}pf014'].isna()).sum()),'help_missing_screen_other':int((q[f'w{w}pf013'].ne(2)&q[f'w{w}pf014'].isna()).sum())} for w in [2,4]}
 else:extras['main_flags']={v:int(p[v].sum()) for v in ['sex_flag','age_flag','union_flag']};extras['main_year_intervals']={str(k):int(v) for k,v in (p.year_follow-p.yearin).value_counts(dropna=False).items()}
 out['models'].append({'cohort':name,'flow':f,'main':main,'sensitivities':sensitivities,'rowwise_max_differences':errors,'saved_result_max_numeric_difference':maxdiff,**extras});print(name,'PASS',f,'maximum numeric difference',maxdiff,flush=True)
# Reconstruct published MI p/CI/variance using saved scalar components, independent of its code.
mi=json.loads((ROOT/'analysis/five_item_MI_results.json').read_text());out['MI_pooling_checks']=[]
for r in mi['models']:
 for k,e in r['pooled_effects'].items():
  m=mi['imputations'];B=e['between_variance'];U=e['within_variance'];T=U+(1+1/m)*B;lam=(1+1/m)*B/T;old_df=(m-1)/lam**2;dfc=r['n']-13;dfobs=(dfc+1)/(dfc+3)*dfc*(1-lam);df=1/(1/old_df+1/dfobs);beta=np.log(e['OR']);ci=np.exp(beta+np.array([-1,1])*t.ppf(.975,df)*np.sqrt(T));pv=2*t.sf(abs(beta/np.sqrt(T)),df)
  assert abs(np.sqrt(T)-e['pooled_SE_logOR'])<1e-12 and abs(df-e['Rubin_Barnard_df'])<1e-8 and np.max(abs(ci-e['CI']))<1e-10
  out['MI_pooling_checks'].append({'cohort':r['cohort'],'scenario':r['scenario'],'marker':k,'n':r['n'],'events':r['events'],'OR':e['OR'],'CI':ci.tolist(),'p':float(pv),'df':float(df),'MC_error_logOR':e['MC_error_logOR'],'MC_error_divided_by_SE':e['MC_error_logOR']/np.sqrt(T)})
(O/'independent_audit_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
# Private checkpoint for rerun of MI; never included in deliverables.
with open(O/'private_rebuilt_frames.pkl','wb') as f:pickle.dump(frames,f)
print('AUDIT PASSED. All complete-case models and saved pooling arithmetic match.',flush=True)
