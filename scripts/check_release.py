"""Check saved main estimates and the independent HC0 solver without raw data."""
from pathlib import Path
import json,ast
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
source=(ROOT/'audit/independent_audit.py').read_text();tree=ast.parse(source)
# Load function definitions only: never import the raw-data audit for this check.
definitions=[x for x in tree.body if isinstance(x,(ast.FunctionDef,ast.Import,ast.ImportFrom))]
namespace={'__file__':str(ROOT/'audit/independent_audit.py')}
exec(compile(ast.Module(body=definitions,type_ignores=[]),'independent_solver_functions','exec'),namespace)
y=np.array([0.]*70+[1.]*30);X=np.ones((100,1));b,V,_,_=namespace['logistic'](X,y)
assert np.allclose(b,[np.log(.3/.7)],atol=1e-10)
assert np.allclose(V,[[1/(100*.3*.7)]],atol=1e-10)
u=namespace['union'](pd.DataFrame([[0,0],[0,np.nan],[1,np.nan],[np.nan,np.nan]]))
assert u.iloc[0]==0 and np.isnan(u.iloc[1]) and u.iloc[2]==1 and np.isnan(u.iloc[3])
data=pd.read_csv(ROOT/'reference_results/audited_results_full_precision.csv');d=data[(data.version=='existing_published_derived_score')&(data.analysis=='five_main')]
assert len(d)==4
saved=json.loads((ROOT/'audit/independent_audit_results.json').read_text())
for _,r in d.iterrows():
 m=next(m['main'] for m in saved['models'] if m['cohort']==r.cohort);e=m['effects'][r.marker]
 assert (r.n,r.events)==(m['n'],m['events'])
 assert np.allclose([r.OR,r.CI_low,r.CI_high,r.p],[e['OR'],*e['CI'],e['p']],atol=1e-10,rtol=0)
print('PASS: four primary effect rows, analytic intercept-only HC0 check, missing/positive outcome logic. No raw-data models rerun.')
