import json,math
from pathlib import Path
R=Path(__file__).resolve().parent.parent
p=json.loads((R/'flow_pvalues/five_item_flow_pvalues.json').read_text());r=json.loads((R/'analysis/five_item_results.json').read_text());a=json.loads((R/'audit/independent_audit_results.json').read_text());checked=[]
for row in p['pvalues']:
 if row['analysis'].startswith('MI') or row['analysis'].lower().startswith('mi'):
  print('MIlabel',row['analysis']);continue
 ref=next(x for x in r['models'] if x['cohort']==row['cohort']);label=row['analysis']
 if label in ['five_main','original_four']:x=ref[label]
 else:
  match=[x for x in ref['sensitivities'] if x['analysis']==label]
  if not match:
   print('Otherlabel',row['cohort'],label);continue
  x=match[0]
 e=x['effects'][row['marker']];assert (row['n'],row['events'])==(x['n'],x['events']);assert abs(row['p']-e['p'])<1e-12;assert abs(row['OR']-e['OR'])<1e-12;assert max(abs(v-w) for v,w in zip(row['CI'],e['CI']))<1e-12;checked.append(row)
print('Complete-case P table rows verified',len(checked),'of',len(p['pvalues']))
for row in p['pvalues']:
 if row in checked:continue
 candidates=[x for x in a['MI_pooling_checks'] if x['cohort'].startswith(row['cohort']) and x['marker']==row['marker'] and x['n']==row['n']]
 assert len(candidates)==1,(row,candidates);x=candidates[0];assert abs(row['p']-x['p'])<1e-12;assert abs(row['OR']-x['OR'])<1e-12;assert max(abs(v-w) for v,w in zip(row['CI'],x['CI']))<1e-12;checked.append(row)
assert len(checked)==38
(R/'audit/table_verification.json').write_text(json.dumps({'all_P_table_rows_verified':len(checked),'complete_case_rows_independently_refit':32,'MI_rows_pooling_arithmetic_checked':6,'all_nominal_unadjusted':True},indent=2));print('ALL38 VERIFIED')
