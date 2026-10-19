"""Figure 1 from archived aggregate estimates or freshly reconstructed models."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter,NullLocator

ROOT=Path(__file__).resolve().parents[1]
def primary_rows(results=None):
 if results:
  obj=json.loads(Path(results).read_text());rows=[]
  for model in obj['models']:
   for marker,e in model['main']['effects'].items():
    rows.append(dict(cohort=model['cohort'],marker=marker,n=model['main']['n'],events=model['main']['events'],OR=e['OR'],CI_low=e['CI'][0],CI_high=e['CI'][1],p=e['p']))
  data=pd.DataFrame(rows)
 else:
  data=pd.read_csv(ROOT/'reference_results/audited_results_full_precision.csv')
  data=data[(data.version=='existing_published_derived_score')&(data.analysis=='five_main')].copy()
 assert len(data)==4 and not data.duplicated(['cohort','marker']).any()
 order=[('memory_low_z','HAALSI'),('memory_low_z','CLHLS'),('mood_z','HAALSI'),('mood_z','CLHLS')]
 return data.set_index(['marker','cohort']).loc[order].reset_index()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,help='Fresh audit/independent_audit_results.json; otherwise use archived aggregate results')
 ap.add_argument('--output-dir',type=Path,default=ROOT/'figures');args=ap.parse_args()
 rows=primary_rows(args.results);out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11.5,'pdf.fonttype':42,'ps.fonttype':42})
 fig=plt.figure(figsize=(9,4.1),facecolor='white');ax=fig.add_axes([.51,.22,.22,.62])
 ax.set_xscale('log');ax.set_xlim(.8,1.85);ax.set_ylim(-.5,3.5)
 ticks=[.8,1,1.25,1.5,1.8];ax.xaxis.set_major_locator(FixedLocator(ticks));ax.xaxis.set_major_formatter(FixedFormatter(['0.8','1','1.25','1.5','1.8']));ax.xaxis.set_minor_locator(NullLocator())
 ax.set_yticks([]);ax.axvline(1,color='#777777',lw=1,ls='--',zorder=1)
 for spine in ['left','right','top']:ax.spines[spine].set_visible(False)
 for x,label,ha in [(.02,'Baseline measure','left'),(.30,'Cohort (n / outcomes)','left'),(.605,'Adjusted OR','center'),(.825,'OR (95% CI)','center'),(.965,'P value','center')]:fig.text(x,.91,label,ha=ha,va='center',weight='bold',fontsize=10)
 for i,r in rows.iterrows():
  y=3-i;fy=.22+.62*(y+.5)/4
  label='Poorer delayed recall' if r.marker=='memory_low_z' else 'Greater psychological\nsymptoms'
  colour='#28556D' if r.cohort=='HAALSI' else '#6E6561';mark='s' if r.cohort=='HAALSI' else 'o'
  ax.errorbar(r.OR,y,xerr=np.array([[r.OR-r.CI_low],[r.CI_high-r.OR]]),fmt=mark,color=colour,ecolor=colour,markersize=6,capsize=3,lw=1.5,zorder=3)
  fig.text(.02,fy,label,va='center',fontsize=11.5)
  fig.text(.30,fy,f'{r.cohort} ({int(r.n)}/{int(r.events)})',va='center',fontsize=10.5)
  fig.text(.825,fy,f'{r.OR:.3f} ({r.CI_low:.3f}–{r.CI_high:.3f})',ha='center',va='center')
  fig.text(.965,fy,'<0.001' if r.p<.001 else f'{r.p:.3f}',ha='center',va='center')
 fig.add_artist(plt.Line2D([.02,.99],[.53,.53],transform=fig.transFigure,color='#DDDDDD',lw=.7))
 fig.text(.605,.065,'Adjusted OR per cohort-specific reference SD (log scale)',ha='center',fontsize=10)
 for ext in ['png','pdf','svg','tiff']:
  options={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
  fig.savefig(out/f'Figure_1_primary_forest.{ext}',dpi=600,facecolor='white',**options)
 plt.close(fig)
 rows.to_csv(out/'Figure_1_source_data.csv',index=False)
 print(f'Figure 1 written to {out}; four cohort-specific estimates; no pooled effect.')
if __name__=='__main__':main()
