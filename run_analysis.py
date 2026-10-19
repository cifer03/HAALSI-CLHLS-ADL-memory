"""Command-line entry point. Input-free checks never claim a raw-data rerun."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,importlib.util
ROOT=Path(__file__).resolve().parent
def run(relative,*args):
 env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
 subprocess.run([sys.executable,str(ROOT/relative),*map(str,args)],cwd=ROOT,env=env,check=True)
def check_inputs(strict=False):
 base=Path(os.environ.get('HAALSI_CLHLS_INPUTS',ROOT/'inputs'));issues=[]
 for r in json.loads((ROOT/'inputs/input_manifest.json').read_text()):
  p=base/r['filename']
  if not p.is_file():issues.append(f'Missing: {p}');continue
  with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
  if digest!=r['sha256']:
   msg=f'Hash differs from archived input: {p.name}. Re-downloaded ZIP containers may differ; verify release and archive members.'
   if strict:issues.append(msg)
   else:print('WARNING:',msg)
 if issues:raise SystemExit('\n'.join(issues))
 print('All four input files are present.')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['check-inputs','verify-archive','main','mi-original','mi-audit','figure','all']);ap.add_argument('--strict-hashes',action='store_true');a=ap.parse_args()
 if a.stage=='check-inputs':check_inputs(a.strict_hashes);return
 if a.stage=='verify-archive':run('audit/verify_tables.py');run('scripts/check_release.py');return
 if a.stage=='figure':run('scripts/plot_forest.py');return
 check_inputs(a.strict_hashes)
 if a.stage in ['mi-original','all'] and importlib.util.find_spec('statsmodels') is None:raise SystemExit('Install requirements-original-mi.txt to rerun the original seed-20261011 MI implementation.')
 # Raw reconstruction is also imported by the MI audit scripts. It writes a private
 # participant-level checkpoint locally; .gitignore excludes this file.
 run('audit/independent_audit.py')
 run('scripts/plot_forest.py','--results',ROOT/'audit/independent_audit_results.json')
 if a.stage in ['mi-original','all']:
  shutil.copy2(ROOT/'audit/private_rebuilt_frames.pkl',ROOT/'analysis/five_item_frames.pkl')
  run('analysis/five_item_MI.py')
 if a.stage in ['mi-audit','all']:
  run('audit/score_quality_and_MI.py');run('audit/rerun_original_HA_MI.py')
 print('Requested raw-data stage finished. Inspect numerical comparisons before replacing archived results.')
if __name__=='__main__':main()
