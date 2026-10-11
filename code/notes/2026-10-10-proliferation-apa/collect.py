from pathlib import Path
import hashlib,json,subprocess
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def main():
 config=json.loads((ROOT/'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json').read_text())
 code='CONFIG = '+repr(config)+'\nPROTOCOL_SHA = '+repr(hashlib.sha256((HERE/'PROTOCOL.zh.md').read_bytes()).hexdigest())+'\n'+(HERE/'extract_remote.py').read_text()
 with (HERE/'inputs.partial.npz').open('wb') as out,(HERE/'extract.log').open('wb') as log:
  subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','213docker','OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/scanpy/bin/python -'],input=code.encode(),stdout=out,stderr=log,check=True)
 (HERE/'inputs.partial.npz').rename(HERE/'inputs.npz');print((HERE/'extract.log').read_text())
if __name__=='__main__':main()
