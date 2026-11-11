"""Compare C++ CSV reader/evaluator to original validated binary inputs."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys
import numpy as np


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--tables',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'converters'/'dirac'))
    from cc_binary import read_pack
    from cc_physics import delta_tail,stopping_at_energy,interpolate_B,bracket
    with args.probe.open(newline='') as f:
        rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    cache={}
    worst=0.;worst_sample=0.;sample_count=0
    for r in rows:
        z=int(r['Z']);zt=int(r['Zt']);n=int(r['N']);e=r['E'];cut=r['cut']
        if z not in cache:
            _,d=read_pack(args.tables/'delta'/f'Z{z:03d}.ccbin')
            _,s=read_pack(args.tables/'stopping'/f'Z{z:03d}.ccbin')
            cache[z]=(d,s)
        d,s=cache[z]
        ti=int(np.flatnonzero(s['target_Z']==zt)[0])
        m=delta_tail(d,e,n,cut,z)
        values={'M0':m[0],'M1':m[1],'M2':m[2],
                'full':stopping_at_energy(d,s,[e],z)[ti,0,n],
                'B':interpolate_B(s,[e])[ti,0]}
        sl=slice(int(s['offsets'][ti]),int(s['offsets'][ti+1]))
        grid=s['energy_MeV_u'][sl]
        _,lo,hi,w=bracket(grid,[e])
        eq=(1-w[0])*s['atima_MeV_cm2'][sl][lo[0]]+w[0]*s['atima_left_MeV_cm2'][sl][hi[0]]
        if e==grid[-1]:eq=s['atima_MeV_cm2'][sl][-1]
        values['eq']=eq
        for k,expected in values.items():
            rel=abs(r[k]-expected)/max(abs(expected),1e-300)
            worst=max(worst,rel)
            if rel>3e-10:
                raise ValueError(f'{k} mismatch Z={z} N={n} E={e} cut={cut}: {r[k]} vs {expected}, rel={rel}')
        if r['samples']:
            mean=m[1]/m[0]
            variance=max(m[2]/m[0]-mean*mean,0)
            sem=math.sqrt(variance/r['samples'])
            score=abs(r['mean_sample']-mean)/max(sem,1e-300)
            worst_sample=max(worst_sample,score);sample_count+=int(r['samples'])
            if score>6:
                raise ValueError(f'Sampler mean discrepancy {score} standard errors')
    result={'passed':True,'rows':len(rows),'maximum_relative_error':worst,
            'sampled_deltas':sample_count,'maximum_sample_mean_standard_errors':worst_sample}
    args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))


if __name__=='__main__':main()
