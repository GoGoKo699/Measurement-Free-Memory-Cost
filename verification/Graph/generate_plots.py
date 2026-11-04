#!/usr/bin/env python3
"""Reproduce finite memory-resource curves from analytical certificates.

Run `python verification/Graph/generate_plots.py --check` from the repository root.
Use --output-dir and --data-dir to keep regenerated files at chosen locations.
All counts use k=1, gamma=epsilon=1/100 and D=1..1000. Dates are omitted from
the vector PDF metadata to make the exported asset reproducible.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.metadata, json, math, platform, tempfile
from fractions import Fraction
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, NullLocator
import numpy as np

SOURCE_COMMIT='3bc227af0dcb359c211dee7caf13e223f7358873'
GAMMA=Fraction(1,100)
EPSILON=Fraction(1,100)
K=1
DMAX=1000

def h2(x): return -x*math.log2(x)-(1-x)*math.log2(1-x)
def dbin(a,p): return a*math.log(a/p)+(1-a)*math.log((1-a)/(1-p))
def floor_fraction(x): return x.numerator//x.denominator
def odd_ceiling(x):
    n=max(1,math.ceil(x))
    return n if n%2 else n+1

def coherent(k,gamma,eps,D):
    tau=floor_fraction(1/(20*gamma)); J=(D+tau-1)//tau
    A=max(1.,math.log(k*J/(10*float(eps)))/math.log(8/3))
    L=math.ceil(math.log2(A)); n=k*5**L
    T=42*k*J*(5**L-1)//4
    err=k*J*math.exp(2**L*math.log(3/8))/10
    assert tau>=1/(40*float(gamma))-1e-12
    assert J<=40*float(gamma)*D+1+1e-12
    assert err<=float(eps)*(1+1e-12)
    return dict(coherent_tau=tau,coherent_J=J,coherent_A=A,coherent_L=L,
                coherent_n=n,coherent_T_upper=T,coherent_error_certificate=err)

def archive(k,gamma,eps,D):
    tau=floor_fraction(1/(16*gamma)); J=(D+tau-1)//tau
    Gamma=-math.expm1(tau*math.log1p(-float(gamma))); Q=1-Gamma
    R0=1-h2(1/8)-math.log2(3)/8; c=dbin(1/16,3*Gamma/4)
    n0=math.ceil(max((k+1)/R0,math.log(2*J/float(eps))/c)); s=n0-k
    dataerr=J*math.exp(-n0*c)
    assert Gamma<=1/16+1e-15 and 3*Gamma/4<1/16 and s>0
    assert R0*n0>=k+1-1e-12 and dataerr<=float(eps)/2*(1+1e-12)
    if J==1:
        B=0; slots=0; recorderr=0.
    else:
        M=s*(J-1); requirement=8/Q**2*math.log(2*M/float(eps))
        B=odd_ceiling(requirement); slots=B*s*(2*3**(J-2)-1)
        recorderr=M*math.exp(-B*Q**2/8)
        assert B%2==1 and B>=requirement and B-2<requirement
        assert recorderr<=float(eps)/2*(1+1e-12)
        # Exact peak includes retained parents and new targets simultaneously.
        assert slots==B*s*(3**(J-2)-1)//2+B*s*(3**(J-1)-1)//2
    N=n0+slots
    assert dataerr+recorderr<=float(eps)*(1+1e-12)
    return dict(archive_tau=tau,archive_J=J,archive_Gamma=Gamma,archive_Q=Q,
                archive_R0=R0,archive_c=c,archive_n0=n0,archive_s=s,archive_B=B,
                archive_N_decimal=str(N),archive_log10_N=math.log10(N),
                archive_data_error_certificate=dataerr,
                archive_record_error_certificate=recorderr,
                archive_error_certificate=dataerr+recorderr,archive_T=0)

def points():
    eps=float(EPSILON); gamma=float(GAMMA); k=K
    kappa=-math.log1p(-gamma); a=2*k*(1-eps-4**(-k))/(1-4**(-k))
    rows=[]
    for D in range(1,DMAX+1):
        co=coherent(k,GAMMA,EPSILON,D); ar=archive(k,GAMMA,EPSILON,D)
        reserve=math.log(2*co['coherent_n']/a)
        lower=a*max(0.,D*kappa-reserve)
        widthlower=max(k,a*math.exp(D*kappa)/2)
        assert lower<=co['coherent_T_upper']+1e-10
        assert 2*co['coherent_n']*(1-gamma)>=a
        assert int(ar['archive_N_decimal'])>=widthlower-1e-10
        rows.append(dict(D=D,k=k,gamma=gamma,epsilon=eps,gamma_D=gamma*D,
                         a=a,kappa=kappa,logarithmic_reserve=reserve,
                         zero_T_width_necessary=widthlower,
                         coherent_T_log_lower=lower,
                         coherent_T_log_lower_integer=math.ceil(lower),**co,**ar))
    r=rows[-1]
    assert r['D']==1000 and r['coherent_n']==125 and r['coherent_T_upper']==260400
    assert r['coherent_L']==3 and r['coherent_J']==200
    assert abs(r['coherent_error_certificate']-.007821321487426758)<1e-15
    return rows

def versions():
    names=('matplotlib','numpy','contourpy','cycler','fonttools','kiwisolver',
           'packaging','pillow','pyparsing','python-dateutil','six')
    return {n:importlib.metadata.version(n) for n in names}

def plot(rows,path):
    plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],
        'mathtext.fontset':'cm','font.size':8.3,'axes.labelsize':9,'axes.titlesize':9,
        'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7.15,
        'legend.frameon':False,'axes.linewidth':.65,'xtick.major.width':.65,
        'ytick.major.width':.65,'lines.linewidth':1.5,'pdf.fonttype':42,
        'ps.fonttype':42,'axes.unicode_minus':False})
    blue,green,rose='#4477AA','#228833','#CC6677'
    x=np.array([r['gamma_D'] for r in rows])
    n=np.array([r['coherent_n'] for r in rows])
    lower=np.array([r['coherent_T_log_lower'] for r in rows])
    upper=np.array([r['coherent_T_upper'] for r in rows])
    fig,axs=plt.subplots(1,2,figsize=(7.12,3.55))
    fig.subplots_adjust(left=.075,right=.987,bottom=.17,top=.90,wspace=.27)
    ax=axs[0]
    ax.plot(x,[r['archive_log10_N'] for r in rows],color=rose,
            label='Archive construction (sufficient)')
    ax.step(x,np.log10(n),where='post',color=green,
            label='Coherent recovery (uses T gates)')
    ax.plot(x,np.log10([r['zero_T_width_necessary'] for r in rows]),color=blue,
            linestyle=(0,(4,2)),label='Zero-T necessary')
    ax.set(xlim=(0,10),ylim=(-1.6,90),xlabel=r'Noise exposure, $\gamma D$',
           ylabel=r'Logarithmic width, $\log_{10} n$')
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.set_title('(a) Physical width',loc='left',pad=8)
    ax.legend(loc='upper left',handlelength=2.1,borderaxespad=.35,labelspacing=.65)
    # The archive certificate spans many decades. Resolve the two small-width
    # curves in the otherwise empty lower-right area, without adding new data.
    detail=ax.inset_axes([.51,.12,.46,.31])
    detail.step(x,np.log10(n),where='post',color=green,linewidth=1.2)
    detail.plot(x,np.log10([r['zero_T_width_necessary'] for r in rows]),
                color=blue,linestyle=(0,(4,2)),linewidth=1.2)
    detail.set(xlim=(0,10),ylim=(0,5),xticks=[0,5,10],yticks=[0,2.5,5])
    detail.set_title('Width detail',fontsize=7.1,pad=3)
    # Coordinates and units are identical to the parent panel.
    detail.tick_params(labelsize=6.5,direction='out',length=2,pad=1.5,width=.5)
    detail.grid(axis='y',color='#e0e0e0',linewidth=.35)
    for spine in detail.spines.values():spine.set_linewidth(.5)
    ax=axs[1]
    first=int(np.flatnonzero(lower>0)[0]); zero_end=x[first-1]
    ax.axvspan(0,zero_end,facecolor='#f0f0f0',edgecolor='none',zorder=-3)
    ax.step(x,upper,where='post',color=green,label='Coherent recovery (sufficient)')
    ax.plot(x,np.ma.masked_where(lower<=0,lower),color=blue,linestyle=(0,(4,2)),
            label='Processing lower bound')
    ax.set_yscale('log'); ax.yaxis.set_minor_locator(NullLocator())
    ax.set(xlim=(0,10),ylim=(.01,1e6),xlabel=r'Noise exposure, $\gamma D$',
           ylabel=r'Internal T/T$^{\dagger}$ count, $m_T$')
    ax.set_title('(b) Processing at the coherent width',loc='left',pad=8)
    ax.text(.04,.10,'Width allowance:\nlogarithmic bound = 0',transform=ax.transAxes,
            fontsize=7.2,va='bottom',color='#555555',linespacing=1.35)
    ax.legend(loc='center',bbox_to_anchor=(.52,.50),handlelength=2.1,
              borderaxespad=.35,labelspacing=.65)
    ax.scatter([10],[260400],color=green,s=17,zorder=5,clip_on=False)
    for ax in axs:
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.tick_params(direction='out',length=3.2,pad=3)
        ax.grid(axis='y',which='major',color='#d9d9d9',linewidth=.45,alpha=.7)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,metadata={'Title':'Finite resource certificates for noisy quantum memory',
        'Author':None,'Creator':'Matplotlib','CreationDate':None,'ModDate':None})
    plt.close(fig)

def generate(output_dir, data_dir):
    rows=points(); data_dir.mkdir(parents=True,exist_ok=True)
    csvpath=data_dir/'resource_comparison.csv'
    with csvpath.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    pdfpath=output_dir/'resource_comparison.pdf'; plot(rows,pdfpath)
    meta={'kind':'Finite analytical certificates, no empirical data or fitted curves',
      'scientific_source_commit':SOURCE_COMMIT,
      'parameters':{'k':K,'gamma':str(GAMMA),'epsilon':str(EPSILON),'D_range':[1,DMAX]},
      'curves':{'zero_T_necessary_width':'max(k,a*exp(D*kappa)/2)',
        'processing_lower_at_coherent_width':'a*max(0,D*kappa-ln(2*n_coherent(D)/a))',
        'coherent':'Finite integer prescription in coherent_recovery.md, section4',
        'archive':'Finite integer prescription in clifford_lifetime.md, equations2-7; complete parent/child overlap'},
      'checked_points':len(rows),'all_error_certificates_le_target':True,
      'all_lower_bounds_le_sufficient_costs':True,'archive_overlap_checked':True,
      'first_D_positive_T_lower':next(r['D'] for r in rows if r['coherent_T_log_lower']>0),
      'D1000_example':rows[-1],
      'limits':['Lower bounds and sufficient construction costs do not give an exact frontier.',
        'The zero-T width bound does not constrain the T-using coherent construction.',
        'Both constructions satisfy half-diamond error<=epsilon, implying the converse Bell-fidelity premise.',
        'All control between noise layers is ideal/global; no hardware-gate-time conclusion follows.',
        'Zero processing lower bounds remain in the CSV and are masked on logarithmic axes.'],
      'python':platform.python_version(),'dependencies':versions(),
      'pdf_sha256':hashlib.sha256(pdfpath.read_bytes()).hexdigest(),
      'csv_sha256':hashlib.sha256(csvpath.read_bytes()).hexdigest(),
      'pdf_dates':'omitted'}
    (data_dir/'resource_comparison.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta


def guard_output(directory, expected_dir, parser):
    """Keep regeneration from modifying the committed reference directory."""
    resolved=directory.resolve()
    if resolved==expected_dir or resolved.is_relative_to(expected_dir):
        parser.error('Output directories must lie outside the committed expected/ references.')


def check_outputs(output_dir, data_dir, expected_dir, strict_provenance=False):
    """Compare newly generated outputs with independent committed references."""
    manifest=json.loads((expected_dir/'checksums.json').read_text())
    generated={'resource_comparison.pdf':output_dir/'resource_comparison.pdf',
               'resource_comparison.csv':data_dir/'resource_comparison.csv',
               'resource_comparison.json':data_dir/'resource_comparison.json'}
    result={}; provenance={}
    for name,path in generated.items():
        reference=expected_dir/name
        expected_bytes=reference.read_bytes()
        reference_sha=hashlib.sha256(expected_bytes).hexdigest()
        declared=manifest['files'][name]
        if reference_sha!=declared['sha256'] or len(expected_bytes)!=declared['size_bytes']:
            raise RuntimeError('Committed reference fails its inventory: '+name)
        actual=path.read_bytes()
        actual_sha=hashlib.sha256(actual).hexdigest()
        same_bytes=actual==expected_bytes
        if name.endswith('.json'):
            expected_values=json.loads(expected_bytes)
            actual_values=json.loads(actual)
            # Runtime labels are provenance rather than analytical values.
            # Report every difference while checking all other fields exactly.
            provenance={'expected':{k:expected_values.pop(k,None) for k in ('python','dependencies')},
                        'generated':{k:actual_values.pop(k,None) for k in ('python','dependencies')}}
            provenance['matches']=provenance['expected']==provenance['generated']
            if actual_values!=expected_values:
                raise RuntimeError('Generated numerical certificate differs: '+name)
        elif not same_bytes:
            raise RuntimeError('Generated output differs from committed reference: '+name)
        result[name]={'reference_sha256':reference_sha,
                      'generated_sha256':actual_sha,
                      'bytes_identical':same_bytes}
    if strict_provenance and not provenance['matches']:
        raise RuntimeError('Analytical outputs match, but strict runtime provenance differs: '+
                           json.dumps(provenance,sort_keys=True))
    return {'status':'PASS','checked_points':DMAX,
            'parameters':{'k':K,'gamma':str(GAMMA),'epsilon':str(EPSILON),
                          'D_range':[1,DMAX]},
            'reference_comparison':result,
            'certificate_fields_match':True,
            'runtime_provenance':provenance,
            'strict_provenance_requested':strict_provenance,
            'warnings':([] if provenance['matches'] else
                        ['Runtime provenance differs; PDF, CSV and analytical certificate fields match.']),
            'scope':'Finite analytical certificates and output reproduction; see the written proofs.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    here=Path(__file__).resolve().parent
    parser.add_argument('--check',action='store_true',
                        help='Regenerate and compare with committed references; temporary output by default.')
    parser.add_argument('--strict-provenance',action='store_true',
                        help='With --check, also require exact Python/dependency-version provenance.')
    parser.add_argument('--output-dir',type=Path,
                        help='Directory for the regenerated PDF (default: generated/Graph beside this script).')
    parser.add_argument('--data-dir',type=Path,
                        help='Directory for regenerated CSV/JSON (default: generated/data beside this script).')
    args=parser.parse_args()
    if args.strict_provenance and not args.check:
        parser.error('--strict-provenance requires --check.')
    expected_dir=(here/'expected').resolve()
    if args.output_dir is not None:guard_output(args.output_dir,expected_dir,parser)
    if args.data_dir is not None:guard_output(args.data_dir,expected_dir,parser)
    if args.check:
        with tempfile.TemporaryDirectory(prefix='memory-cost-graph-') as temporary:
            base=Path(temporary)
            output_dir=args.output_dir if args.output_dir is not None else base/'Graph'
            data_dir=args.data_dir if args.data_dir is not None else base/'data'
            generate(output_dir,data_dir)
            print(json.dumps(check_outputs(output_dir,data_dir,expected_dir,args.strict_provenance),indent=2))
    else:
        output_dir=args.output_dir if args.output_dir is not None else here/'generated/Graph'
        data_dir=args.data_dir if args.data_dir is not None else here/'generated/data'
        guard_output(output_dir,expected_dir,parser);guard_output(data_dir,expected_dir,parser)
        meta=generate(output_dir,data_dir)
        print(json.dumps({'pdf':str(output_dir/'resource_comparison.pdf'),
                          'pdf_sha256':meta['pdf_sha256'],'D1000':meta['D1000_example']},indent=2))


if __name__=='__main__': main()

