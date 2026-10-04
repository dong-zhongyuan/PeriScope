"""Current-input extension of the copied 20261002 benchmark. See frozen protocol."""
import argparse,copy,json,os,time,sys
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.special import logsumexp
from scipy.spatial.distance import cdist,pdist
from scipy.stats import spearmanr
from run_benchmark import A,W,O,dump,sha,mmd,measures,load_model
from ccwm import sinkhorn_divergence
from cellot_icnns import ICNN

I=A/'interim/rescreen_20261002'
BS=['cDC','classical_mono','nonclassical_mono','pDC'];RS=['astro','microglia_mhc2']
AXES=[b+'__'+r for b in BS for r in RS]
DEVICE=os.environ.get('BENCHMARK_DEVICE','cuda')
RUN_SEED=int(os.environ.get('BENCHMARK_SEED','42'))
FOCUS_AXES=AXES
RESULT_DIR=O/'results'/f'seed_{RUN_SEED}'
MODEL_DIR=O/'models' if RUN_SEED==42 else O/'models'/f'seed_{RUN_SEED}'
RESULT_DIR.mkdir(parents=True,exist_ok=True)
MODEL_DIR.mkdir(parents=True,exist_ok=True)

def sample(z,mask,n,rng,min_cells=25):
    ids=np.flatnonzero(mask);ds=[d for d in np.unique(z['donor'][ids]) if np.sum(z['donor'][ids]==d)>=min_cells]
    if not ds: raise ValueError('No adequately covered donor')
    ds=np.resize(rng.permutation(ds),n)
    return np.array([rng.choice(ids[z['donor'][ids]==d]) for d in ds])

def condition(n,c,r,ns):
    out=np.zeros((n,ns+4),np.float32);out[:,c]=1;out[:,2+r]=1;out[:,-1]=1;return out

def encode(model,x,brain=False,seed=410):
    mu=[];lv=[];u=[]
    with torch.no_grad():
        for lo in range(0,len(x),256):
            a,b=(model.enc_brain if brain else model.enc_blood)(torch.as_tensor(x[lo:lo+256],device=DEVICE))
            mu.append(a.cpu().numpy());lv.append(b.cpu().numpy())
            if not brain:u.append(model.protein_head(a).cpu().numpy())
    mu=np.concatenate(mu);lv=np.concatenate(lv)
    z=(mu+np.exp(.5*lv)*np.random.default_rng(seed).normal(size=mu.shape)).astype(np.float32)
    return z,mu,(None if brain else np.concatenate(u))

def decode(model,z):
    with torch.no_grad():return np.concatenate([model.dec_brain(torch.as_tensor(z[lo:lo+256],device=DEVICE,dtype=torch.float32)).cpu().numpy() for lo in range(0,len(z),256)])

def prepare(seed):
    O.mkdir(exist_ok=True,parents=True);(O/'cache').mkdir(exist_ok=True);(O/'models').mkdir(exist_ok=True)
    manifest=O/'prepared_eight.json'
    if manifest.exists():return
    b=dict(np.load(I/'blood.npz'));r=dict(np.load(I/'brain_train.npz'));t=dict(np.load(I/'brain_locked_test.npz'))
    states=list(b['state_vocab'].astype(str));model=load_model(seed,DEVICE);rng=np.random.default_rng(seed)
    assert np.array_equal(b['genes'],r['genes']) and np.array_equal(r['genes'],t['genes'])
    bt=~(b['is_val']|b['is_test']);rt=~r['is_val']
    assert not set(b['donor'][bt])&set(b['donor'][b['is_test']])
    assert not set(r['donor'][rt])&set(t['donor'])
    meta={'seed':seed,'n_genes':len(b['genes']),'n_proteins':model.cfg.n_proteins,'axes':AXES,'sources':{p.name:sha(p) for p in I.glob('*.npz')},'checkpoint_sha256':sha(A/f'results/optimization_20261003/models/seed_{seed}/model.pt'),'protocol_sha256':sha(W/'protocol.txt'),'blood_test_donors':b['donor_names'][np.unique(b['donor'][b['is_test']])].tolist(),'test_blood_per_condition':1}
    for axis in AXES:
        if (O/'cache'/f'{axis}.npz').exists():continue
        bs,rs=axis.split('__');bi=states.index(bs);ri=states.index(rs);cache={}
        for c in (0,1):
            for prefix,maskb,maskr,n in [('train',bt,rt,2048),('val',b['is_val'],r['is_val'],192)]:
                xb=sample(b,maskb&(b['state']==bi)&(b['cond2']==c),n,rng)
                yr=sample(r,maskr&(r['state']==ri)&(r['cond2']==c),n,rng)
                cache[f'{prefix}_x{c}']=b['X'][xb];cache[f'{prefix}_y{c}']=r['X'][yr]
                z,mu,u=encode(model,b['X'][xb],seed=seed+100*c)
                cache[f'{prefix}_z{c}']=z;cache[f'{prefix}_u{c}']=u
                cache[f'{prefix}_yz{c}']=encode(model,r['X'][yr],brain=True,seed=seed+200*c)[1]
            test_mask=b['is_test']&(b['state']==bi)&(b['cond2']==c)
            ids=sample(b,test_mask,128,rng,min_cells=1)
            cache[f'test_available_cells{c}']=np.array(int(test_mask.sum()))
            cache[f'test_unique_draws{c}']=np.array(len(np.unique(ids)))
            cache[f'test_x{c}']=b['X'][ids];z,mu,u=encode(model,b['X'][ids],seed=seed+100*c)
            cache[f'test_z{c}']=z;cache[f'test_u{c}']=u
            cache[f'test_blood_donor{c}']=b['donor'][ids]
            ds=[];obs=[];profiles=[]
            for d in np.unique(t['donor'][(t['state']==ri)&(t['cond2']==c)]):
                ii=np.flatnonzero((t['state']==ri)&(t['cond2']==c)&(t['donor']==d))
                if len(ii)<25:continue
                ds.append(d);obs.append(t['X'][rng.choice(ii,128,replace=len(ii)<128)]);profiles.append(t['X'][ii].mean(0))
            cache[f'obs{c}']=np.stack(obs);cache[f'profiles{c}']=np.stack(profiles);cache[f'donors{c}']=np.array(ds)
            ids=np.flatnonzero(rt&(r['state']==ri)&(r['cond2']==c));cache[f'mean{c}']=np.stack([r['X'][ids[r['donor'][ids]==d]].mean(0) for d in np.unique(r['donor'][ids]) if np.sum(r['donor'][ids]==d)>=25]).mean(0)
            cache[f'cond{c}']=condition(1,c,ri,len(states))[0]
        ref=np.concatenate([cache[f'train_y{c}'][:128] for c in (0,1)])
        cache['bw']=np.median(pdist(ref,'sqeuclidean'))
        cache['ot_scale']=np.array(model.cfg.ot_scale)
        cache['genes']=b['genes'];cache['blood_state']=np.array(bs);cache['brain_state']=np.array(rs)
        np.savez(O/'cache'/f'{axis}.npz',**cache)
        print('PREPARED',axis,flush=True)
    dump(manifest,meta)

def read(axis):
    k=dict(np.load(O/'cache'/f'{axis}.npz'))
    if RUN_SEED!=42:
        k.update(dict(np.load(O/'representations'/f'seed_{RUN_SEED}'/f'{axis}.npz')))
    return k

def assess(name,pred,k,seed,training=None):
    rows=[]
    for c in (0,1):
        for j,d in enumerate(k[f'donors{c}']):
            row=measures(pred[c][c],k[f'obs{c}'][j],float(k['bw']))
            row.update(condition=c,donor=int(d));rows.append(row)
    obs=k['profiles1'].mean(0)-k['profiles0'].mean(0);dirs={}
    for task,delta in {'condition_aware':pred[1][1].mean(0)-pred[0][0].mean(0),'fixed_control':pred[0][1].mean(0)-pred[0][0].mean(0),'fixed_PD':pred[1][1].mean(0)-pred[1][0].mean(0)}.items():
        rho=None if np.ptp(delta)<1e-9 else float(spearmanr(delta,obs).statistic)
        dirs[task]={'rho':rho,'delta_mse':float(np.mean((delta-obs)**2)),'predicted_rms':float(np.sqrt(np.mean(delta**2))),'observed_rms':float(np.sqrt(np.mean(obs**2))),'n_blood_PD':1,'n_blood_control':1,'n_brain_PD':len(k['donors1']),'n_brain_control':len(k['donors0'])}
    result={'method':name,'seed':seed,'blood_state':str(k['blood_state']),'brain_state':str(k['brain_state']),'distribution_by_donor':rows,'distribution_average':{key:float(np.mean([r[key] for r in rows])) for key in ['mean_mse','energy_per_sqrt_gene','mmd','variance_ratio']},'direction':dirs,'training':training}
    tag=f'{name}__{k["blood_state"]}__{k["brain_state"]}'
    dump(RESULT_DIR/f'{tag}.json',result)
    np.savez_compressed(RESULT_DIR/f'{tag}.npz',**{f'map{c}_blood{bc}':pred[c][bc] for c in (0,1) for bc in (0,1)})
    print('EVALUATED',tag,flush=True)

def full_predict(model,k,c,bc,prefix='test'):
    z=torch.as_tensor(k[f'{prefix}_z{bc}'],device=DEVICE);u=torch.as_tensor(k[f'{prefix}_u{bc}'],device=DEVICE)
    cc=torch.as_tensor(np.repeat(k[f'cond{c}'][None],len(z),0),device=DEVICE)
    with torch.enable_grad():_,q=model.solve_equilibrium(z,z.clone(),u,cc)
    return decode(model,q.detach().cpu().numpy())

def ot_plan(x,y):
    cost=cdist(x,y,'sqeuclidean');cost/=max(float(np.median(cost)),1e-8);lk=-cost/.1
    f=np.zeros(len(x));g=np.zeros(len(y));la=-np.log(len(x));lb=-np.log(len(y))
    for step in range(300):
        f=la-logsumexp(lk+g[None,:],axis=1);g=lb-logsumexp(lk+f[:,None],axis=0)
        if step%20==0:
            mat=np.exp(lk+f[:,None]+g[None,:])
            if max(abs(mat.sum(1)-1/len(x)).max(),abs(mat.sum(0)-1/len(y)).max())<1e-6:break
    mat=np.exp(lk+f[:,None]+g[None,:]);err=float(max(abs(mat.sum(1)-1/len(x)).max(),abs(mat.sum(0)-1/len(y)).max()))
    if not np.isfinite(mat).all() or err>1e-4:raise ValueError('OT marginal convergence failed')
    return mat,err

def baselines(seed,axis):
    model=load_model(seed,DEVICE);k=read(axis);preds={a:{c:{} for c in (0,1)} for a in ['shared_label_mean','shared_label_marginal','ot_ridge','ot_ridge_residual','periscope']};records=[]
    for c in (0,1):
        z=k[f'train_z{c}'];yz=k[f'train_yz{c}'];plan,err=ot_plan(z,yz)
        target=(plan@yz)/plan.sum(1)[:,None]
        feat=np.concatenate([z,k[f'train_u{c}']],axis=1).astype(float);mean=feat.mean(0);sd=np.maximum(feat.std(0),.01);f=(feat-mean)/sd;tm=target.mean(0)
        ev,V=np.linalg.eigh(f.T@f);cross=V.T@(f.T@(target-tm));scores=[];best=None
        vf=(np.concatenate([k[f'val_z{c}'],k[f'val_u{c}']],axis=1)-mean)/sd
        for alpha in [.1,1.,10.,100.,1000.]:
            coef=V@(cross/(np.maximum(ev,0)+alpha)[:,None]);val=decode(model,vf@coef+tm);score=mmd(val,k[f'val_y{c}'],float(k['bw']));scores.append({'alpha':alpha,'validation_mmd':score})
            if best is None or score<best[0]:best=(score,alpha,coef.copy())
        _,alpha,coef=best
        # Sample source index conditional on each target from training OT coupling.
        rng=np.random.default_rng(seed+10*c);indices=np.array([rng.choice(len(z),p=plan[:,j]/plan[:,j].sum()) for j in range(len(yz))]);residual=yz-(f[indices]@coef+tm);residual-=residual.mean(0)
        records.append({'condition':c,'alpha':alpha,'validation_path':scores,'ot_marginal_error':err,'pairing':'training-only OT pseudo-pairs'})
        np.savez(MODEL_DIR/f'ridge_{axis}_{c}.npz',coef=coef,mean=mean,sd=sd,target_mean=tm,residual=residual)
        eps=residual[rng.integers(len(residual),size=128)]
        for bc in (0,1):
            xx=k[f'test_x{bc}'];n=len(xx);fx=(np.concatenate([k[f'test_z{bc}'],k[f'test_u{bc}']],axis=1)-mean)/sd
            preds['shared_label_mean'][c][bc]=np.repeat(decode(model,yz).mean(0)[None],n,0)
            preds['shared_label_marginal'][c][bc]=decode(model,yz[:n])
            preds['ot_ridge'][c][bc]=decode(model,fx@coef+tm);preds['ot_ridge_residual'][c][bc]=decode(model,fx@coef+tm+eps)
            preds['periscope'][c][bc]=full_predict(model,k,c,bc)
    for name,pred in preds.items():assess(name,pred,k,seed,records if 'ridge' in name else None)


class ConditionalMap(torch.nn.Module):
    def __init__(self,kind,dim,ncontext):
        super().__init__();self.kind=kind;self.dim=dim
        if kind=='conditional_mlp':self.net=torch.nn.Sequential(torch.nn.Linear(dim+ncontext,256),torch.nn.Softplus(),torch.nn.Linear(256,256),torch.nn.Softplus(),torch.nn.Linear(256,dim))
        else:self.net=ICNN(dim+ncontext,[64]*4,kernel_init_fxn=lambda w:torch.nn.init.uniform_(w,0,.01))
    def forward(self,z,ctx):
        if self.kind=='conditional_mlp':return z+self.net(torch.cat([z,ctx],1))
        if not z.requires_grad:z=z.detach().requires_grad_(True)
        potential=self.net(torch.cat([z,ctx],1));grad=torch.autograd.grad(potential.sum(),z,create_graph=True)[0]
        return grad

def shared(seed,kind,max_steps=12000):
    model=load_model(seed,DEVICE)
    for param in model.parameters():param.requires_grad_(False)
    ks=[read(a) for a in AXES];all_z=np.concatenate([k[f'train_z{c}'] for k in ks for c in (0,1)]);all_u=np.concatenate([k[f'train_u{c}'] for k in ks for c in (0,1)])
    zm=all_z.mean(0);zs=np.maximum(all_z.std(0),.1);um=all_u.mean(0);us=np.maximum(all_u.std(0),.1)
    def features(k,c,prefix,bc=None):
        if bc is None:bc=c
        z=torch.as_tensor((k[f'{prefix}_z{bc}']-zm)/zs,device=DEVICE);u=(k[f'{prefix}_u{bc}']-um)/us
        ctx=torch.as_tensor(np.concatenate([u,np.repeat(k[f'cond{c}'][None],len(u),0)],1),device=DEVICE);return z,ctx
    zmean=torch.as_tensor(zm,device=DEVICE);zscale=torch.as_tensor(zs,device=DEVICE)
    train=[];val=[]
    for k in ks:
        for c in (0,1):
            train.append((*features(k,c,'train'),torch.as_tensor(k[f'train_y{c}'],device=DEVICE)))
            val.append((*features(k,c,'val'),k[f'val_y{c}'],float(k['bw'])))
    torch.manual_seed(seed);net=ConditionalMap(kind,model.cfg.z_dim,model.cfg.n_proteins+model.cfg.cond_dim).to(DEVICE)
    opt=torch.optim.Adam(net.parameters(),lr=2e-4);rng=np.random.default_rng(seed);best=float('inf');bad=0;beststate=None;history=[];start=time.time();reason='max_budget'
    ck=MODEL_DIR/f'{kind}_eight_axis.pt';log=MODEL_DIR/f'{kind}_eight_axis.json'
    if ck.exists() and log.exists() and json.loads(log.read_text()).get('complete'):
        net.load_state_dict(torch.load(ck,map_location=DEVICE,weights_only=False)['state_dict']);stat=json.loads(log.read_text())
    else:
        for step in range(1,max_steps+1):
            z,ctx,y=train[rng.integers(len(train))];ix=rng.integers(len(z),size=128);iy=rng.integers(len(y),size=128)
            opt.zero_grad();out=model.dec_brain(net(z[ix],ctx[ix])*zscale+zmean)
            loss=sinkhorn_divergence(out,y[iy],eps=.1,iters=20,scale=model.cfg.ot_scale);loss.backward();opt.step()
            if kind!='conditional_mlp':net.net.clamp_w()
            if step%500==0 or step==max_steps:
                scores=[]
                for vz,vc,vy,bw in val:
                    q=model.dec_brain(net(vz,vc)*zscale+zmean).detach().cpu().numpy();scores.append(mmd(q,vy,bw))
                score=float(np.mean(scores))
                if not np.isfinite(score):raise ValueError('nonfinite shared map validation')
                history.append({'step':step,'validation_mmd':score,'per_stratum':scores,'train_loss':float(loss.detach()),'elapsed_s':time.time()-start})
                if score<best-1e-5:best=score;beststep=step;beststate=copy.deepcopy(net.state_dict());bad=0
                else:bad+=1
                dump(log,{'complete':False,'history':history,'best_step':beststep});print(kind,step,score,flush=True)
                if step>=3000 and bad>=6:reason='validation_plateau';break
        net.load_state_dict(beststate);torch.save({'state_dict':beststate,'zm':zm,'zs':zs,'um':um,'us':us},ck)
        stat={'complete':True,'history':history,'best_step':beststep,'stop_reason':reason,'budget':max_steps,'training_axes':AXES,'representation':'frozen current PeriScope modules','algorithm':'conditional residual MLP' if kind=='conditional_mlp' else 'official ICNN architecture; conditioned gradient map; matched Sinkhorn loss, not original CellOT dual training'};dump(log,stat)
    for k in ks:
        pred={c:{} for c in (0,1)}
        for c in (0,1):
            for bc in (0,1):
                z,ctx=features(k,c,'test',bc);pred[c][bc]=model.dec_brain(net(z,ctx)*zscale+zmean).detach().cpu().numpy()
        assess(kind,pred,k,seed,stat)



if __name__=='__main__':
    raise SystemExit('Use shared_pipeline.py to prepare and evaluate the shared-representation benchmark.')
