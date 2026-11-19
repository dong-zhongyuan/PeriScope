"""Donor-balanced nested validation for frozen image embeddings."""
from dataclasses import dataclass
import numpy as np

ALPHAS=(1.,10.,100.,1000.,10000.,100000.)


def donor_weights(groups):
    groups=np.asarray(groups)
    unique,counts=np.unique(groups,return_counts=True)
    sizes=dict(zip(unique,counts))
    return np.array([len(groups)/(len(unique)*sizes[g]) for g in groups])


def constrain_composition(raw,fallback):
    result=np.maximum(raw,0)
    sums=result.sum(axis=1)
    good=sums>0
    result[good]/=sums[good,None]
    result[~good]=fallback
    return result


@dataclass
class RidgeModel:
    mean: np.ndarray
    scale: np.ndarray
    y_mean: np.ndarray
    coefficient: np.ndarray

    def predict(self,x):
        raw=((np.asarray(x,dtype=np.float64)-self.mean)/self.scale)@self.coefficient+self.y_mean
        return constrain_composition(raw,self.y_mean),raw


def fit_grid(x,y,groups,alphas=ALPHAS):
    x=np.asarray(x,dtype=np.float64);y=np.asarray(y,dtype=np.float64);groups=np.asarray(groups)
    if x.ndim!=2 or y.ndim!=2 or len(x)!=len(y) or len(groups)!=len(x):raise ValueError('Invalid calibration dimensions')
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not len(x):raise ValueError('Invalid calibration values')
    if np.min(y)<0 or not np.allclose(y.sum(1),1,atol=1e-10):raise ValueError('Targets must be RNA compositions')
    if any(a<=0 for a in alphas):raise ValueError('Ridge penalties must be positive')
    weights=donor_weights(groups);mean=np.average(x,axis=0,weights=weights)
    scale=np.sqrt(np.average((x-mean)**2,axis=0,weights=weights));scale[scale<1e-12]=1
    z=(x-mean)/scale;y_mean=np.average(y,axis=0,weights=weights)
    gram=z.T@(z*weights[:,None]);rhs=z.T@((y-y_mean)*weights[:,None])
    eig,vec=np.linalg.eigh(gram);eig=np.maximum(eig,0);rotated=vec.T@rhs
    return {float(alpha):RidgeModel(mean,scale,y_mean,vec@(rotated/(eig[:,None]+alpha))) for alpha in alphas}


def choose_alpha(x,y,groups,alphas=ALPHAS):
    groups=np.asarray(groups);donors=sorted(np.unique(groups));rows=[]
    if len(donors)<3:raise ValueError('At least three inner-fold donors required')
    for held in donors:
        train=groups!=held;test=~train
        models=fit_grid(x[train],y[train],groups[train],alphas)
        for alpha,model in models.items():
            predicted,_=model.predict(x[test])
            score=float(np.sqrt(np.mean((predicted-y[test])**2,axis=0)).mean())
            rows.append(dict(inner_test_donor=str(held),alpha=alpha,macro_state_RMSE=score,training_donors=[str(g) for g in donors if g!=held]))
    scores={float(a):float(np.mean([r['macro_state_RMSE'] for r in rows if r['alpha']==a])) for a in alphas}
    best=min(scores.values());chosen=max(a for a,s in scores.items() if abs(s-best)<=1e-12)
    return chosen,rows,scores


def outer_fold(x,y,groups,held,alphas=ALPHAS):
    x=np.asarray(x);y=np.asarray(y);groups=np.asarray(groups)
    train=groups!=held;test=groups==held
    if not np.any(test) or len(np.unique(groups[train]))<3:raise ValueError('Invalid outer donor holdout')
    # No held-out target is passed to preprocessing or parameter selection.
    alpha,inner,scores=choose_alpha(x[train],y[train],groups[train],alphas)
    model=fit_grid(x[train],y[train],groups[train],[alpha])[alpha]
    predictions,raw=model.predict(x[test])
    return dict(test_indices=np.flatnonzero(test),prediction=predictions,raw_prediction=raw,baseline=np.tile(model.y_mean,(test.sum(),1)),selected_alpha=alpha,inner=inner,inner_scores=scores,training_donors=[str(g) for g in sorted(np.unique(groups[train]))],test_donor=str(held),model=model)
