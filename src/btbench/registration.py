"""Data-only forward registration and validation of complex IQ references.

The fitted map sends native/lab coordinates into the reference image. It is
inverted once to obtain the reference->lab field used by the filtering API.
Geometry and source truth are never inputs to fitting or acceptance decisions.
"""
from dataclasses import asdict, dataclass
import hashlib
import numpy as np
from scipy import linalg, ndimage, optimize


@dataclass(frozen=True)
class MotionGate:
    max_translation: float = 1.25
    max_shape_coefficient: float = .35
    max_displacement: float = 1.5
    minimum_motion_rms: float = .03
    identity_residual: float = .01
    minimum_relative_improvement: float = .15
    minimum_absolute_improvement: float = .002
    maximum_validation_error: float = .12
    maximum_condition: float = 1e4
    minimum_singular: float = .005
    minimum_fit_pixels: int = 6
    minimum_validation_pixels: int = 4
    minimum_valid_fraction: float = .95
    minimum_jacobian: float = .45
    maximum_jacobian: float = 1.8
    maximum_inverse_error: float = .05
    smooth_relative_improvement: float = .1
    smooth_absolute_improvement: float = .003
    shape_ridge: float = .01
    robust_scale: float = .03
    max_evaluations: int = 30


def _hash(a):
    a=np.ascontiguousarray(a)
    return hashlib.sha256(str(a.shape).encode()+a.dtype.str.encode()+a.tobytes()).hexdigest()


def _basis(shape):
    z,x=np.meshgrid(np.linspace(-1,1,shape[0]),np.linspace(-1,1,shape[1]),indexing='ij')
    return np.stack([np.ones(shape),z,x,np.sin(np.pi*z)*np.sin(np.pi*x)])


def linear_sample(image, coordinates, gradients=False):
    """Complex bilinear interpolation with continuous zero extension.

    Return exact within-cell derivatives when requested. Each neighbour outside
    the array is zero; crossing an edge by epsilon does not discard a whole pixel.
    No reflection, wrapping or truth-dependent support mask is used.
    """
    image=np.asarray(image);q=np.asarray(coordinates);z0=np.floor(q[0]).astype(int);x0=np.floor(q[1]).astype(int)
    dz=q[0]-z0;dx=q[1]-x0;shape=image.shape[-2:]
    def at(z,x):
        valid=(z>=0)&(z<shape[0])&(x>=0)&(x<shape[1])
        return image[:,np.clip(z,0,shape[0]-1),np.clip(x,0,shape[1]-1)]*valid[None]
    a,b,c,d=at(z0,x0),at(z0+1,x0),at(z0,x0+1),at(z0+1,x0+1)
    value=(1-dz)*(1-dx)*a+dz*(1-dx)*b+(1-dz)*dx*c+dz*dx*d
    if not gradients:return value
    return value,np.stack([(1-dx)*(b-a)+dx*(d-c),(1-dz)*(c-a)+dz*(d-b)])


def inverse_map(displacement, iterations=40):
    grid=np.indices(displacement.shape[1:],dtype=float);q=grid.copy()
    for _ in range(iterations):
        q=grid-np.stack([ndimage.map_coordinates(a,q,order=1,mode='nearest') for a in displacement])
    return q


def projected_residual(parameters, reference, points, observed, basis, ridge=None):
    """Variable-projection residual/Jacobian; nuisance complex gain per view.

    Nuisance gains are fitted only on the supplied registration pixels. They
    are used to assess geometry, never to rescale a filtered output or blood truth.
    """
    theta=parameters.reshape(2,len(basis));coordinates=points+np.einsum('dk,kp->dp',theta,basis)
    template,gradient=linear_sample(reference,coordinates,gradients=True)
    derivative=np.stack([gradient[d]*b[None] for d in range(2) for b in basis],axis=-1)
    denominator=np.maximum(np.sum(abs(template)**2,axis=1),1e-30)
    gain=np.sum(template.conj()*observed,axis=1)/denominator
    dnum=np.sum(derivative.conj()*observed[:,:,None],axis=1)
    dden=2*np.real(np.sum(template.conj()[:,:,None]*derivative,axis=1))
    dg=(dnum-gain[:,None]*dden)/denominator[:,None]
    jac=gain[:,None,None]*derivative+template[:,:,None]*dg[:,None]
    scale=np.maximum(np.sqrt(np.mean(abs(observed)**2,axis=1)),1e-12)
    residual=(gain[:,None]*template-observed)/scale[:,None]
    jac=jac/scale[:,None,None]
    residual=np.r_[residual.real.ravel(),residual.imag.ravel()]
    jac=np.concatenate([jac.real.reshape(-1,len(parameters)),jac.imag.reshape(-1,len(parameters))])
    if ridge is not None:
        weight=np.sqrt(len(residual))*np.asarray(ridge)
        residual=np.r_[residual,weight*parameters];jac=np.vstack([jac,np.diag(weight)])
    return residual,jac,gain


def _fit(reference, points, observed, basis, initial, options):
    n=len(basis);bounds=np.tile([options.max_translation]+[options.max_shape_coefficient]*(n-1),2)
    ridge=np.tile([0.]+[np.sqrt(options.shape_ridge)]*(n-1),2)
    cache={}
    def both(p):
        if 'p' not in cache or not np.array_equal(p,cache['p']):
            cache.update(p=p.copy(),value=projected_residual(p,reference,points,observed,basis,ridge))
        return cache['value']
    result=optimize.least_squares(lambda p:both(p)[0],np.clip(initial,-bounds+1e-8,bounds-1e-8),
        jac=lambda p:both(p)[1],bounds=(-bounds,bounds),loss='soft_l1',f_scale=options.robust_scale,
        max_nfev=options.max_evaluations,ftol=1e-7,xtol=1e-7,gtol=1e-7)
    _,jac,gain=projected_residual(result.x,reference,points,observed,basis)
    singular=linalg.svdvals(jac)/np.sqrt(jac.shape[0])
    return dict(parameters=result.x,gain=gain,success=bool(result.success),evaluations=result.nfev,
        cost=float(result.cost),minimum_singular=float(singular[-1]),
        condition=float(singular[0]/max(singular[-1],1e-30)))


def _validate(candidate, reference, points, observed, basis, image_shape, dense_basis, fit_points, fit_basis):
    theta=candidate['parameters'].reshape(2,len(basis))
    coordinates=points+np.einsum('dk,kp->dp',theta,basis)
    prediction=candidate['gain'][:,None]*linear_sample(reference,coordinates)
    candidate['validation_error']=float(np.linalg.norm(prediction-observed)/max(np.linalg.norm(observed),1e-30))
    full_inverse=np.einsum('dk,kzx->dzx',theta,dense_basis[:len(basis)])
    grid=np.indices(image_shape,dtype=float);forward=inverse_map(full_inverse)-grid
    candidate['field']=forward;candidate['rms']=float(np.sqrt(np.mean(forward**2)))
    candidate['maximum_displacement']=float(np.max(abs(forward)))
    zz,zx=np.gradient(forward[0]);xz,xx=np.gradient(forward[1]);det=(1+zz)*(1+xx)-zx*xz
    candidate['minimum_jacobian']=float(det.min());candidate['maximum_jacobian']=float(det.max())
    combined=np.concatenate([coordinates,fit_points+np.einsum('dk,kp->dp',theta,fit_basis)],axis=1)
    valid=(combined[0]>=0)&(combined[0]<=image_shape[0]-1)&(combined[1]>=0)&(combined[1]<=image_shape[1]-1)
    candidate['valid_fraction']=float(valid.mean())
    # The estimated inverse and its sampled forward inverse must agree.
    mapped=grid+forward
    closure=forward+np.stack([ndimage.map_coordinates(d,mapped,order=1,mode='nearest') for d in full_inverse])
    candidate['inverse_error']=float(np.sqrt(np.mean(closure**2)))
    return candidate


def _reason(candidate, baseline, options):
    if not candidate['success']:return 'solver_limit'
    if candidate['minimum_singular']<options.minimum_singular or candidate['condition']>options.maximum_condition:return 'unobservable'
    if candidate['maximum_displacement']>options.max_displacement or candidate['valid_fraction']<options.minimum_valid_fraction:return 'support'
    if (candidate['minimum_jacobian']<options.minimum_jacobian or candidate['maximum_jacobian']>options.maximum_jacobian
            or candidate['inverse_error']>options.maximum_inverse_error):return 'invalid_map'
    if candidate['rms']<options.minimum_motion_rms:return 'small_motion'
    improvement=baseline-candidate['validation_error']
    if (candidate['validation_error']>options.maximum_validation_error
            or improvement<options.minimum_absolute_improvement
            or improvement<options.minimum_relative_improvement*baseline):return 'validation'
    return 'accepted'


def register_displacement(views, reference, roi, mode='rigid', options=None, reference_views=None):
    """Fit on fixed tissue pixels and gate on disjoint pixels in each frame.

    Mode smooth uses hierarchical rigid->smooth validation. The explicit
    *_ungated modes retain valid numerical fits as diagnostic ablations.
    Fallback is exactly zero displacement. Each frame uses only its observation
    and the disjoint-time calibration reference, with no future trajectory fit.
    """
    options=options or MotionGate();gate=not mode.endswith('_ungated');model=mode.removesuffix('_ungated')
    if model not in ('none','rigid','smooth'):raise ValueError('Unknown registration mode')
    n,_,nz,nx=views.shape
    if model=='none':return np.zeros((n,2,nz,nx)),dict(registration='none',gate_enabled=False)
    if reference_views is None:
        ref=np.asarray(reference)[None].astype(complex);images=np.mean(abs(views),axis=1,keepdims=True).astype(complex)
        representation='envelope fallback'
    else:
        ref=np.asarray(reference_views,dtype=complex);images=np.asarray(views)
        if ref.shape!=images.shape[1:]:raise ValueError('Reference views must match IQ views')
        representation='complex calibrated views'
    grid=np.indices((nz,nx),dtype=float);zz,xx=np.indices((nz,nx));eligible=np.asarray(roi,bool).copy()
    eligible[[0,-1]]=False;eligible[:,[0,-1]]=False
    validation=eligible&((zz+2*xx)%3==0);training=eligible&~validation
    points=grid[:,training];check_points=grid[:,validation];dense=_basis((nz,nx))
    basis=dense[:,training];check_basis=dense[:,validation]
    labels=['reference_consistent','small_motion','accepted','validation','unobservable','solver_limit','support','invalid_map','insufficient_support','ungated']
    columns=['model','reason','baseline_validation_nrmse','rigid_validation_nrmse','smooth_validation_nrmse',
        'applied_displacement_rms','raw_displacement_rms','condition','minimum_jacobian','native_valid_fraction','solver_evaluations']
    rows=[];out=[];raw_parameters=[]
    def finite(x):return float(x) if x is not None and np.isfinite(x) else None
    for frame in images:
        chosen=None;rigid=None;smooth=None;baseline=None;reason='insufficient_support'
        if training.sum()>=options.minimum_fit_pixels and validation.sum()>=options.minimum_validation_pixels:
            y=frame[:,training];check_y=frame[:,validation]
            initial=np.zeros(2)
            _,j0,g0=projected_residual(initial,ref,points,y,basis[:1])
            baseline=float(np.linalg.norm(g0[:,None]*linear_sample(ref,check_points)-check_y)/max(np.linalg.norm(check_y),1e-30))
            sv=linalg.svdvals(j0)/np.sqrt(j0.shape[0])
            if gate and sv[-1]<options.minimum_singular:reason='unobservable'
            elif gate and baseline<options.identity_residual:reason='reference_consistent'
            else:
                candidates=[]
                for dz in (-1.,-.5,0.,.5,1.):
                    for dx in (-1.,-.5,0.,.5,1.):
                        p=np.array([dz,dx]);res,_,_=projected_residual(p,ref,points,y,basis[:1])
                        candidates.append((float(np.mean(res**2)),p))
                candidates.sort(key=lambda c:c[0]);starts=[candidates[0][1]]
                starts+=next(([p] for _,p in candidates if np.linalg.norm(p-starts[0])>=.5),[])
                fits=[_fit(ref,points,y,basis[:1],p,options) for p in starts]
                rigid=min(fits,key=lambda c:c['cost'])
                _validate(rigid,ref,check_points,check_y,check_basis[:1],(nz,nx),dense,points,basis[:1])
                raw=rigid
                reason=_reason(rigid,baseline,options)
                if reason=='accepted' or not gate:chosen=rigid
                if model=='smooth' and (not gate or chosen is None or chosen['validation_error']>=options.smooth_absolute_improvement):
                    initial=np.zeros((2,4));initial[:,0]=rigid['parameters']
                    smooth=_fit(ref,points,y,basis,initial.ravel(),options)
                    _validate(smooth,ref,check_points,check_y,check_basis,(nz,nx),dense,points,basis)
                    raw=smooth;sreason=_reason(smooth,baseline,options)
                    better=(chosen is None or (chosen['validation_error']-smooth['validation_error']>=options.smooth_absolute_improvement
                            and smooth['validation_error']<=(1-options.smooth_relative_improvement)*chosen['validation_error']))
                    if not gate or (sreason=='accepted' and better):chosen=smooth;reason=sreason
                    elif chosen is None:reason=sreason
                if not gate:reason='ungated'
                if chosen is None and raw['rms']<options.minimum_motion_rms:reason='small_motion'
        raw=smooth or rigid
        field=np.zeros((2,nz,nx)) if chosen is None else chosen['field']
        out.append(field);p=np.zeros(8)
        if raw is not None:
            if len(raw['parameters'])==2:p[[0,4]]=raw['parameters']
            else:p[:]=raw['parameters']
        raw_parameters.append(p)
        rows.append([0 if chosen is None else (1 if len(chosen['parameters'])==2 else 2),labels.index(reason),
            finite(baseline),finite(rigid['validation_error'] if rigid else None),finite(smooth['validation_error'] if smooth else None),
            float(np.sqrt(np.mean(field**2))),finite(raw['rms'] if raw else 0.),
            finite(raw['condition'] if raw else None),finite(raw['minimum_jacobian'] if raw else None),
            finite(raw['valid_fraction'] if raw else None),int(sum(c['evaluations'] for c in (rigid,smooth) if c))])
    out=np.asarray(out);counts={label:sum(row[1]==i for row in rows) for i,label in enumerate(labels)}
    return out,dict(registration=mode,representation=representation,gate_enabled=gate,
        model_coordinates='native to calibration reference; inverted for filtering',
        fit_pixels=int(training.sum()),validation_pixels=int(validation.sum()),
        fit_pixels_sha256=_hash(training),validation_pixels_sha256=_hash(validation),
        input_sha256=_hash(views),parameters_sha256=_hash(np.asarray(raw_parameters)),
        displacement_sha256=_hash(out),accepted_frames=sum(row[0]!=0 for row in rows),frames=n,
        reason_counts=counts,options=asdict(options),audit_columns=columns,audit_reason_labels=labels,audit_rows=rows)
