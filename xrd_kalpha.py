"""Conservative Cu K-alpha doublet flagging and COD export provenance.

This module is independent of PyQt6 and never alters measured peak coordinates or
intensities. It labels candidate Kalpha2 peaks ONLY when two wavelengths originate
from the input data and match the wavelength used by the current dataset.
"""
from __future__ import annotations

import math
from typing import Any


def predicted_kalpha2_position(two_theta_kalpha1_deg: float, lambda1_A: float, lambda2_A: float) -> float | None:
    try:
        t = float(two_theta_kalpha1_deg)
        l1, l2 = float(lambda1_A), float(lambda2_A)
        if not all(math.isfinite(v) for v in (t,l1,l2)) or not (0 < t < 179 and 0 < l1 < l2):
            return None
        arg = (l2/l1)*math.sin(math.radians(t/2))
        if not (0 < arg < 1):
            return None
        return 2*math.degrees(math.asin(arg))
    except (ValueError, TypeError, OverflowError):
        return None


def verified_wavelength_pair(dataset: Any) -> tuple[float,float] | None:
    """No guessed second wavelength: require source data and unchanged primary wavelength."""
    md = getattr(dataset,'metadata',{}) or {}
    try:
        l1 = float(md['wavelength_kalpha1_A'])
        l2 = float(md['wavelength_kalpha2_A'])
        current = float(dataset.wavelength_A)
    except (KeyError,TypeError,ValueError):
        return None
    if not (0.1 < l1 < l2 < 3.0 and abs(current-l1) < 0.0001):
        return None
    return l1,l2


def find_resolved_kalpha2(peaks: list[Any], lambda1_A: float, lambda2_A: float,
                           angle_tolerance_deg: float=0.055) -> dict[int,dict]:
    """Return candidate Kalpha2 peak indices, with Kalpha1 parent and residual.

    Conservative heuristics: wavelength ratio, resolved angular separation,
    plausible secondary intensity ratio, and a one-to-one pairing. Flags are
    *possible* doublets, never confirmed phase identities. Existing peak data
    are not modified. Peak indices are list positions, not peak.index channels.
    """
    pairs=[]
    for i,p in enumerate(peaks):
        try:
            theta=float(p.two_theta); intensity=float(p.height)
            width=float(p.fwhm_deg)
        except (ValueError,TypeError,AttributeError):
            continue
        if not (math.isfinite(theta) and math.isfinite(intensity) and intensity>0):
            continue
        target=predicted_kalpha2_position(theta,lambda1_A,lambda2_A)
        if target is None:
            continue
        gap=target-theta
        if gap < 0.13 or (math.isfinite(width) and width>0 and gap < 1.05*width):
            # Unresolved/inseparable split: no independent Kalpha2 peak to exclude.
            continue
        for j,q in enumerate(peaks):
            if i==j:
                continue
            try:
                angle=float(q.two_theta); height=float(q.height)
            except (ValueError,TypeError,AttributeError):
                continue
            if not (angle>theta and math.isfinite(height) and height>0):
                continue
            ratio=height/intensity
            if not 0.12 <= ratio <= 0.90:
                continue
            residual=angle-target
            if abs(residual) <= min(angle_tolerance_deg,0.25*gap):
                pairs.append((abs(residual),i,j,residual,ratio,target))
    pairs.sort(key=lambda r:(r[0],r[1],r[2]))
    used=set(); flagged={}
    for err,i,j,residual,ratio,target in pairs:
        if i in used or j in used:
            continue
        used.add(i); used.add(j)
        flagged[j]={'parent_index':i, 'predicted_2theta':target,
                    'delta_2theta':residual, 'height_ratio':ratio, 'status':'Kalpha2 candidate'}
    return flagged


def dataset_kalpha2_candidates(dataset: Any) -> dict[int,dict]:
    pair=verified_wavelength_pair(dataset)
    if not pair:
        return {}
    return find_resolved_kalpha2(list(getattr(dataset,'peaks',[]) or []),*pair)


def cod_assignment_for_peak(dataset: Any,peak: Any) -> dict:
    """Trace only unambiguous accepted or automatically assigned COD sources."""
    md=getattr(dataset,'metadata',{}) or {}
    tt=float(peak.two_theta)
    assignments=md.get('manual_cod_peak_assignments',[]) or []
    if isinstance(assignments,list):
        best=None; best_err=math.inf
        for r in assignments:
            if not isinstance(r,dict):
                continue
            try:
                err=abs(tt-float(r['peak_two_theta']))
            except (ValueError,TypeError,KeyError):
                continue
            if err<=0.035 and err<best_err:
                best,best_err=r,err
        if best and str(getattr(peak,'phase_source','')).lower()=='cod':
            return dict(best)
    if str(getattr(peak,'phase_source','')).lower()!='cod':
        return {}
    # When automatic identification produced a unique same-named COD fit, its
    # provenance can be shown; never guess from a single unrelated COD fit.
    phase=str(getattr(peak,'phase','') or '').strip().lower()
    matches=[]
    for r in md.get('phase_fits',[]) or []:
        if not isinstance(r,dict) or r.get('cod_id') in ('',None):
            continue
        names={str(r.get(k,'') or '').strip().lower() for k in ('phase','plot_phase')}
        if phase and phase in names:
            matches.append(r)
    if len(matches)==1:
        return dict(matches[0])
    return {}


def traceability_rows(dataset: Any) -> tuple[list[str],list[list]]:
    """Build a full measured-peak ledger without requiring a GUI or an XLSX writer."""
    headers=['2theta zmierzone [deg]','d [A]','Irel [%]','Faza (opis)','hkl oryginalne',
             'Źródło przypisania','COD ID','Grupa przestrzenna','hkl COD',
             'COD 2theta [deg]','Delta 2theta [deg]','Ocena dubletu','Kalpha1 2theta [deg]',
             'Reszta dubletu [deg]','Status weryfikacji']
    pairs=dataset_kalpha2_candidates(dataset)
    rows=[]
    for idx,p in enumerate(getattr(dataset,'peaks',[]) or []):
        cod=cod_assignment_for_peak(dataset,p)
        db=pairs.get(idx,{})
        parent=(float(dataset.peaks[db['parent_index']].two_theta) if db else '')
        rows.append([float(p.two_theta),float(p.d_A),float(p.rel_intensity),
                     str(p.phase or ''),str(p.hkl or ''),str(p.phase_source or 'nieprzypisany'),
                     cod.get('cod_id',''),cod.get('space_group',''),cod.get('hkl',''),
                     cod.get('ref_two_theta',''),cod.get('delta_2theta',''),
                     'PRAWDOPODOBNY Kalpha2 — zweryfikuj' if db else '',
                     parent,db.get('delta_2theta',''),
                     'przypisanie ręczne COD' if cod and cod.get('peak_two_theta') is not None else 'wymaga oceny operatora'])
    return headers,rows
