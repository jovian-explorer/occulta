"""Data sources: where each mission's light curves live, which stars a unit holds, and
how star identifiers are encoded.

Star ids are integers in one namespace: TESS stars keep their TIC number, K2 stars are
2e10 + EPIC, Kepler stars are 3e10 + KIC. Units (one sector, campaign or quarter) are
integer codes: QLP sector s -> s, SPOC sector s -> 1000 + s, K2 campaign c -> 2000 + c,
Kepler quarter file q -> 3000 + q, TESS-SPOC full-frame light curves of sector s -> 4000 + s,
ZTF light curves of filter zg, zr, zi -> 5001, 5002, 5003 (5000 + filter id; one unit per filter because zero points differ). All times are stored as BTJD (BJD - 2457000).
"""
import io, re
import numpy as np
import requests

K2_OFF, KEP_OFF = 20_000_000_000, 30_000_000_000
KEPLER_TO_BTJD = 2454833.0 - 2457000.0

QLP_URL = ('https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HLSP/qlp/s{s:04d}/{a}/{b}/{c}/{d}/'
           'hlsp_qlp_tess_ffi_s{s:04d}-{tic:016d}_tess_{v}_llc.fits')
QLP_TARGETS = 'https://archive.stsci.edu/hlsps/qlp/target_lists/s{s:04d}.csv'
SPOC_SCRIPT = 'https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{s}_lc.sh'
K2SFF_SCRIPT = 'https://archive.stsci.edu/hlsps/k2sff/download_scripts/hlsp_k2sff_k2_lightcurve_c{c:02d}_kepler_v1_bulk-dl-llc.sh'
TESS_SPOC_URL = ('https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HLSP/tess-spoc/s{s:04d}/target/{a}/{b}/{c}/{d}/'
                 'hlsp_tess-spoc_tess_phot_{tic:016d}-s{s:04d}_tess_v1_lc.fits')
TESS_SPOC_TARGETS = 'https://archive.stsci.edu/hlsps/tess-spoc/target_lists/s{s:04d}.csv'
KEPLER_URL = 'https://archive.stsci.edu/missions/kepler/lightcurves/{k4}/{k9}/kplr{k9}-{ts}_llc.fits'
# long-cadence file stamps, Q0 .. Q17 (Q4 is split in two)
KEPLER_STAMPS = ['2009131105131', '2009166043257', '2009259160929', '2009350155506', '2010009091648',
                 '2010078095331', '2010174085026', '2010265121752', '2010355172524', '2011073133259',
                 '2011177032512', '2011271113734', '2012004120508', '2012088054726', '2012179063303',
                 '2012277125453', '2013011073258', '2013098041711', '2013131215648']
K2_CAMPAIGNS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 102, 111, 112, 12, 13, 14, 15, 16, 17, 18, 19, 91, 92]
QLP_SECTORS = list(range(104, 0, -1))
ZTF_FILTERS = {'zg': 5001, 'zr': 5002, 'zi': 5003}    # unit = 5000 + filter id; ZTF is searched by position, not from a target list


def mission_of(sid):
    return 'Kepler' if sid >= KEP_OFF else 'K2' if sid >= K2_OFF else 'TESS'


def name_of(sid):
    m = mission_of(sid)
    return f'KIC {sid - KEP_OFF}' if m == 'Kepler' else f'EPIC {sid - K2_OFF}' if m == 'K2' else f'TIC {sid}'


def unit_name(u):
    if u >= 5000: return 'ZTF ' + {v: k for k, v in ZTF_FILTERS.items()}.get(u, f'filter {u - 5000}')
    if u >= 4000: return f'TESS S{u - 4000} FFI SPOC'
    if u >= 3000: return f'Kepler file {u - 3000}'
    if u >= 2000: return f'K2 C{u - 2000}'
    if u >= 1000: return f'TESS S{u - 1000} 2-min'
    return f'TESS S{u} FFI'


def fmt_of_unit(u):
    return 'ztf' if u >= 5000 else 'spoc' if u >= 4000 else 'kepler' if u >= 3000 else 'k2sff' if u >= 2000 else 'spoc' if u >= 1000 else 'qlp'


def units(mission):
    if mission == 'QLP':
        return QLP_SECTORS
    if mission == 'SPOC':
        return [1000 + s for s in QLP_SECTORS]
    if mission == 'TESS-SPOC':
        return [4000 + s for s in QLP_SECTORS]
    if mission == 'ZTF':
        return list(ZTF_FILTERS.values())
    if mission == 'K2':
        return [2000 + c for c in K2_CAMPAIGNS]
    if mission == 'Kepler':
        return [3000 + i for i in range(len(KEPLER_STAMPS))]
    raise ValueError(mission)


def _script_urls(url, pattern):
    r = requests.get(url, timeout=900)
    if r.status_code != 200:
        return None
    out = {}
    for line in r.text.splitlines():
        m = re.search(pattern, line)
        u = re.findall(r'https://[^\'" ]+', line)
        if m and u:
            out[int(m.group(1))] = u[-1]
    return out


def unit_targets(u, kepler_ids=None):
    """{star id: url or None} for every star with a light curve in this unit (None = build from template)."""
    if u >= 5000:
        return None        # ZTF has no per-unit target list (IRSA asks for HATS above ~1000 light curves): use `occulta ztf`
    if u < 1000 or u >= 4000:
        r = requests.get((QLP_TARGETS if u < 1000 else TESS_SPOC_TARGETS).format(s=u % 1000 if u >= 4000 else u), timeout=900)
        if r.status_code != 200:
            return None
        ids = np.loadtxt(io.StringIO(r.text), delimiter=',', skiprows=1, usecols=0, dtype=np.int64, ndmin=1)
        return {int(i): None for i in ids}
    if 1000 <= u < 2000:
        d = _script_urls(SPOC_SCRIPT.format(s=u - 1000), r'-s\d{4}-(\d{16})-')
        return None if d is None else {int(k): v for k, v in d.items()}
    if 2000 <= u < 3000:
        d = _script_urls(K2SFF_SCRIPT.format(c=u - 2000), r'lightcurve_(\d+)-c\d+_kepler')
        return None if d is None else {K2_OFF + k: v for k, v in d.items()}
    # Kepler: every catalogued target is tried; quarters it was not observed in return 404
    return {KEP_OFF + k: None for k in kepler_ids}


def url_for(u, sid, versions=('v01', 'v02')):
    if u < 1000:
        z = f'{sid:016d}'
        return [QLP_URL.format(s=u, a=z[:4], b=z[4:8], c=z[8:12], d=z[12:], tic=sid, v=v) for v in versions]
    if u >= 5000:
        raise ValueError('ZTF light curves are fetched by position, see occulta.ztf')
    if u >= 4000:
        z = f'{sid:016d}'
        return [TESS_SPOC_URL.format(s=u - 4000, a=z[:4], b=z[4:8], c=z[8:12], d=z[12:], tic=sid)]
    if u >= 3000:
        k9 = f'{sid - KEP_OFF:09d}'
        return [KEPLER_URL.format(k4=k9[:4], k9=k9, ts=KEPLER_STAMPS[u - 3000])]
    raise ValueError('SPOC and K2 urls come from the unit script')


def kepler_targets():
    """All Kepler targets with DR25 stellar properties (NASA Exoplanet Archive)."""
    r = requests.get('https://exoplanetarchive.ipac.caltech.edu/TAP/sync', params={
        'query': 'select kepid, kepmag, teff, radius, mass, ra, dec, dist from q1_q17_dr25_stellar', 'format': 'csv'},
        timeout=900)
    import csv
    rows = list(csv.DictReader(io.StringIO(r.text)))
    return rows
