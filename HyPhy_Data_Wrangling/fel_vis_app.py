#!/usr/bin/env python3
"""
Streamlit visualization app for HyPhy FEL output.

Launch:
  make vis
  # or directly:
  .venv/bin/streamlit run fel_vis_app.py
"""

import io
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
import streamlit.components.v1 as components

sys.path.insert(0, str(Path(__file__).parent))
from parse_fel import parse_fel_output
from map_codons import build_mapping

# ── Paths ──────────────────────────────────────────────────────────────────────
FEL_INPUT      = Path(__file__).parent.parent / 'HyPhy' / 'HyPhy_FEL_output.txt'
ALIGNMENT      = Path(__file__).parent.parent / 'HyPhy' / 'codon_align.faa'
CIF_PATH       = Path(__file__).parent.parent / 'StructuralModeling' / '8adg.cif'
HIST_COMPONENT = Path(__file__).resolve().parent / 'hist_component'

st.set_page_config(page_title='BamA FEL Explorer', layout='wide')
st.title('BamA Site Selection Explorer')


@st.cache_data
def load_fel(path: str):
    return parse_fel_output(Path(path))

@st.cache_data
def load_mapping(alignment: str) -> dict:
    return build_mapping(Path(alignment))

@st.cache_data
def load_cif(path: str) -> str:
    return Path(path).read_text()


df      = load_fel(str(FEL_INPUT))
mapping = load_mapping(str(ALIGNMENT))
cif_str = load_cif(str(CIF_PATH))

# ── Hover lookup: cif_residue -> FEL data (all sites in the FEL table) ────────
@st.cache_data
def build_hover_data(fel_path: str, alignment: str) -> dict:
    _df  = parse_fel_output(Path(fel_path))
    _map = build_mapping(Path(alignment))
    codon_to_cif = {codon: cif for codon, cif in _map.items() if cif is not None}
    out = {}
    for _, row in _df.iterrows():
        cif = codon_to_cif.get(int(row['Codon']))
        if cif is None:
            continue
        a, b = float(row['alpha']), float(row['beta'])
        w = b / a if a > 0 else None
        out[cif] = {
            'codon':     int(row['Codon']),
            'alpha':     round(a, 4),
            'beta':      round(b, 4),
            'omega':     round(w, 4) if w is not None else None,
            'pvalue':    round(float(row['pvalue']), 4),
            'selection': row['Selection'],
        }
    return out

hover_lookup = build_hover_data(str(FEL_INPUT), str(ALIGNMENT))
hover_data_js = json.dumps(hover_lookup)

# ── Derived: dN/dS per site ────────────────────────────────────────────────────
omega        = np.where(df['alpha'] > 0, df['beta'] / df['alpha'], np.nan)
omega_finite = omega[np.isfinite(omega)]
omega_p95    = round(float(np.percentile(omega_finite, 95)), 2) if len(omega_finite) else 1.0
omega_max    = max(omega_p95, 1.0)

# ── Session state ──────────────────────────────────────────────────────────────
if 'omega_cutoff' not in st.session_state:
    st.session_state.omega_cutoff = omega_max
if 'omega_precise_input' not in st.session_state:
    st.session_state['omega_precise_input'] = float(omega_max)
if 'last_drag_t' not in st.session_state:
    st.session_state['last_drag_t'] = -1
if 'highlight_color' not in st.session_state:
    st.session_state['highlight_color'] = '#F05E5E'
if 'bg_color' not in st.session_state:
    st.session_state['bg_color'] = '#DFDEDE'
if 'ribbon_color' not in st.session_state:
    st.session_state['ribbon_color'] = '#000000'

# ── Sidebar filters ────────────────────────────────────────────────────────────
st.sidebar.header('Filter Controls')

selection_type = st.sidebar.radio('Selection type', ['both', 'Neg', 'Pos'])

alpha_p95 = float(np.percentile(df['alpha'], 95))
alpha_range = st.sidebar.slider(
    'Alpha range (synonymous rate)',
    min_value=0.0, max_value=round(alpha_p95, 2),
    value=(0.0, round(alpha_p95, 2)), step=0.01,
)
st.sidebar.caption(f'Full alpha max: {df["alpha"].max():.1f}')

beta_max_val = round(float(df['beta'].max()), 2)
beta_range = st.sidebar.slider(
    'Beta range (non-synonymous rate)',
    min_value=0.0, max_value=beta_max_val,
    value=(0.0, beta_max_val), step=0.01,
)

lrt_max_val = round(float(df['LRT'].max()), 1)
lrt_range = st.sidebar.slider(
    'LRT range',
    min_value=0.0, max_value=lrt_max_val,
    value=(0.0, lrt_max_val), step=0.1,
)

pvalue_cutoff = st.sidebar.slider(
    'P-value cutoff',
    min_value=0.0, max_value=1.0,
    value=0.1, step=0.001,
    format='%.3f',
)

# ── dN/dS histogram component ──────────────────────────────────────────────────
omega_hist_comp = components.declare_component('omega_hist', path=str(HIST_COMPONENT))

st.subheader('dN/dS Distribution — click or drag to set cutoff')
hist_result = omega_hist_comp(
    values=omega_finite.tolist(),
    cutoff=st.session_state.omega_cutoff,
    xaxis_title='dN/dS (ω)',
    key='omega_hist',
)
new_omega = None
if isinstance(hist_result, dict):
    drag_t = hist_result.get('t', -1)
    if drag_t != st.session_state['last_drag_t']:
        new_omega = float(hist_result['cutoff'])
        st.session_state['last_drag_t'] = drag_t
        st.session_state.omega_cutoff = new_omega
        st.session_state['omega_precise_input'] = new_omega

omega_cutoff = st.session_state.omega_cutoff

# ── Precise cutoff input ───────────────────────────────────────────────────────
omega_input = st.number_input(
    'Precise dN/dS cutoff:',
    min_value=0.0, max_value=float(omega_max),
    step=0.0001, format='%.6f',
    key='omega_precise_input',
)
if new_omega is None and abs(omega_input - st.session_state.omega_cutoff) > 1e-9:
    st.session_state.omega_cutoff = omega_input
    st.rerun()

# ── Blurb ──────────────────────────────────────────────────────────────────────
st.markdown(
    'The **α** and **β** values here come directly from `HyPhy_FEL_output.txt` — a table '
    'HyPhy produced by running its Fixed Effects Likelihood (FEL) method on a codon alignment '
    'of 99 bacterial BamA sequences. For each of the 1,191 codon positions in that alignment, '
    'HyPhy estimated two per-site substitution rates using maximum likelihood across the full '
    'phylogenetic tree: **α** is the synonymous rate (DNA changes that do not alter the amino '
    'acid, e.g. AAA→AAG, both Lys) and **β** is the non-synonymous rate (changes that do, '
    'e.g. AAA→GAA, Lys→Glu). Only sites where β significantly differed from α (LRT p ≤ 0.1) '
    'were written into that output file — those are the only sites visible here. '
    '**This app then computes ω = β/α** from those stored values. A site where the amino acid '
    'is identical across all 99 species has β ≈ 0, so ω ≈ 0 — the signature of strong '
    '**purifying selection**. Dragging the cutoff left isolates the residues most intolerant '
    'of amino acid change across the sampled bacterial lineages.'
)

# ── Apply filters ──────────────────────────────────────────────────────────────
mask = (
    (df['pvalue'] <= pvalue_cutoff)
    & (df['alpha'] >= alpha_range[0]) & (df['alpha'] <= alpha_range[1])
    & (df['beta']  >= beta_range[0])  & (df['beta']  <= beta_range[1])
    & (df['LRT']   >= lrt_range[0])   & (df['LRT']   <= lrt_range[1])
    & ((omega <= omega_cutoff) | ~np.isfinite(omega))
)
if selection_type != 'both':
    mask &= df['Selection'] == selection_type

filtered = df[mask]
neg      = filtered[filtered['Selection'] == 'Neg']
pos      = filtered[filtered['Selection'] == 'Pos']

def to_cif(codons):
    return [mapping[c] for c in codons if mapping.get(c) is not None]

all_res = to_cif(filtered['Codon'].tolist())

st.markdown(
    f'**{len(filtered)}** sites selected &nbsp;|&nbsp; '
    f'{len(neg)} negative &nbsp;|&nbsp; {len(pos)} positive'
)

# ── 3Dmol viewer ───────────────────────────────────────────────────────────────
highlight_color = st.session_state['highlight_color']
bg_color        = st.session_state['bg_color']
ribbon_color    = st.session_state['ribbon_color']

cif_escaped     = cif_str.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
all_res_js      = str(all_res)
bg_color_js     = bg_color.replace("'", "\\'")
ribbon_color_js = ribbon_color.replace("'", "\\'")
hl_color_js     = highlight_color.replace("'", "\\'")

viewer_html = f"""
<!DOCTYPE html>
<html>
<head>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    body {{ margin: 0; background: {bg_color}; }}
    #wrap {{ position: relative; width: 100%; height: 600px; }}
    #viewer {{ width: 100%; height: 100%; }}
    #save-btn {{
      position: absolute; top: 10px; right: 10px; z-index: 100;
      background: rgba(255,255,255,0.12); color: #ddd;
      border: 1px solid rgba(255,255,255,0.25);
      padding: 5px 12px; border-radius: 4px;
      cursor: pointer; font-size: 13px; font-family: sans-serif;
      backdrop-filter: blur(4px);
    }}
    #save-btn:hover {{ background: rgba(255,255,255,0.22); color: #fff; }}
    #tooltip {{
      position: absolute; display: none; z-index: 200;
      background: rgba(10,10,20,0.92); color: #e8e8e8;
      border: 1px solid rgba(255,255,255,0.18);
      padding: 8px 11px; border-radius: 6px;
      font-family: monospace; font-size: 12px; line-height: 1.6;
      pointer-events: none; white-space: nowrap;
      box-shadow: 0 2px 12px rgba(0,0,0,0.5);
    }}
    #tooltip .tt-header {{ font-size: 13px; font-weight: bold; color: #fff; margin-bottom: 3px; }}
    #tooltip .tt-sel-neg  {{ color: #7bcfff; }}
    #tooltip .tt-sel-pos  {{ color: #ffaa55; }}
    #tooltip .tt-dim {{ color: #888; font-size: 11px; }}
  </style>
</head>
<body>
  <div id="wrap">
    <div id="viewer"></div>
    <div id="tooltip"></div>
    <button id="save-btn" onclick="savePNG()">Save PNG</button>
  </div>
  <script>
    const cifData    = `{cif_escaped}`;
    const allRes     = {all_res_js};
    const hoverData  = {hover_data_js};
    const tooltip    = document.getElementById('tooltip');
    const wrap       = document.getElementById('wrap');

    let viewer = $3Dmol.createViewer('viewer', {{ backgroundColor: '{bg_color_js}' }});
    viewer.addModel(cifData, 'cif');

    viewer.setStyle({{ chain: 'A' }}, {{ cartoon: {{ color: '{ribbon_color_js}', opacity: 0.6 }} }});
    viewer.setStyle({{ not: {{ chain: 'A' }} }}, {{}});

    if (allRes.length > 0) {{
      viewer.addStyle(
        {{ chain: 'A', resi: allRes }},
        {{ cartoon: {{ color: '{hl_color_js}', opacity: 1.0 }},
           sphere:  {{ radius: 0.4, color: '{hl_color_js}' }} }}
      );
    }}

    viewer.setHoverable({{}}, true,
      function(atom, v, event) {{
        if (atom.chain !== 'A') return;
        const d = hoverData[atom.resi];
        let html = `<div class="tt-header">${{atom.resn}} ${{atom.resi}}</div>`;
        if (d) {{
          const selClass = d.selection === 'Neg' ? 'tt-sel-neg' : 'tt-sel-pos';
          const selLabel = d.selection === 'Neg' ? '&#8595; Purifying' : '&#8593; Positive';
          const omegaStr = d.omega !== null ? d.omega.toFixed(4) : '<span class="tt-dim">undef (α=0)</span>';
          html += `<span class="tt-dim">codon ${{d.codon}}</span><br>`;
          html += `&#945; ${{d.alpha.toFixed(4)}} &nbsp; &#946; ${{d.beta.toFixed(4)}}<br>`;
          html += `&#969; ${{omegaStr}} &nbsp; p ${{d.pvalue.toFixed(4)}}<br>`;
          html += `<span class="${{selClass}}">${{selLabel}}</span>`;
        }} else {{
          html += `<span class="tt-dim">no FEL data</span>`;
        }}
        tooltip.innerHTML = html;
        tooltip.style.display = 'block';
        const wrapRect = wrap.getBoundingClientRect();
        let x = event.clientX - wrapRect.left + 14;
        let y = event.clientY - wrapRect.top  + 14;
        // keep tooltip inside wrap
        if (x + 180 > wrap.offsetWidth)  x = event.clientX - wrapRect.left - 190;
        if (y + 120 > wrap.offsetHeight) y = event.clientY - wrapRect.top  - 130;
        tooltip.style.left = x + 'px';
        tooltip.style.top  = y + 'px';
      }},
      function() {{
        tooltip.style.display = 'none';
      }}
    );

    viewer.zoomTo({{ chain: 'A' }});
    viewer.render();

    function savePNG() {{
      const uri = viewer.pngURI();
      const a = document.createElement('a');
      a.href = uri;
      a.download = 'bama_structure.png';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    }}
  </script>
</body>
</html>
"""

components.html(viewer_html, height=620, scrolling=False)

# ── Colors & Export ────────────────────────────────────────────────────────────
def make_hist_png(vals, cutoff, bar_color, line_color, background):
    fig, ax = plt.subplots(figsize=(8, 3))
    fig.patch.set_facecolor(background)
    ax.set_facecolor(background)
    ax.hist(vals, bins=60, color=bar_color, edgecolor='none')
    ax.axvline(cutoff, color=line_color, linewidth=2)
    ax.text(cutoff, ax.get_ylim()[1] * 0.97, f' ω = {cutoff:.3f}',
            color=line_color, fontsize=9, va='top')
    ax.set_xlabel('dN/dS (ω)', color='#cccccc')
    ax.set_ylabel('count', color='#cccccc')
    ax.tick_params(colors='#cccccc')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('#444444')
    ax.spines['left'].set_color('#444444')
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()

with st.expander('Colors & Export'):
    cp1, cp2, cp3 = st.columns(3)
    with cp1:
        st.color_picker('Highlighted residues', highlight_color, key='highlight_color')
    with cp2:
        st.color_picker('Background', bg_color, key='bg_color')
    with cp3:
        st.color_picker('Ribbon', ribbon_color, key='ribbon_color')

    ex1, ex2 = st.columns(2)
    with ex1:
        hist_png = make_hist_png(
            omega_finite.tolist(), omega_cutoff,
            bar_color=ribbon_color, line_color=highlight_color, background=bg_color,
        )
        st.download_button(
            'Export Histogram PNG', hist_png,
            file_name='dnds_histogram.png', mime='image/png',
        )
    with ex2:
        st.caption('Use the **Save PNG** button in the top-right of the 3D viewer to export the structure.')

# ── Codon position output ──────────────────────────────────────────────────────
st.subheader('Codon positions')
codon_str = ','.join(str(c) for c in sorted(filtered['Codon'].tolist()))
st.text_area(
    label=f'{len(filtered)} positions — copy and paste into a text file',
    value=codon_str,
    height=120,
)
