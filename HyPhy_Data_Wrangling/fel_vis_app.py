#!/usr/bin/env python3
"""
Streamlit visualization app for HyPhy FEL output.

Launch:
  make vis
  # or directly:
  .venv/bin/streamlit run fel_vis_app.py
"""

import io
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
if 'highlight_color' not in st.session_state:
    st.session_state['highlight_color'] = '#ff4444'
if 'bg_color' not in st.session_state:
    st.session_state['bg_color'] = '#1a1a2e'
if 'ribbon_color' not in st.session_state:
    st.session_state['ribbon_color'] = '#888888'

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
new_omega = omega_hist_comp(
    values=omega_finite.tolist(),
    cutoff=st.session_state.omega_cutoff,
    xaxis_title='dN/dS (ω)',
    key='omega_hist',
)
if new_omega is not None:
    st.session_state.omega_cutoff = new_omega
    st.session_state['omega_precise_input'] = float(new_omega)

omega_cutoff = st.session_state.omega_cutoff

# ── Precise cutoff input ───────────────────────────────────────────────────────
omega_input = st.number_input(
    'Precise dN/dS cutoff:',
    min_value=0.0, max_value=float(omega_max),
    step=0.001, format='%.3f',
    key='omega_precise_input',
)
if abs(omega_input - st.session_state.omega_cutoff) > 1e-9:
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
  </style>
</head>
<body>
  <div id="wrap">
    <div id="viewer"></div>
    <button id="save-btn" onclick="savePNG()">Save PNG</button>
  </div>
  <script>
    const cifData = `{cif_escaped}`;
    const allRes  = {all_res_js};

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
