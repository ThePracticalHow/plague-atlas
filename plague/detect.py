"""
Pathogen detection from expression data. No BAM needed. Works on ANY count matrix.

Three scores per cell:
    ALARM   Fungal detection genes (chitinases, Dectin-1, TLR2/4, CARD9)
    SPORE   Lysosomal expansion genes (cathepsins, LAMP, TFEB, V-ATPase)
    MELANIN Melanin machinery genes (TYR, PMEL, DCT)

Usage:
    plague detect sample.h5ad
    plague detect sample.h5ad --group cell_type
"""
import json
import time
import numpy as np
from collections import Counter

ALARM_GENES = {
    'CHI3L1': 'chitinase-3-like-1 (chitin detector)',
    'CHI3L2': 'chitinase-3-like-2',
    'CHIT1': 'chitotriosidase (chitin degrader)',
    'CLEC7A': 'Dectin-1 (beta-glucan receptor)',
    'CLEC4E': 'Mincle (trehalose receptor)',
    'TLR2': 'toll-like receptor 2 (fungal PAMP)',
    'TLR4': 'toll-like receptor 4 (fungal mannan)',
    'CARD9': 'CARD9 (Dectin downstream)',
    'IL17A': 'IL-17A (anti-fungal cytokine)',
    'IL17F': 'IL-17F (anti-fungal cytokine)',
    'DEFB1': 'beta-defensin 1 (anti-fungal peptide)',
    'PTX3': 'pentraxin 3 (Aspergillus opsonization)',
    'MBL2': 'mannose-binding lectin',
}

SPORE_GENES = {
    'LAMP1': 'lysosomal membrane 1',
    'LAMP2': 'lysosomal membrane 2',
    'CTSA': 'cathepsin A', 'CTSB': 'cathepsin B', 'CTSC': 'cathepsin C',
    'CTSD': 'cathepsin D', 'CTSE': 'cathepsin E', 'CTSF': 'cathepsin F',
    'CTSG': 'cathepsin G', 'CTSH': 'cathepsin H', 'CTSK': 'cathepsin K',
    'CTSL': 'cathepsin L', 'CTSS': 'cathepsin S', 'CTSZ': 'cathepsin Z',
    'ATP6V0A1': 'V-ATPase a1', 'ATP6V0B': 'V-ATPase b',
    'ATP6V1A': 'V-ATPase A', 'ATP6V1B2': 'V-ATPase B2',
    'HEXA': 'hexosaminidase A', 'HEXB': 'hexosaminidase B',
    'GBA': 'glucocerebrosidase', 'GAA': 'acid alpha-glucosidase',
    'NPC1': 'cholesterol transport', 'NPC2': 'cholesterol transport 2',
    'TFEB': 'lysosome master transcription factor',
    'CD63': 'lysosomal/melanosomal membrane',
    'MCOLN1': 'TRPML1 (lysosomal calcium)',
}

MELANIN_GENES = {
    'PMEL': 'premelanosome protein',
    'TYR': 'tyrosinase (melanin synthesis)',
    'TYRP1': 'tyrosinase-related protein 1',
    'DCT': 'dopachrome tautomerase',
    'MLANA': 'melan-A',
    'SLC45A2': 'melanin transporter',
    'OCA2': 'melanocyte-specific transporter',
    'MC1R': 'melanocortin 1 receptor',
}

AUTO_GROUP_KEYS = ['cell_type', 'condition', 'leiden', 'louvain', 'sample', 'group']


def _score_panel(X, var_names, panel_genes, sparse=False):
    """Score each cell on a gene panel. Returns per-cell sum and per-gene stats."""
    genes_found = {}
    gene_indices = []
    for g in panel_genes:
        if g in var_names:
            idx = var_names.index(g)
            genes_found[g] = idx
            gene_indices.append(idx)

    if not gene_indices:
        return np.zeros(X.shape[0]), genes_found, {}

    idx_arr = np.array(gene_indices)
    panel_expr = np.asarray(X[:, idx_arr].toarray()) if sparse else X[:, idx_arr]
    per_cell_score = panel_expr.sum(axis=1).flatten()

    gene_stats = {}
    for i, gene in enumerate(genes_found):
        col = panel_expr[:, i]
        gene_stats[gene] = {
            'pct_cells': round(float((col > 0).mean() * 100), 1),
            'mean_expression': round(float(col.mean()), 3),
            'description': panel_genes[gene],
        }

    return np.asarray(per_cell_score).flatten(), genes_found, gene_stats


def detect(h5ad_path, group_key=None, output_path=None, quiet=False):
    """Score every cell for pathogen presence indicators.

    Args:
        h5ad_path: path to .h5ad file.
        group_key: column in obs to group by.
        output_path: save JSON results here.
        quiet: suppress print output.

    Returns:
        dict with panel scores, per-group breakdown, and gene stats.
    """
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError(
            "scanpy required for h5ad files. "
            "Install with: pip install plague-atlas[h5ad]"
        )

    t0 = time.time()
    if not quiet:
        print(f'\n  Plague Detector')
        print(f'  Input: {h5ad_path}')

    adata = sc.read_h5ad(h5ad_path)
    n_cells, n_genes = adata.shape
    if not quiet:
        print(f'  Shape: {n_cells:,} cells x {n_genes:,} genes')

    var_names = list(adata.var_names)
    sparse = hasattr(adata.X, 'toarray')

    alarm_scores, alarm_found, alarm_stats = _score_panel(adata.X, var_names, ALARM_GENES, sparse)
    spore_scores, spore_found, spore_stats = _score_panel(adata.X, var_names, SPORE_GENES, sparse)
    melanin_scores, melanin_found, melanin_stats = _score_panel(adata.X, var_names, MELANIN_GENES, sparse)

    if not quiet:
        print(f'\n  Panel coverage:')
        print(f'    ALARM  (fungal detection):   {len(alarm_found):2d}/{len(ALARM_GENES)} genes')
        print(f'    SPORE  (lysosomal expansion): {len(spore_found):2d}/{len(SPORE_GENES)} genes')
        print(f'    MELANIN (melanin machinery): {len(melanin_found):2d}/{len(MELANIN_GENES)} genes')

    def _norm(arr):
        mx = arr.max()
        return arr / mx if mx > 0 else arr

    composite = _norm(alarm_scores) + _norm(spore_scores) + _norm(melanin_scores)

    alarm_positive = alarm_scores > 0
    spore_positive = (spore_scores > np.median(spore_scores[spore_scores > 0])
                      if (spore_scores > 0).any() else spore_scores > 0)
    high_composite = composite > np.percentile(composite, 75)

    if not quiet:
        print(f'\n  Results ({n_cells:,} cells):')
        print(f'    ALARM positive: {alarm_positive.sum():,} ({100*alarm_positive.mean():.1f}%)')
        print(f'    SPORE high:     {spore_positive.sum():,} ({100*spore_positive.mean():.1f}%)')
        print(f'    Composite top25%: {high_composite.sum():,} ({100*high_composite.mean():.1f}%)')

    if group_key and group_key in adata.obs.columns:
        groups = adata.obs[group_key].astype(str).values
    else:
        for key in AUTO_GROUP_KEYS:
            if key in adata.obs.columns:
                groups = adata.obs[key].astype(str).values
                group_key = key
                break
        else:
            groups = np.array(['ALL'] * n_cells)
            group_key = 'ALL'

    group_results = {}
    gcounts = Counter(groups)

    if not quiet:
        print(f'\n  Per-group scores (by {group_key}):')
        print(f'  {"Group":30s}  {"Cells":>7}  {"Alarm%":>7}  {"Spore":>7}  {"Composite":>9}')
        print(f'  {"-"*70}')

    for grp in sorted(gcounts, key=lambda g: -gcounts[g]):
        mask = groups == grp
        n = mask.sum()
        if n < 10:
            continue
        alarm_pct = 100 * alarm_positive[mask].mean()
        spore_mean = spore_scores[mask].mean()
        composite_mean = composite[mask].mean()
        group_results[grp] = {
            'n_cells': int(n),
            'alarm_pct': round(float(alarm_pct), 1),
            'spore_mean': round(float(spore_mean), 2),
            'composite_mean': round(float(composite_mean), 3),
        }
        if not quiet:
            print(f'  {grp[:29]:30s}  {n:>7,}  {alarm_pct:>6.1f}%  {spore_mean:>7.1f}  {composite_mean:>9.3f}')

    elapsed = time.time() - t0

    results = {
        'tool': 'plague_detector',
        'version': '1.0',
        'source': str(h5ad_path),
        'n_cells': n_cells,
        'n_genes': n_genes,
        'panels': {
            'alarm': {'genes_found': len(alarm_found), 'genes_total': len(ALARM_GENES), 'gene_stats': alarm_stats},
            'spore': {'genes_found': len(spore_found), 'genes_total': len(SPORE_GENES), 'gene_stats': spore_stats},
            'melanin': {'genes_found': len(melanin_found), 'genes_total': len(MELANIN_GENES), 'gene_stats': melanin_stats},
        },
        'summary': {
            'alarm_positive_pct': round(float(100 * alarm_positive.mean()), 1),
            'spore_high_pct': round(float(100 * spore_positive.mean()), 1),
            'composite_high_pct': round(float(100 * high_composite.mean()), 1),
        },
        'per_group': group_results,
        'compute_time': round(elapsed, 1),
    }

    if output_path:
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        if not quiet:
            print(f'\n  Saved: {output_path}')

    return results


def detect_from_matrix(X, gene_names, groups=None, group_key=None, quiet=False):
    """Score cells from a raw count matrix (no scanpy needed).

    Args:
        X: (n_cells, n_genes) array or sparse matrix.
        gene_names: list of gene name strings.
        groups: optional group labels array.
        group_key: name for grouping.
        quiet: suppress output.

    Returns:
        dict with panel scores per group.
    """
    var_names = list(gene_names)
    sparse = hasattr(X, 'toarray')
    n_cells = X.shape[0]

    alarm_scores, alarm_found, alarm_stats = _score_panel(X, var_names, ALARM_GENES, sparse)
    spore_scores, spore_found, spore_stats = _score_panel(X, var_names, SPORE_GENES, sparse)
    melanin_scores, melanin_found, melanin_stats = _score_panel(X, var_names, MELANIN_GENES, sparse)

    def _norm(arr):
        mx = arr.max()
        return arr / mx if mx > 0 else arr

    composite = _norm(alarm_scores) + _norm(spore_scores) + _norm(melanin_scores)
    alarm_positive = alarm_scores > 0

    if groups is None:
        groups = np.array(['ALL'] * n_cells)
        group_key = group_key or 'ALL'
    else:
        groups = np.asarray(groups, dtype=str)
        group_key = group_key or 'group'

    group_results = {}
    for grp in sorted(set(groups)):
        mask = groups == grp
        n = mask.sum()
        if n < 10:
            continue
        group_results[grp] = {
            'n_cells': int(n),
            'alarm_pct': round(float(100 * alarm_positive[mask].mean()), 1),
            'spore_mean': round(float(spore_scores[mask].mean()), 2),
            'composite_mean': round(float(composite[mask].mean()), 3),
        }

    return {
        'tool': 'plague_detector',
        'version': '1.0',
        'n_cells': n_cells,
        'panels': {
            'alarm': {'genes_found': len(alarm_found), 'genes_total': len(ALARM_GENES)},
            'spore': {'genes_found': len(spore_found), 'genes_total': len(SPORE_GENES)},
            'melanin': {'genes_found': len(melanin_found), 'genes_total': len(MELANIN_GENES)},
        },
        'per_group': group_results,
    }
