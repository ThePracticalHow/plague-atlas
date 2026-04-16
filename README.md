# plague

Pathogen detection engine. Fungal scanner, gene panel diagnosis, dietary war scoring.

## What it does

`plague` detects fungal pathogen indicators in single-cell expression data and scores foods by their antifungal/profungal properties.

### Detection panels

| Panel | Genes | What it detects |
|-------|-------|-----------------|
| **ALARM** | CHI3L1, CLEC7A, TLR2, CARD9, ... | Host fungal detection system activation |
| **SPORE** | LAMP1/2, cathepsins, TFEB, V-ATPase, ... | Lysosomal expansion (pathogen growth signature) |
| **MELANIN** | TYR, PMEL, DCT, TYRP1, ... | Melanin machinery (pathogen pigmentation) |

### Food war scoring

Every food scored by published antifungal compound data:
- **Garlic** (+17): allicin, ajoene (membrane disruption)
- **Oregano** (+15): carvacrol, thymol (TOR inhibition)
- **Coconut oil** (+14): caprylic, capric, lauric acid
- **White sugar** (-8): pure carbon for fungi
- **Beer** (-7): immunosuppressive, feeds yeast

## Install

```bash
pip install plague-atlas

# For h5ad support:
pip install plague-atlas[h5ad]
```

## Usage

### Command line

```bash
# Detect pathogen indicators in expression data
plague detect data.h5ad
plague detect data.h5ad --group cell_type

# Score a food
plague food garlic
plague food "white sugar"
plague food --report

# Show gene panels
plague panels

# Full pipeline: coupling tensor + plague detection
plague full data.h5ad --group cell_type
```

### Python

```python
from plague import detect, FoodAtlas

# Gene panel detection (requires scanpy)
results = detect("data.h5ad", group_key="cell_type")

# From raw matrix (no scanpy needed)
from plague.detect import detect_from_matrix
results = detect_from_matrix(X, gene_names, groups=labels)

# Food scoring
atlas = FoodAtlas()
atlas.score_food("garlic")       # {'score': 17, ...}
atlas.score_food("white sugar")  # {'score': -8, ...}
atlas.war_report()               # Full report
```

## Dependencies

- `numpy`, `scipy` (always)
- `coherence-tensor` (for combined tensor+detection pipeline)
- `scanpy` (optional, for h5ad files)

## License

MIT
