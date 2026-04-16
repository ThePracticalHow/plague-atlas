"""
plague — Pathogen detection engine.

    plague detect data.h5ad          # Score cells for fungal presence
    plague food garlic               # Score a food's war alignment
    plague scan data.bam             # Scan BAM for fungal transcripts

Three detection systems:
    ALARM   Fungal detection genes (CHI3L1, CLEC7A, TLR2, CARD9, ...)
    SPORE   Lysosomal expansion genes (cathepsins, LAMP, TFEB, ...)
    MELANIN Melanin machinery genes (TYR, PMEL, DCT, ...)

Plus: food war scoring (antifungal vs profungal compounds).
"""

__version__ = "0.1.0"

from .detect import detect, ALARM_GENES, SPORE_GENES, MELANIN_GENES
from .food import FoodAtlas
