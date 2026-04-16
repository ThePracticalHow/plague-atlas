"""
Food Atlas — Score every food by whose side it's on.

For each food compound: antifungal (helps the host) or profungal (helps the pathogen)?
Based on published MIC/IC50 data, membrane disruption studies, and immune modulation.

    from plague import FoodAtlas
    atlas = FoodAtlas()
    atlas.score_food("garlic")       # high positive (allicin)
    atlas.score_food("white sugar")  # negative (pure carbon for fungi)
"""

from typing import Dict, List, Tuple

ANTIFUNGAL_COMPOUNDS = {
    'allicin': {
        'mechanism': 'Membrane disruption. Kills 20+ Candida + Aspergillus species.',
        'targets': ['Candida', 'Aspergillus', 'Cryptococcus'],
        'potency': 9,
        'food_sources': ['garlic', 'onion', 'leek', 'shallot', 'chives'],
        'citation': 'Frontiers Microbiol 2022; PubMed 15743425',
    },
    'carvacrol': {
        'mechanism': 'Calcium stress mimicry. TOR pathway inhibition. Membrane disruption.',
        'targets': ['Candida', 'Aspergillus', 'broad spectrum'],
        'potency': 8,
        'food_sources': ['oregano', 'thyme', 'savory', 'marjoram'],
        'citation': 'PMC2981246',
    },
    'thymol': {
        'mechanism': 'Membrane disruption. Synergistic with carvacrol.',
        'targets': ['Candida', 'Aspergillus', 'dermatophytes'],
        'potency': 7,
        'food_sources': ['thyme', 'oregano', 'basil'],
    },
    'caprylic_acid': {
        'mechanism': 'Lipid bilayer integration. Membrane permeabilization.',
        'targets': ['Candida'],
        'potency': 7,
        'food_sources': ['coconut oil', 'palm kernel oil', 'butter', 'human breast milk'],
    },
    'lauric_acid': {
        'mechanism': 'Membrane disruption. Most abundant MCFA in coconut oil.',
        'targets': ['Candida', 'broad spectrum'],
        'potency': 7,
        'food_sources': ['coconut oil', 'palm kernel oil', 'human breast milk'],
    },
    'curcumin': {
        'mechanism': 'Ergosterol synthesis disruption. ROS generation.',
        'targets': ['Candida', 'Aspergillus', 'Penicillium'],
        'potency': 6,
        'food_sources': ['turmeric', 'curry powder', 'mustard'],
    },
    'cinnamaldehyde': {
        'mechanism': 'Membrane + cell wall disruption.',
        'targets': ['Candida', 'Aspergillus'],
        'potency': 7,
        'food_sources': ['cinnamon', 'cassia'],
    },
    'eugenol': {
        'mechanism': 'Membrane disruption. Ergosterol binding.',
        'targets': ['Candida', 'Aspergillus', 'Penicillium'],
        'potency': 6,
        'food_sources': ['cloves', 'allspice', 'basil', 'cinnamon', 'nutmeg'],
    },
    'acetic_acid': {
        'mechanism': 'pH disruption. Membrane permeabilization.',
        'targets': ['Candida', 'broad spectrum'],
        'potency': 4,
        'food_sources': ['vinegar', 'fermented vegetables', 'kombucha'],
    },
    'tannins': {
        'mechanism': 'Protein precipitation. Iron chelation.',
        'targets': ['Candida', 'broad spectrum'],
        'potency': 5,
        'food_sources': ['black tea', 'green tea', 'red wine', 'pomegranate',
                         'cranberry', 'walnut', 'dark chocolate'],
    },
    'lactoferrin': {
        'mechanism': 'Iron sequestration. Starves fungi of essential iron.',
        'targets': ['Candida', 'Aspergillus', 'broad spectrum'],
        'potency': 7,
        'food_sources': ['breast milk', 'colostrum', 'raw milk', 'whey'],
    },
    'berberine': {
        'mechanism': 'AMPK activator. Membrane disruption. DNA intercalation.',
        'targets': ['Candida', 'Cryptococcus'],
        'potency': 7,
        'food_sources': ['goldenseal', 'Oregon grape', 'barberry'],
    },
    'beta_glucan': {
        'mechanism': 'Immune stimulation. Activates Dectin-1 (anti-fungal receptor).',
        'targets': ['immune boost against all fungi'],
        'potency': 5,
        'food_sources': ['oats', 'barley', 'mushrooms', 'nutritional yeast'],
    },
    'sulforaphane': {
        'mechanism': 'NRF2 activation. Some antifungal activity.',
        'targets': ['Candida', 'broad spectrum'],
        'potency': 4,
        'food_sources': ['broccoli', 'broccoli sprouts', 'cauliflower', 'kale',
                         'brussels sprouts', 'cabbage'],
    },
}

PROFUNGAL_COMPOUNDS = {
    'glucose': {
        'mechanism': 'Primary carbon source for fungi.',
        'potency': -8,
        'food_sources': ['white sugar', 'candy', 'soft drinks', 'fruit juice',
                         'white bread', 'white rice', 'corn syrup', 'honey'],
    },
    'fructose': {
        'mechanism': 'Fermentable sugar. Feeds Candida.',
        'potency': -7,
        'food_sources': ['fruit juice', 'high-fructose corn syrup', 'agave',
                         'soft drinks', 'dried fruit'],
    },
    'refined_starch': {
        'mechanism': 'Rapidly converted to glucose.',
        'potency': -5,
        'food_sources': ['white bread', 'white rice', 'pasta', 'crackers',
                         'processed cereals', 'pastries'],
    },
    'ethanol': {
        'mechanism': 'Immunosuppressive. Gut barrier disruption. Feeds yeast.',
        'potency': -7,
        'food_sources': ['beer', 'wine', 'spirits', 'cocktails'],
    },
}


class FoodAtlas:
    """Score foods by whose side they're on in the host-pathogen war."""

    def __init__(self):
        self.antifungal = ANTIFUNGAL_COMPOUNDS
        self.profungal = PROFUNGAL_COMPOUNDS
        self._index = self._build_index()

    def _build_index(self) -> Dict[str, Dict]:
        index = {}
        for name, data in self.antifungal.items():
            for food in data['food_sources']:
                fl = food.lower()
                if fl not in index:
                    index[fl] = {'antifungal': [], 'profungal': [], 'score': 0}
                index[fl]['antifungal'].append((name, data['potency']))
                index[fl]['score'] += data['potency']

        for name, data in self.profungal.items():
            for food in data['food_sources']:
                fl = food.lower()
                if fl not in index:
                    index[fl] = {'antifungal': [], 'profungal': [], 'score': 0}
                index[fl]['profungal'].append((name, data['potency']))
                index[fl]['score'] += data['potency']
        return index

    def score_food(self, name: str) -> Dict:
        """Score a food: positive = fights fungi, negative = feeds fungi."""
        nl = name.lower()
        if nl in self._index:
            return {**self._index[nl], 'name': nl}

        matches = [k for k in self._index if nl in k or k in nl]
        if matches:
            best = max(matches, key=lambda k: abs(self._index[k]['score']))
            return {**self._index[best], 'name': best, 'matched_from': name}

        return {'antifungal': [], 'profungal': [], 'score': 0, 'name': name, 'unknown': True}

    def top_antifungal(self, n: int = 20) -> List[Tuple[str, int]]:
        ranked = sorted(self._index.items(), key=lambda x: -x[1]['score'])
        return [(name, d['score']) for name, d in ranked[:n] if d['score'] > 0]

    def top_profungal(self, n: int = 20) -> List[Tuple[str, int]]:
        ranked = sorted(self._index.items(), key=lambda x: x[1]['score'])
        return [(name, d['score']) for name, d in ranked[:n] if d['score'] < 0]

    def war_report(self):
        """Print the full food war report."""
        print("=" * 60)
        print("  FOOD WAR REPORT")
        print("=" * 60)

        print("\n  ANTIFUNGAL FOODS (fight the pathogen):")
        for name, score in self.top_antifungal(15):
            bar = '+' * score
            print(f"    {score:+3d}  {name:30s}  {bar}")

        print("\n  PROFUNGAL FOODS (feed the pathogen):")
        for name, score in self.top_profungal(15):
            bar = '-' * abs(score)
            print(f"    {score:+3d}  {name:30s}  {bar}")

        total = len(self._index)
        anti = sum(1 for d in self._index.values() if d['score'] > 0)
        pro = sum(1 for d in self._index.values() if d['score'] < 0)
        print(f"\n  Total indexed: {total} foods ({anti} antifungal, {pro} profungal)")
