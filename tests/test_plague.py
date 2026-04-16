"""Smoke tests for plague detection engine."""

import numpy as np
import pytest


def _make_synthetic_with_alarm(n_cells=300, seed=42):
    """Build synthetic matrix with known alarm/spore genes."""
    rng = np.random.RandomState(seed)

    alarm = ['CHI3L1', 'CLEC7A', 'TLR2', 'CARD9', 'IL17A']
    spore = ['LAMP1', 'LAMP2', 'CTSA', 'CTSB', 'CTSD', 'TFEB', 'HEXA', 'NPC1']
    melanin = ['PMEL', 'TYR', 'TYRP1', 'DCT']
    other = [f'GENE{i}' for i in range(1, 201)]

    genes = alarm + spore + melanin + other
    X = rng.poisson(2, (n_cells, len(genes))).astype(float)

    # First half: elevated alarm. Second half: zero alarm.
    for i in range(len(alarm)):
        X[:n_cells // 2, i] *= 5
        X[n_cells // 2:, i] = 0

    return X, genes


class TestDetectFromMatrix:
    def test_basic_detection(self):
        from plague.detect import detect_from_matrix
        X, genes = _make_synthetic_with_alarm()
        results = detect_from_matrix(X, genes, quiet=True)
        assert results['tool'] == 'plague_detector'
        assert results['n_cells'] == 300
        assert results['panels']['alarm']['genes_found'] > 0
        assert results['panels']['spore']['genes_found'] > 0

    def test_grouping(self):
        from plague.detect import detect_from_matrix
        X, genes = _make_synthetic_with_alarm()
        groups = np.array(['infected'] * 150 + ['healthy'] * 150)
        results = detect_from_matrix(X, genes, groups=groups, quiet=True)
        assert 'infected' in results['per_group']
        assert 'healthy' in results['per_group']
        assert results['per_group']['infected']['alarm_pct'] > results['per_group']['healthy']['alarm_pct']


class TestFoodAtlas:
    def test_init(self):
        from plague.food import FoodAtlas
        atlas = FoodAtlas()
        assert len(atlas._index) > 0

    def test_score_garlic(self):
        from plague.food import FoodAtlas
        atlas = FoodAtlas()
        result = atlas.score_food('garlic')
        assert result['score'] > 0
        assert len(result['antifungal']) > 0

    def test_score_sugar(self):
        from plague.food import FoodAtlas
        atlas = FoodAtlas()
        result = atlas.score_food('white sugar')
        assert result['score'] < 0

    def test_unknown_food(self):
        from plague.food import FoodAtlas
        atlas = FoodAtlas()
        result = atlas.score_food('xyznonexistent')
        assert result.get('unknown') is True

    def test_top_lists(self):
        from plague.food import FoodAtlas
        atlas = FoodAtlas()
        anti = atlas.top_antifungal(5)
        pro = atlas.top_profungal(5)
        assert len(anti) > 0
        assert len(pro) > 0
        assert all(score > 0 for _, score in anti)
        assert all(score < 0 for _, score in pro)


class TestCLI:
    def test_help(self):
        import subprocess
        result = subprocess.run(
            ['python', '-m', 'plague', '--help'],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0
        assert 'plague' in result.stdout.lower()

    def test_panels(self):
        import subprocess
        result = subprocess.run(
            ['python', '-m', 'plague', 'panels'],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0
        assert 'ALARM' in result.stdout

    def test_food_report(self):
        import subprocess
        result = subprocess.run(
            ['python', '-m', 'plague', 'food', '--report'],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0
        assert 'garlic' in result.stdout.lower()

    def test_food_score(self):
        import subprocess
        result = subprocess.run(
            ['python', '-m', 'plague', 'food', 'oregano'],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0
        assert 'Score:' in result.stdout

    def test_version(self):
        import subprocess
        result = subprocess.run(
            ['python', '-m', 'plague', 'version'],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0
        assert '0.1.0' in result.stdout
