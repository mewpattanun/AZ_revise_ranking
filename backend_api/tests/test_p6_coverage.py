"""
Unit tests for Priority 6: separate coverage from score.

AZ slide 17: "Never interpret missing protein as measured zero."

Covers:
- Step 2: RNA-only / Protein-only scores are un-capped (full biology magnitude)
- Step 3: API response carries modality_coverage as a separate field
- Step 4: require_protein filter excludes RNA-only cells
- Step 5: in a LOW off tier with mixed protein coverage, cells are tied
"""

import pytest

from app.services.ranking_service import run_ranking

HIGH = {"hugo": "G1", "ensembl_id": "ENSG_G1", "direction": "high", "weight": 1.0}
LOW = {"hugo": "G1", "ensembl_id": "ENSG_G1", "direction": "low", "weight": 1.0}


# Step 2: un-cap 

def test_rna_only_score_is_not_capped(fake_data):
    """
    P6 Step 2: an RNA-only cell with the strongest RNA signal should score
    higher than cells whose only strength is protein. Before P6, RNA-only
    scores were multiplied by w_rna (= capped). Now they run at full scale.
    """
    # A is top on RNA (DepMap 100, HPA 90). B is weak on RNA but has protein.
    fake_data["protein"] = {"B": 50.0}

    results = run_ranking([HIGH], w_rna=0.7, w_protein=0.3, top_n=10)
    top = next(r for r in results if r["ach_id"] == "A")

    # Full-scale: A's score is the biology magnitude itself, not w_rna * mag.
    # Pre-P6 cap would force RNA-only <= 0.7; now it reaches 1.0.
    assert top["scenario"] == "RNA only"
    assert top["score"] > 0.7, f"RNA-only score was capped: {top['score']}"


# Step 3: modality_coverage field 

def test_modality_coverage_reports_available_modalities(fake_data):
    """
    P6 Step 3: modality_coverage is the fraction of the user's requested
    modalities actually present for the cell line — not baked into score.
      RNA + Protein -> 1.0
      RNA only      -> w_rna
      Protein only  -> w_protein
    """
    # Protein for A and C only (B and D are RNA-only). rank_protein needs
    # at least 2 cells in the protein series to compute robust z.
    fake_data["protein"] = {"A": 50.0, "C": 30.0}

    results = run_ranking([HIGH], w_rna=0.7, w_protein=0.3, top_n=10)
    by_ach = {r["ach_id"]: r for r in results}

    assert by_ach["A"]["scenario"] == "RNA+Protein"
    assert by_ach["A"]["modality_coverage"] == 1.0

    assert by_ach["B"]["scenario"] == "RNA only"
    assert by_ach["B"]["modality_coverage"] == pytest.approx(0.7)


# Step 4: require_protein filter 

def test_require_protein_filter_drops_rna_only_cells(fake_data):
    """P6 Step 4: require_protein=True removes cells with no protein data."""
    fake_data["protein"] = {"A": 50.0, "C": 20.0}   # B and D have no protein

    results = run_ranking([HIGH], top_n=10, require_protein=True)
    scenarios = [r["scenario"] for r in results]
    ach_ids = [r["ach_id"] for r in results]

    assert all("Protein" in s for s in scenarios)
    assert "B" not in ach_ids and "D" not in ach_ids
    assert set(ach_ids) == {"A", "C"}


def test_require_protein_default_is_off(fake_data):
    """P6 Step 4: require_protein defaults to False — RNA-only cells kept."""
    fake_data["protein"] = {"A": 50.0}   # only A has protein

    results = run_ranking([HIGH], top_n=10)
    scenarios = [r["scenario"] for r in results]
    assert "RNA only" in scenarios   # B, C, D still in results


# Step 5: LOW off tier mixed-coverage 

def test_low_off_tier_is_tied_under_mixed_coverage(fake_data):
    """
    P6 Step 5: the known LOW-search bug.

    In an off tier where some cells have protein and some don't, the cells
    without protein previously got combined_score = 0 while cells with low
    protein scored > 0 — so RNA-only cells unfairly ranked first. The fix
    (apply_low_expression_to_protein, Case B) ties all off cells.

    Prab's Case A (all off cells have protein) is preserved by a separate
    test in test_low_expression.py.
    """
    # All four cells are off (TPM < 1 in both sources).
    fake_data["expression"] = {
        "fact_expression_depmap": {"A": 0.1, "B": 0.1, "C": 0.1, "D": 0.1},
        "fact_expression_hpa":    {"A": 0.1, "B": 0.1, "C": 0.1, "D": 0.1},
    }
    # Only B and D have protein — mixed coverage.
    fake_data["protein"] = {"B": 5.0, "D": 2.0}

    results = run_ranking([LOW], low_expression_mode="off", top_n=10)
    top4 = results[:4]

    # All four scores should be equal (tied in the off tier).
    scores = [r["score"] for r in top4]
    assert len(set(scores)) == 1, f"expected tied scores, got {scores}"

    # Tiebreak by ach_id, so A, B, C, D in alphabetical order.
    assert [r["ach_id"] for r in top4] == ["A", "B", "C", "D"]
