"""
Core ranking endpoint — POST /api/rank
"""

from fastapi import APIRouter, HTTPException

from app.config import EXPR_DETECTION_TPM, EXPR_MODERATE_TPM, EXPR_OFF_TPM, SCORING_VERSION
from app.models import RankRequest, RankResponse
from app.services import data_service, ranking_service
from app.services import sensitivity_service

router = APIRouter(prefix="/api", tags=["ranking"])


@router.post("/rank", response_model=RankResponse)
def rank_cell_lines(req: RankRequest):
    """Run the full ranking pipeline and return ranked cell lines."""
    # Resolve each HUGO symbol → Ensembl ID
    resolved_genes = []
    for gt in req.genes:
        gene = data_service.resolve_gene(gt.hugo)
        if gene is None:
            raise HTTPException(status_code=404, detail=f"Gene '{gt.hugo}' not found")
        resolved_genes.append({
            "hugo": gene["hugo_symbol"],
            "ensembl_id": gene["ensembl_id"],
            "direction": gt.direction,
            "weight": gt.weight,
        })

    # GEO is evidence only. Old clients may still send "geo"; it is dropped
    # (with a warning) so the echoed sources match what was actually scored.
    warnings = []
    if "geo" in req.sources:
        warnings.append("GEO is not used for scoring; it is shown as per-study evidence in the cell line detail only.")
    scoring_sources = [s for s in req.sources if s in ranking_service.SCORING_SOURCE_NAMES] \
        or ranking_service.SCORING_SOURCE_NAMES

    notes: list[str] = []
    results = ranking_service.run_ranking(
        genes=resolved_genes,
        w_rna=req.w_rna,
        w_protein=req.w_protein,
        mutation_mode=req.mutation_mode,
        fusion_mode=req.fusion_mode,
        msi_max=req.msi_max,
        cin_max=req.cin_max,
        exclude_metabolite=req.exclude_metabolite,
        metabolite_threshold=req.metabolite_threshold,
        exclude_mirna=req.exclude_mirna,
        mirna_threshold=req.mirna_threshold,
        disease_filter=req.disease_filter,
        lineage_filter=req.lineage_filter,
        subtype_filter=req.subtype_filter,
        target_disease=req.target_disease,
        target_lineage=req.target_lineage,
        target_subtype=req.target_subtype,
        q6_boost_enabled=req.q6_boost_enabled,
        core_only=req.core_only,
        require_protein=req.require_protein,
        top_n=req.top_n,
        scoring_method=req.scoring_method,
        sources=scoring_sources,
        assay_type=req.assay_type,
        multi_gene_method=req.multi_gene_method,
        low_expression_mode=req.low_expression_mode,
        rna_method=req.rna_method,
        notes_out=notes,
    )

    return RankResponse(
        query={
            "scoring_version": SCORING_VERSION,
            "genes": [{"hugo": g["hugo"], "direction": g["direction"], "weight": g["weight"]} for g in resolved_genes],
            "w_rna": req.w_rna,
            "w_protein": req.w_protein,
            "mutation_mode": req.mutation_mode,
            "fusion_mode": req.fusion_mode,
            "msi_max": req.msi_max,
            "cin_max": req.cin_max,
            "exclude_metabolite": req.exclude_metabolite,
            "metabolite_threshold": req.metabolite_threshold,
            "exclude_mirna": req.exclude_mirna,
            "mirna_threshold": req.mirna_threshold,
            "disease_filter": req.disease_filter,
            "lineage_filter": req.lineage_filter,
            "subtype_filter": req.subtype_filter,
            "target_disease": req.target_disease,
            "target_lineage": req.target_lineage,
            "target_subtype": req.target_subtype,
            "q6_boost_enabled": req.q6_boost_enabled,
            "core_only": req.core_only,
            "require_protein": req.require_protein,
            "top_n": req.top_n,
            "scoring_method": req.scoring_method,
            "multi_gene_method": req.multi_gene_method,
            "rna_method": req.rna_method,
            "low_expression": {
                "mode": req.low_expression_mode,
                "detection_tpm": EXPR_DETECTION_TPM,
                "off_tpm": EXPR_OFF_TPM,
                "moderate_tpm": EXPR_MODERATE_TPM,
            },
            "sources": scoring_sources,
            "assay_type": req.assay_type,
        },
        total_results=len(results),
        results=results,
        warnings=warnings,
        notes=notes,
    )




@router.post("/sensitivity")
def sensitivity_analysis(req: RankRequest):
    """
    Weight-sensitivity report for the same request as /api/rank.

    Runs the ranking at 5 RNA:Protein weight settings (RNA-only,
    90:10, 70:30, 50:50, 30:70) and returns top-K overlap, rank
    correlation, and a list of weight-sensitive candidates against
    the 70:30 default. Does not change the ranking algorithm or the
    default weight — AZ priority 7 is about visibility.
    """
    # 1. Resolve HUGO -> Ensembl (same as /api/rank)
    resolved_genes = []
    for gt in req.genes:
        gene = data_service.resolve_gene(gt.hugo)
        if gene is None:
            raise HTTPException(status_code=404, detail=f"Gene '{gt.hugo}' not found")
        resolved_genes.append({
            "hugo": gene["hugo_symbol"],
            "ensembl_id": gene["ensembl_id"],
            "direction": gt.direction,
            "weight": gt.weight,
        })

    # 2. GEO is evidence only (same handling as /api/rank).
    scoring_sources = [s for s in req.sources if s in ranking_service.SCORING_SOURCE_NAMES] \
        or ranking_service.SCORING_SOURCE_NAMES

    # 3. Build kwargs for run_ranking — identical to /api/rank except
    # w_rna / w_protein, which the sensitivity sweep overrides per setting.
    base_kwargs = {
        "mutation_mode": req.mutation_mode,
        "fusion_mode": req.fusion_mode,
        "msi_max": req.msi_max,
        "cin_max": req.cin_max,
        "exclude_metabolite": req.exclude_metabolite,
        "metabolite_threshold": req.metabolite_threshold,
        "exclude_mirna": req.exclude_mirna,
        "mirna_threshold": req.mirna_threshold,
        "disease_filter": req.disease_filter,
        "lineage_filter": req.lineage_filter,
        "subtype_filter": req.subtype_filter,
        "target_disease": req.target_disease,
        "target_lineage": req.target_lineage,
        "target_subtype": req.target_subtype,
        "q6_boost_enabled": req.q6_boost_enabled,
        "core_only": req.core_only,
        "require_protein": req.require_protein,
        "top_n": req.top_n,
        "scoring_method": req.scoring_method,
        "sources": scoring_sources,
        "assay_type": req.assay_type,
        "multi_gene_method": req.multi_gene_method,
        "low_expression_mode": req.low_expression_mode,
    }

    # 4. Run sensitivity sweep
    report = sensitivity_service.compute_sensitivity(resolved_genes, base_kwargs)

    # 5. Echo back the query context so the UI can display what was tested
    report["query"] = {
        "scoring_version": SCORING_VERSION,
        "genes": [{"hugo": g["hugo"], "direction": g["direction"], "weight": g["weight"]} for g in resolved_genes],
    }
    return report