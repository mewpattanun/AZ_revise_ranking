"""
Pydantic models for request/response validation.
"""

from pydantic import BaseModel, Field, model_validator


# ── Request models ────────────────────────────────────────────

class GeneTarget(BaseModel):
    hugo: str = Field(..., description="HUGO gene symbol, e.g. EGFR")
    direction: str = Field("high", pattern="^(high|low)$", description="high = want high expression, low = want low")
    weight: float = Field(1.0, gt=0, description="Importance of this gene in a multi-gene search (default: equal weights)")


class RankRequest(BaseModel):
    genes: list[GeneTarget] = Field(..., min_length=1, description="Target genes with direction")
    w_rna: float = Field(0.7, ge=0.0, le=1.0, description="Weight for RNA expression")
    w_protein: float = Field(0.3, ge=0.0, le=1.0, description="Weight for protein")
    mutation_mode: str = Field("ignore", pattern="^(include|exclude|ignore)$")
    fusion_mode: str = Field("ignore", pattern="^(include|exclude|ignore)$")
    msi_max: float | None = Field(None, description="Exclude cell lines with MSIScore above this")
    cin_max: float | None = Field(None, description="Exclude cell lines with CIN above this")
    exclude_metabolite: str | None = Field(None, description="Metabolite name to filter on")
    metabolite_threshold: float | None = Field(None, description="Exclude above this level; required if exclude_metabolite is set")
    exclude_mirna: str | None = Field(None, description="miRNA id to filter on, e.g. hsa-miR-21")
    mirna_threshold: float | None = Field(None, description="Exclude above this level; required if exclude_mirna is set")
    disease_filter: str | None = Field(None, description="Filter by primary_disease")
    lineage_filter: str | None = Field(None, description="Filter by lineage")
    subtype_filter: str | None = Field(None, description="Filter by Subtype")
    target_disease: str | None = Field(None, description="Preferred disease for Q6 boost")
    target_lineage: str | None = Field(None, description="Preferred lineage for Q6 boost")
    target_subtype: str | None = Field(None, description="Preferred subtype for Q6 boost")
    q6_boost_enabled: bool = Field(True, description="Enable Q6 soft boost (default on)")
    core_only: bool = Field(False, description="Only cell lines with data in all expression sources")
    require_protein: bool = Field(False, description="P6: only cell lines with protein data (exclude RNA-only)")
    top_n: int = Field(20, ge=1, le=200, description="Number of results to return")
    scoring_method: str = Field(
        "robust_z",
        pattern="^(robust_z|rrf|zscore|percentile)$",
        description="Scoring method: robust_z (score from the size of the RNA consensus robust z, default), "
                    "or the legacy rank-based rrf / zscore / percentile",
    )
    multi_gene_method: str = Field(
        "geometric",
        pattern="^(geometric|summed_rrf)$",
        description="How several genes are combined: geometric = weighted geometric mean (default), summed_rrf = legacy",
    )
    low_expression_mode: str = Field(
        "off",
        pattern="^(off|measurable|score)$",
        description="What 'low' means for LOW genes: off = 'TPM below 1', below 1 TPM in every source first (default); "
                    "measurable = 'TPM greater than 1', 1-3 TPM first; score = legacy, sort by score only",
    )
    off_method: str = Field(
        "fixed",
        pattern="^(fixed|percent_of_reference|mixture|distance_from_on|strict_mixture|custom)$",
        description="How a gene is judged 'off' in a cell line. fixed = below 1 TPM (default). The others can "
                    "also call a line off at 1 TPM or more when it is low for that gene: percent_of_reference "
                    "(< 10% of the gene's reference level), mixture (low group of a 2-group split), "
                    "distance_from_on (> 3 SD below the expressing group), strict_mixture (mixture, only genes "
                    "with a clear on/off split), custom (below off_cutoff_tpm). Below 1 TPM is always off.",
    )
    off_cutoff_tpm: float | None = Field(
        None, ge=1.0, le=100000.0,
        description="Cutoff in TPM for off_method 'custom' (a line below it counts as off); at least 1",
    )
    rna_method: str = Field(
        "robust_z",
        pattern="^(robust_z|standard_z|log2_tpm)$",
        description="How each RNA source is put on a common scale before the sources are averaged: "
                    "robust_z = (x - median) / (1.4826 MAD) (default), standard_z = (x - mean) / SD "
                    "(both use only the cell lines at or above 1 TPM for the reference), "
                    "log2_tpm = log2(TPM + 1) with no z-score",
    )
    sources: list[str] = Field(
        default=["depmap", "hpa"],
        description="RNA expression sources to include in ranking: depmap, hpa. Protein is always scored separately (see w_protein)."
    )
    assay_type: str | None = Field(
        None,
        description="Planned assay type for compatibility warnings: adherent_screen, 3d_spheroid, suspension_screen, flexible"
    )

    @model_validator(mode="after")
    def _some_weight(self):
        if self.w_rna + self.w_protein <= 0:
            raise ValueError("w_rna and w_protein cannot both be 0")
        if self.off_method == "custom" and self.off_cutoff_tpm is None:
            raise ValueError("off_method 'custom' needs off_cutoff_tpm")
        return self


class CompareRequest(BaseModel):
    ach_ids: list[str] = Field(..., min_length=2, description="Cell line IDs to compare")
    gene: str = Field(..., description="Gene to compare across cell lines")


# ── Response models ───────────────────────────────────────────

class GeneInfo(BaseModel):
    hugo_symbol: str
    ensembl_id: str


class GeneSearchResult(BaseModel):
    results: list[GeneInfo]


class RankedCellLine(BaseModel):
    rank: int
    ach_id: str
    cell_line_name: str | None = None
    primary_disease: str | None = None
    lineage: str | None = None
    score: float
    # Share of the RNA/protein weight this cell line has data for (1.0 = both).
    # Replaces the old "confidence" field, which held the same number.
    modality_coverage: float
    scenario: str
    growth_pattern: str | None = None
    is_core: bool = False
    mutations: list[dict] | None = None
    fusions: list[dict] | None = None
    q6_score: float | None = None
    match_level: str | None = None
    base_score: float | None = None
    per_gene_rank: dict[str, int] | None = None
    q7_status: str | None = None
    q7_warning: str | None = None
    score_flags: list[str] | None = None
    geo: dict[str, dict] | None = None
    # Per gene: {"category", "tpm", "tpm_by_source", "sources_disagree"}
    expression: dict[str, dict] | None = None


class RankResponse(BaseModel):
    query: dict
    total_results: int
    results: list[RankedCellLine]
    warnings: list[str] = []
    # Plain-language notes about the result (shown above the results table)
    notes: list[str] = []


class SourceScore(BaseModel):
    source: str
    tpm: float | None = None
    z_score: float | None = None
    percentile: float | None = None
    z_rank: int | None = None
    pct_rank: int | None = None


class CellLineDetail(BaseModel):
    ach_id: str
    cell_line_name: str | None = None
    primary_disease: str | None = None
    lineage: str | None = None
    sex: str | None = None
    growth_pattern: str | None = None
    gene: str
    expression_sources: list[SourceScore]
    protein_intensity: float | None = None
    protein_z_score: float | None = None
    mutations: list[dict]
    fusions: list[dict]
    score: float | None = None
    confidence: float | None = None
    modality_coverage: float | None = None
    scenario: str | None = None


class FilterOptions(BaseModel):
    values: list[str]
    count: int
