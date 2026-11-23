import { Fragment, useState, useEffect } from "react";
import { getCellLineEvidence, getSimilarCellLines } from "../../api";
import styles from "./DetailPanel.module.css";
import { categoryInfo, formatTpm, expressionTooltip } from "../../expression";
import { methodLabel as rnaMethodLabel } from "../../compareMethods";

export default function DetailPanel({ achId, genes, resultRow, scoringMethod = "robust_z", rnaMethod = "robust_z", maxScore = 0, wRna = 0.7, wProtein = 0.3, mutationMode = "ignore", fusionMode = "ignore", onClose }) {
  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(true);
  const [similar, setSimilar] = useState(null);
  const [similarLoading, setSimilarLoading] = useState(false);
  const [similarOpen, setSimilarOpen] = useState(false);

  useEffect(() => {
    if (!achId || genes.length === 0) return;
    setLoading(true);
    setSimilar(null);
    setSimilarOpen(false);
    getCellLineEvidence(achId, genes, wRna, wProtein)
      .then(setData)
      .catch((err) => {
        console.error("Evidence fetch error:", err);
        setData(null);
      })
      .finally(() => setLoading(false));
  }, [achId, genes, wRna, wProtein]);

  function loadSimilar() {
    if (similarLoading || similar) return;
    setSimilarOpen(true);
    setSimilarLoading(true);
    getSimilarCellLines(achId, 5)
      .then(setSimilar)
      .catch((err) => {
        console.error("Similar cell lines fetch error:", err);
        setSimilar(null);
      })
      .finally(() => setSimilarLoading(false));
  }

  if (!achId) return null;

  const score = resultRow?.score;
  const confidence = resultRow?.confidence;
  const scenario = resultRow?.scenario;

  // Colour the score using RAW score value (same thresholds as the results table).
  // Raw score gives a colour that reflects absolute quality, so it doesn't
  // mislead when the top result in a filtered list happens to have a low score.
  const scoreColor =
    score == null ? "#e05a4d"
    : score >= 0.7 ? "#2f9e6f"
    : score >= 0.4 ? "#f0b429"
    : "#e05a4d";

  // Per-source combined scores. For single gene: { DepMap: 0.89, HPA: 0.85 }
  // For multi-gene: { EGFR: { DepMap: 0.89 }, SEC61G: { DepMap: 0.92 } }
  const rawSourceScores = data?.source_scores || {};
  // Detect if nested (multi-gene) or flat (single-gene)
  const isNestedScores = genes.length > 1 && typeof Object.values(rawSourceScores)[0] === "object";

  // Per-gene combined scores for multi-gene searches: { EGFR: 0.98, SEC61G: 0.95 }
  const perGeneScores = data?.per_gene_scores || {};

  return (
    <div className={styles.overlay} onClick={onClose}>
      <div className={styles.panel} onClick={(e) => e.stopPropagation()}>
        {/* ── Header ── */}
        <div className={styles.header}>
          <div className={styles.headerInfo}>
            <h2 className={styles.title}>
              {data?.cell_line_name || resultRow?.cell_line_name || achId}
              <span className={styles.subtitle}>
                {" — Multi-Omics Evidence Breakdown"}
              </span>
            </h2>
            <span className={styles.achMeta}>
              {achId}
              {data?.primary_disease && ` · ${data.primary_disease}`}
              {data?.lineage && ` · ${data.lineage}`}
              {data?.sex && ` · ${data.sex}`}
              {data?.growth_pattern && ` · ${data.growth_pattern}`}
            </span>
          </div>
          <div className={styles.headerRight}>
            {score != null && (
              <div className={styles.scoreStack}>
                <div className={styles.detailScore} style={{ background: scoreColor }}>{score.toFixed(2)}</div>
                {resultRow?.modality_coverage != null && (
                  <div className={styles.coverageLabel}>
                    Coverage: {resultRow.modality_coverage.toFixed(2)}
                  </div>
                )}
              </div>
            )}
            <button className={styles.closeBtn} onClick={onClose}>×</button>
          </div>
        </div>

        {loading ? (
          <div className={styles.loading}>
            <div className={styles.spinner} />
            Computing evidence breakdown...
          </div>
        ) : !data ? (
          <div className={styles.loading}>No evidence data available</div>
        ) : (
          <div className={styles.body}>
            {rnaMethod !== "robust_z" && (
              <div className={styles.methodNote}>
                The results table uses <b>{rnaMethodLabel(rnaMethod)}</b>. The z-scores and ranks below are always robust z.
              </div>
            )}
            {/* ── Per-gene expression cards ── */}
            <div className={styles.geneCards}>
              {genes.map((g) => {
                const gd = data.genes?.[g.hugo];
                if (!gd) return null;
                const dirLabel = gd.direction === "high" ? "HIGH ↑" : "LOW ↓";
                const dirClass = gd.direction === "high" ? styles.dirHigh : styles.dirLow;
                return (
                  <div key={g.hugo} className={styles.geneCard}>
                    <div className={styles.geneCardHeader}>
                      <span className={styles.geneCardLabel}>
                        {g.hugo} Expression
                      </span>
                      <span className={`${styles.dirBadge} ${dirClass}`}>
                        {dirLabel}
                      </span>
                    </div>
                    {perGeneScores[g.hugo] != null && (
                      <div className={styles.geneCardScore}>
                        Score: <strong>{perGeneScores[g.hugo].toFixed(4)}</strong>
                      </div>
                    )}
                    {gd.expression_rank == null && (
                      <div className={styles.geneCardRank}>No expression data</div>
                    )}
                    {resultRow?.expression?.[g.hugo] && (() => {
                      const info = resultRow.expression[g.hugo];
                      const cat = categoryInfo(info.category);
                      return (
                        <div className={styles.geneCardSub} title={expressionTooltip(info)}>
                          Expression: <strong>{formatTpm(info)}</strong>{" "}
                          <span className={styles.exprPill} style={{ color: cat.color, background: cat.bg }}>
                            {cat.label}
                          </span>
                          {info.sources_disagree && (
                            <span style={{ color: "#856404" }}> ⚠ DepMap and HPA disagree</span>
                          )}
                        </div>
                      );
                    })()}
                    <div className={styles.geneCardSub}>
                      {scoringMethod === "robust_z"
                        ? `Score from the size of one RNA consensus robust z (${gd.source_count} source${gd.source_count !== 1 ? "s" : ""})`
                        : scoringMethod === "rrf"
                        ? `RRF ensemble (Z-Score + Percentile) on one RNA consensus from ${gd.source_count} source${gd.source_count !== 1 ? "s" : ""}`
                        : scoringMethod === "zscore"
                        ? `Z-Score ranking on one RNA consensus from ${gd.source_count} source${gd.source_count !== 1 ? "s" : ""}`
                        : `Percentile ranking on one RNA consensus from ${gd.source_count} source${gd.source_count !== 1 ? "s" : ""}`}
                    </div>
                    <div className={styles.geneCardSub}>
                      Sources:{" "}
                      {["DepMap", "HPA"].map((src) => {
                        const has = gd.sources.some((s) => s.source === src);
                        return (
                          <span key={src} className={has ? styles.srcYes : styles.srcNo}>
                            {src} {has ? "✓" : "✗"}
                          </span>
                        );
                      })}
                    </div>
                    {gd.rna_consensus && (
                      <div className={styles.geneCardSub}>
                        RNA consensus (robust z): <strong>{gd.rna_consensus.z.toFixed(2)}</strong>
                        {" · "}
                        Discordance (QC):{" "}
                        <strong>
                          {gd.rna_consensus.discordance != null
                            ? gd.rna_consensus.discordance.toFixed(2)
                            : "n/a (one source)"}
                        </strong>
                      </div>
                    )}
                    {gd.rna_consensus?.mad_zero_sources?.length > 0 && (
                      <div className={styles.geneCardSub} style={{ color: "#856404" }}>
                        ⚠ MAD = 0 in {gd.rna_consensus.mad_zero_sources.join(" + ")}: most cell lines at or above
                        1 TPM have the same value, so the standard deviation was used for the z-score.
                      </div>
                    )}
                  </div>
                );
              })}

              {/* Protein summary card — always shown */}
              <div className={styles.geneCard}>
                <div className={styles.geneCardHeader}>
                  <span className={styles.geneCardLabel}>Protein Detection</span>
                </div>
                {data.protein ? (
                  <>
                    <div className={styles.geneCardRank}>
                      Rank #{data.protein.z_rank}{" "}
                      <span className={styles.rankTotal}>/ {data.protein.total?.toLocaleString()}</span>
                    </div>
                    <div className={styles.geneCardSub}>
                      {genes[0]?.hugo}: Intensity {data.protein.intensity} · Gygi proteomics
                    </div>
                    <div className={`${styles.geneCardSub} ${styles.proteinConfirm}`}>
                      Confirmed at protein level ✓
                    </div>
                  </>
                ) : (
                  <>
                    <div className={styles.geneCardRank}>Not measured</div>
                    <div className={styles.geneCardSub}>
                      {genes[0]?.hugo}: not profiled in Gygi proteomics
                      <br />
                      <span style={{ color: "#999", fontSize: "11px" }}>
                        "Not measured" ≠ "measured as zero" — no information available
                      </span>
                    </div>
                  </>
                )}
              </div>

            </div>

            {/* ── Source-Level Scores per gene ── */}
            {genes.map((g) => {
              const gd = data.genes?.[g.hugo];
              if (!gd || (gd.sources.length === 0 && !gd.geo)) return null;
              return (
                <div key={g.hugo} className={styles.sourceSection}>
                  <h4 className={styles.sectionLabel}>
                    Source-Level Scores — {g.hugo}
                  </h4>
                  <div className={styles.sourceGrid}>
                    {gd.sources.map((src) => {
                      const combined = isNestedScores
                        ? rawSourceScores[g.hugo]?.[src.source]
                        : rawSourceScores[src.source];
                      return (
                        <div key={src.source} className={styles.sourceCard}>
                          <h5 className={styles.sourceTitle}>{src.source}</h5>
                          <div className={styles.sourceRow}>
                            <span>TPM</span>
                            <span className={styles.val}>{src.tpm}</span>
                          </div>
                          {scoringMethod === "robust_z" && src.robust_z != null && (
                            <div className={styles.sourceRow}>
                              <span>Robust Z</span>
                              <span className={styles.val}>
                                {src.robust_z > 0 ? "+" : ""}{src.robust_z}
                              </span>
                            </div>
                          )}
                          {(scoringMethod === "rrf" || scoringMethod === "zscore") && (
                            <div className={styles.sourceRow}>
                              <span>Z-Score</span>
                              <span className={styles.val}>
                                {src.z_score > 0 ? "+" : ""}{src.z_score}
                              </span>
                            </div>
                          )}
                          {(scoringMethod === "rrf" || scoringMethod === "percentile") && (
                            <div className={styles.sourceRow}>
                              <span>Percentile</span>
                              <span className={styles.val}>{src.percentile}%</span>
                            </div>
                          )}
                          {(scoringMethod === "rrf" || scoringMethod === "zscore") && (
                            <div className={styles.sourceRow}>
                              <span>Z-Score Rank</span>
                              <span className={styles.val}>#{src.z_rank}</span>
                            </div>
                          )}
                          {(scoringMethod === "rrf" || scoringMethod === "percentile") && (
                            <div className={styles.sourceRow}>
                              <span>Percentile Rank</span>
                              <span className={styles.val}>#{src.pct_rank}</span>
                            </div>
                          )}
                          {combined != null && (
                            <div className={`${styles.sourceRow} ${styles.sourceScoreRow}`}>
                              <span>Score</span>
                              <span className={styles.val}>{combined.toFixed(4)}</span>
                            </div>
                          )}
                        </div>
                      );
                    })}

                    {/* Protein source card — always shown */}
                    <div className={styles.sourceCard}>
                      <h5 className={styles.sourceTitle}>Protein (Gygi)</h5>
                      {data.protein ? (
                        <>
                          <div className={styles.sourceRow}>
                            <span>Intensity</span>
                            <span className={styles.val}>{data.protein.intensity}</span>
                          </div>
                          {scoringMethod === "robust_z" && data.protein.robust_z != null && (
                            <div className={styles.sourceRow}>
                              <span>Robust Z</span>
                              <span className={styles.val}>
                                {data.protein.robust_z > 0 ? "+" : ""}{data.protein.robust_z}
                              </span>
                            </div>
                          )}
                          {(scoringMethod === "rrf" || scoringMethod === "zscore") && (
                            <div className={styles.sourceRow}>
                              <span>Z-Score</span>
                              <span className={styles.val}>
                                {data.protein.z_score > 0 ? "+" : ""}{data.protein.z_score}
                              </span>
                            </div>
                          )}
                          {(scoringMethod === "rrf" || scoringMethod === "percentile") && (
                            <div className={styles.sourceRow}>
                              <span>Percentile</span>
                              <span className={styles.val}>{data.protein.percentile}%</span>
                            </div>
                          )}
                          {(scoringMethod === "rrf" || scoringMethod === "zscore") && (
                            <div className={styles.sourceRow}>
                              <span>Z-Score Rank</span>
                              <span className={styles.val}>#{data.protein.z_rank}</span>
                            </div>
                          )}
                          {(scoringMethod === "rrf" || scoringMethod === "percentile") && (
                            <div className={styles.sourceRow}>
                              <span>Percentile Rank</span>
                              <span className={styles.val}>#{data.protein.pct_rank}</span>
                            </div>
                          )}
                          {scoringMethod === "robust_z" && data.protein.mad_zero && (
                            <div className={styles.sourceRow} style={{ color: "#856404" }}>
                              <span>⚠ MAD = 0: SD used for the z-score</span>
                            </div>
                          )}
                        </>
                      ) : (
                        <div className={styles.sourceRow}>
                          <span>No protein data</span>
                          <span className={styles.val}>—</span>
                        </div>
                      )}
                    </div>

                    {/* GEO support card: extra evidence only, not part of the score or rank.
                        Only present when the GSM -> study map has been built. */}
                    {gd.geo && (
                      <div className={`${styles.sourceCard} ${styles.geoCard}`}>
                        <h5 className={styles.sourceTitle}>GEO support</h5>
                        <div className={styles.sourceRow}>
                          <span>Status</span>
                          <span
                            className={styles.val}
                            style={{
                              color:
                                gd.geo.status === "supportive" ? "#2d6a4f"
                                : gd.geo.status === "contradictory" ? "#b02a37"
                                : gd.geo.status === "mixed" ? "#856404"
                                : "#999",
                            }}
                          >
                            {gd.geo.status === "unavailable"
                              ? "Unavailable"
                              : `${gd.geo.status.charAt(0).toUpperCase()}${gd.geo.status.slice(1)} ${gd.geo.supportive}/${gd.geo.total}`}
                          </span>
                        </div>
                        {gd.geo.total === 0 && (
                          <div className={styles.sourceRow}>
                            <span>No usable GEO study</span>
                            <span className={styles.val}>—</span>
                          </div>
                        )}
                        {gd.geo.studies?.length > 0 && (
                          // One aligned row per study: study | samples | z | verdict
                          <div className={styles.geoTable}>
                            <span className={styles.geoHead}>Study</span>
                            <span className={`${styles.geoHead} ${styles.geoNum}`}>Samples</span>
                            <span className={`${styles.geoHead} ${styles.geoNum}`}>Robust Z</span>
                            <span className={`${styles.geoHead} ${styles.geoNum}`}>Verdict</span>
                            {gd.geo.studies.map((s) => {
                              const verdict =
                                s.verdict === "supportive" ? { mark: "✓", label: "Supports", color: "#2d6a4f" }
                                : s.verdict === "contradictory" ? { mark: "✗", label: "Contradicts", color: "#b02a37" }
                                : { mark: "~", label: "Neutral", color: "#856404" };
                              return (
                                <Fragment key={s.gse_id}>
                                  <span className={styles.geoStudy}>{s.gse_id}</span>
                                  <span className={styles.geoNum}>{s.n_samples}</span>
                                  <span className={`${styles.geoNum} ${styles.val}`}>
                                    {s.z > 0 ? "+" : ""}{s.z}
                                  </span>
                                  <span className={styles.geoNum} style={{ color: verdict.color, fontWeight: 600 }}>
                                    {verdict.mark} {verdict.label}
                                  </span>
                                </Fragment>
                              );
                            })}
                          </div>
                        )}
                        <div className={styles.geneCardSub} style={{ marginTop: 6 }}>
                          Not used in the score or rank.
                        </div>
                      </div>
                    )}
                  </div>
                </div>
              );
            })}

            {/* ── Mutation detail (only when mutation filter is active) ── */}
            {mutationMode === "include" && data.mutations && data.mutations.length > 0 && (
              <div className={styles.sourceSection}>
                <h4 className={styles.sectionLabel}>
                  Mutations in {genes.map((g) => g.hugo).join(", ")}
                </h4>
                <div className={styles.mutTableWrap}>
                  <table className={styles.mutTable}>
                    <thead>
                      <tr>
                        <th>Variant</th>
                        <th>Type</th>
                        <th>Protein</th>
                        <th className={styles.thCenter}>Driver</th>
                        <th className={styles.thCenter}>Hotspot</th>
                        <th className={styles.thCenter}>Damaging</th>
                        <th className={styles.thCenter}>Loss of function</th>
                        <th className={styles.thCenter}>Mutation %</th>
                      </tr>
                    </thead>
                    <tbody>
                      {data.mutations.map((m, i) => (
                        <tr key={i}>
                          <td className={styles.mono}>
                            {m.chrom}:{m.pos} {m.ref}→{m.alt}
                          </td>
                          <td>{m.variant_type || "—"}</td>
                          <td className={styles.proteinChange}>
                            {m.protein_change || "—"}
                          </td>
                          <td className={styles.tdCenter}>
                            <span className={`${styles.pill} ${m.is_driver ? styles.pillPass : styles.pillFail}`}>
                              {m.is_driver ? "Yes" : "No"}
                            </span>
                          </td>
                          <td className={styles.tdCenter}>
                            <span className={`${styles.pill} ${m.is_hotspot ? styles.pillPass : styles.pillFail}`}>
                              {m.is_hotspot ? "Yes" : "No"}
                            </span>
                          </td>
                          <td className={styles.tdCenter}>
                            <span className={`${styles.pill} ${m.is_damaging ? styles.pillPass : styles.pillFail}`}>
                              {m.is_damaging ? "Yes" : "No"}
                            </span>
                          </td>
                          <td className={styles.tdCenter}>
                            <span className={`${styles.pill} ${m.is_lof ? styles.pillPass : styles.pillFail}`}>
                              {m.is_lof ? "Yes" : "No"}
                            </span>
                          </td>
                          <td className={styles.tdCenter}>
                            {m.allele_freq != null
                              ? `${(m.allele_freq * 100).toFixed(1)}%`
                              : "—"}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            )}

            {/* ── Fusions (only when fusion filter is active) ── */}
            {fusionMode === "include" && data.fusions && data.fusions.length > 0 && (
              <div className={styles.sourceSection}>
                <h4 className={styles.sectionLabel}>
                  Fusions ({data.fusions.length})
                </h4>
                <div className={styles.mutTableWrap}>
                  <table className={styles.mutTable}>
                    <thead>
                      <tr>
                        <th>Fusion</th>
                        <th>Confidence</th>
                        <th>Frame</th>
                        <th>Fusion Fragments (per million)</th>
                      </tr>
                    </thead>
                    <tbody>
                      {data.fusions.map((f, i) => (
                        <tr key={i}>
                          <td className={styles.mono}>
                            {f.fusion_name || `${f.gene1_hugo}—${f.gene2_hugo}`}
                          </td>
                          <td>{f.confidence}</td>
                          <td>{f.reading_frame || "—"}</td>
                          <td>{f.ffpm != null ? f.ffpm.toFixed(3) : "—"}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            )}

            {/* ── Similar Alternatives ── */}
            <div className={styles.sourceSection}>
              <div
                className={styles.similarHeader}
                onClick={() => { if (!similarOpen) loadSimilar(); else setSimilarOpen(!similarOpen); }}
                style={{ cursor: "pointer" }}
              >
                <h4 className={styles.sectionLabel} style={{ margin: 0 }}>
                  {similarOpen ? "▾" : "▸"} Similar Alternatives
                </h4>
                <span className={styles.similarHint}>
                  {similarOpen ? "" : "Click to find similar cell lines"}
                </span>
              </div>
              {similarOpen && similarLoading ? (
                <div className={styles.similarLoading}>Finding similar cell lines...</div>
              ) : similarOpen && similar?.similar?.length > 0 ? (
                <div className={styles.similarGrid}>
                  {similar.similar.map((s) => {
                    const simPct = (s.similarity * 100).toFixed(1);
                    const simColor =
                      s.similarity >= 0.8 ? "#2f9e6f"
                      : s.similarity >= 0.5 ? "#f0b429"
                      : "#e05a4d";
                    return (
                      <div key={s.ach_id} className={styles.similarCard}>
                        <div className={styles.similarCardHeader}>
                          <div className={styles.similarCardName}>
                            {s.cell_line_name || s.ach_id}
                            {s.is_derivative && (
                              <span className={styles.derivBadge} title="Shares a parent cell line">
                                Derivative
                              </span>
                            )}
                          </div>
                          <div className={styles.similarScore} style={{ background: simColor }}>
                            {simPct}%
                          </div>
                        </div>
                        <div className={styles.similarMeta}>
                          {s.primary_disease && <span>{s.primary_disease}</span>}
                          {s.lineage && <span> · {s.lineage}</span>}
                        </div>
                        <div className={styles.similarLayers}>
                          {s.expression_similarity != null && (
                            <span className={styles.layerTag}>
                              Expr {(s.expression_similarity * 100).toFixed(0)}%
                            </span>
                          )}
                          {s.protein_similarity != null && (
                            <span className={styles.layerTag}>
                              Prot {(s.protein_similarity * 100).toFixed(0)}%
                            </span>
                          )}
                          {s.metabolomics_similarity != null && (
                            <span className={styles.layerTag}>
                              Metab {(s.metabolomics_similarity * 100).toFixed(0)}%
                            </span>
                          )}
                          {s.mirna_similarity != null && (
                            <span className={styles.layerTag}>
                              miRNA {(s.mirna_similarity * 100).toFixed(0)}%
                            </span>
                          )}
                        </div>
                        {s.scored_on && (
                          <div className={styles.similarScoredOn}>
                            Scored on: {s.scored_on.join(", ")}
                          </div>
                        )}
                      </div>
                    );
                  })}
                </div>
              ) : similarOpen && similar ? (
                <div className={styles.similarLoading}>No similar cell lines found</div>
              ) : null}
            </div>

            {/* ── Scoring method info ── */}
            <div className={styles.methodInfo}>
              <h4 className={styles.methodTitle}>Scoring Method</h4>
              <p className={styles.methodText}>
                {scoringMethod === "robust_z" && (
                  <>Expression scored from the <code>size</code> of one RNA consensus robust z (DepMap + HPA, TPM, median/MAD), compressed with <code>asinh</code> and scaled to 0–1. Rank is the sorted score, so a cell line far above the rest scores clearly higher than one just above the median.</>
                )}
                {scoringMethod === "rrf" && (
                  <>Expression ranked via <code>Z-Score</code> + <code>Percentile</code> on a single RNA consensus (robust z, DepMap + HPA), combined with <code>Reciprocal Rank Fusion</code> (k=60).</>
                )}
                {scoringMethod === "zscore" && (
                  <>Expression ranked via <code>Z-Score</code> only on a single RNA consensus (robust z, DepMap + HPA), combined with <code>Reciprocal Rank Fusion</code> (k=60).</>
                )}
                {scoringMethod === "percentile" && (
                  <>Expression ranked via <code>Percentile</code> only on a single RNA consensus (robust z, DepMap + HPA), combined with <code>Reciprocal Rank Fusion</code> (k=60).</>
                )}
                {data.protein && (
                  <>
                    {scoringMethod === "robust_z"
                      ? " Protein is scored the same way, from the size of its own robust z, "
                      : " Protein scored separately, "}
                    combined with RNA at the selected weighting.
                  </>
                )}
                {genes.length > 1 && (
                  <> Several genes are combined with a weighted geometric mean of the per-gene scores.</>
                )}
                {" "}Per-source scores use the selected RNA/Protein weights on the
                same scale as the overall score.
              </p>
            </div>

            {/* ── Data coverage bar ── */}
            <div className={styles.omicsBar}>
              <div>
                <span className={styles.omicsLabel}>Data Coverage</span>
                <br />
                <span className={styles.omicsText}>
                  RNA-seq ({data.data_coverage.rna_sources}/2 sources)
                  {data.data_coverage.has_protein && " · Proteomics"}
                  {data.data_coverage.has_mutations && " · Mutations"}
                  {data.data_coverage.has_fusions && " · Fusions"}
                </span>
              </div>
              <span className={styles.omicsCount}>
                {data.data_coverage.total_types}/{data.data_coverage.max_types} DATA TYPES
              </span>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
