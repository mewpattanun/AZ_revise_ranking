import styles from "./SensitivityModal.module.css";

/**
 * Weight-sensitivity modal (AZ priority 7).
 *
 * Visualises how the ranking changes across 5 RNA:Protein weight settings
 * (RNA only, 90:10, 70:30 default, 50:50, 30:70). Shows:
 *   1. A top-10 overlap heatmap vs the 70:30 default.
 *   2. A rank-correlation row vs the default.
 *   3. A list of weight-sensitive candidates (default rank <= 50 whose
 *      rank changes by >= 10 positions across settings), with a mini
 *      table of rank per setting.
 *
 * The ranking algorithm and the default weight are unchanged — this is a
 * pure explainability layer.
 */

const SETTING_LABELS = {
  rna_only: "RNA only",
  "90_10":  "90 : 10",
  "70_30":  "70 : 30",
  "50_50":  "50 : 50",
  "30_70":  "30 : 70",
};

const HEATMAP_ORDER = ["rna_only", "90_10", "70_30", "50_50", "30_70"];

// Colour thresholds for overlap and correlation (same green/yellow/red bands
// used in the AZ slide 19 heatmap mockup).
function overlapClass(v) {
  if (v >= 0.85) return styles.good;
  if (v >= 0.60) return styles.warn;
  return styles.bad;
}

function correlationClass(v) {
  if (v >= 0.90) return styles.good;
  if (v >= 0.70) return styles.warn;
  return styles.bad;
}


export default function SensitivityModal({ data, loading, onClose }) {
  return (
    <div className={styles.overlay} onClick={onClose}>
      <div className={styles.modal} onClick={(e) => e.stopPropagation()}>
        {/* Header */}
        <div className={styles.header}>
          <h2 className={styles.title}>Weight Sensitivity Analysis</h2>
          <button className={styles.closeBtn} onClick={onClose}>×</button>
        </div>

        <div className={styles.body}>
          {loading && (
            <div className={styles.loading}>
              <div className={styles.spinner} />
              <p>Running the sweep across 5 weight settings… (~5–10 s)</p>
            </div>
          )}

          {!loading && !data && (
            <div className={styles.empty}>
              No sensitivity data yet. Run a search first.
            </div>
          )}

          {!loading && data && (
            <>
              {/* Intro */}
              <p className={styles.intro}>
                The ranking is computed at <strong>5 RNA : Protein weight
                settings</strong> and compared to the default (70 : 30). High
                overlap and correlation mean the candidate list is robust to
                the weight choice.
              </p>

              {/* Heatmap */}
              <h3 className={styles.sectionTitle}>Top-10 Overlap vs Default</h3>
              <div className={styles.heatmap}>
                <div className={styles.heatmapRow}>
                  <div className={styles.heatmapLabel}>Overlap</div>
                  {HEATMAP_ORDER.map((key) => {
                    const v = data.top_10_overlap?.[key] ?? 0;
                    return (
                      <div key={key} className={styles.heatmapCell}>
                        <div className={styles.heatmapSettingLabel}>
                          {SETTING_LABELS[key]}
                          {key === data.default_setting && (
                            <span className={styles.defaultTag}> (default)</span>
                          )}
                        </div>
                        <div className={`${styles.heatmapValue} ${overlapClass(v)}`}>
                          {(v * 100).toFixed(0)}%
                        </div>
                      </div>
                    );
                  })}
                </div>

                <div className={styles.heatmapRow}>
                  <div className={styles.heatmapLabel}>Rank corr.</div>
                  {HEATMAP_ORDER.map((key) => {
                    const v = data.rank_correlation?.[key] ?? 0;
                    return (
                      <div key={key} className={styles.heatmapCell}>
                        <div className={styles.heatmapSettingLabel}>&nbsp;</div>
                        <div className={`${styles.heatmapValue} ${correlationClass(v)}`}>
                          {v.toFixed(2)}
                        </div>
                      </div>
                    );
                  })}
                </div>
              </div>

              <p className={styles.legend}>
                <span className={`${styles.swatch} ${styles.good}`} /> stable ≥ 85 % / ρ ≥ 0.90 &nbsp;
                <span className={`${styles.swatch} ${styles.warn}`} /> caution &nbsp;
                <span className={`${styles.swatch} ${styles.bad}`} /> sensitive
              </p>

              {/* Sensitive candidates */}
              <h3 className={styles.sectionTitle}>
                Weight-Sensitive Candidates ({data.sensitive_candidates?.length ?? 0})
              </h3>
              {(!data.sensitive_candidates || data.sensitive_candidates.length === 0) ? (
                <p className={styles.noSensitive}>
                  No candidate in the default top-50 changes rank by more
                  than 10 positions across the sweep — the ranking is stable
                  for this query.
                </p>
              ) : (
                <div className={styles.tableWrap}>
                  <table className={styles.candidatesTable}>
                    <thead>
                      <tr>
                        <th>Cell line (ACH-ID)</th>
                        <th>Default rank</th>
                        <th>Max Δ rank</th>
                        {HEATMAP_ORDER.map((key) => (
                          <th key={key}>{SETTING_LABELS[key]}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {data.sensitive_candidates.map((c) => (
                        <tr key={c.ach_id}>
                          <td className={styles.mono}>{c.ach_id}</td>
                          <td>#{c.default_rank}</td>
                          <td className={styles.warnText}>±{c.max_rank_change}</td>
                          {HEATMAP_ORDER.map((key) => (
                            <td key={key}>
                              {c.ranks_per_setting?.[key] != null
                                ? `#${c.ranks_per_setting[key]}`
                                : "—"}
                            </td>
                          ))}
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {/* Footer note */}
              <p className={styles.footerNote}>
                Weight selection is a methodological decision. The ranking
                algorithm and the default weight are unchanged — this view
                exists so the dependence is visible and the default is
                defensible rather than hidden.
              </p>
            </>
          )}
        </div>
      </div>
    </div>
  );
}
