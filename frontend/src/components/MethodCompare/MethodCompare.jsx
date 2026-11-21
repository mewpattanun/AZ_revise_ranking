import { useEffect, useMemo, useState } from "react";
import { rankCellLines } from "../../api";
import {
  RNA_METHODS,
  LIST_LIMIT,
  buildComparison,
  comparisonCsv,
  formatRnaValue,
  methodLabel,
  rankText,
} from "../../compareMethods";
import styles from "./MethodCompare.module.css";

function formatTpm(tpm) {
  if (tpm == null) return "-";
  return tpm < 10 ? tpm.toFixed(2) : tpm.toFixed(1);
}

// Better (+) or worse (-) than the robust z rank, e.g. "up 3"
function RankChange({ rank, baseline }) {
  if (rank == null || baseline == null || rank === baseline) return null;
  const up = rank < baseline;
  return (
    <span className={up ? styles.up : styles.down} title="Compared with the robust z rank">
      {up ? "▲" : "▼"}
      {Math.abs(rank - baseline)}
    </span>
  );
}

export default function MethodCompare({ genes, params, onClose }) {
  const [byMethod, setByMethod] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const k = Math.min(params.top_n || 20, LIST_LIMIT);

  useEffect(() => {
    let cancelled = false;
    setLoading(true);
    setError(null);
    setByMethod(null);
    const request = (rna_method) => ({ ...params, rna_method, top_n: LIST_LIMIT });
    Promise.all(
      RNA_METHODS.map((m) => rankCellLines(request(m.value)).then((res) => [m.value, res.results]))
    )
      .then((entries) => {
        if (!cancelled) setByMethod(Object.fromEntries(entries));
      })
      .catch((err) => {
        if (!cancelled) setError(err.message);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, [params]);

  useEffect(() => {
    const onKey = (e) => e.key === "Escape" && onClose();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [onClose]);

  const comparison = useMemo(
    () => (byMethod ? buildComparison(byMethod, k, genes.map((g) => g.hugo)) : null),
    [byMethod, k, genes]
  );
  const singleGene = genes.length === 1;

  function downloadCsv() {
    const blob = new Blob([comparisonCsv(comparison, genes)], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `rna_method_comparison_${genes.map((g) => g.hugo).join("_")}.csv`;
    a.click();
    URL.revokeObjectURL(url);
  }

  return (
    <div className={styles.overlay} onClick={onClose}>
      <div className={styles.modal} onClick={(e) => e.stopPropagation()}>
        <div className={styles.header}>
          <div>
            <h2 className={styles.title}>Compare RNA normalisation methods</h2>
            <div className={styles.subtitle}>
              {genes.map((g) => `${g.hugo} (${g.direction})`).join(", ")} · top {k} of each method
            </div>
          </div>
          <button className={styles.closeBtn} onClick={onClose} aria-label="Close">
            ×
          </button>
        </div>

        <div className={styles.body}>
          {loading && <div className={styles.status}>Ranking with the 3 methods...</div>}
          {error && <div className={styles.error}>{error}</div>}

          {comparison && (
            <>
              <div className={styles.overlaps}>
                {comparison.overlaps.map((o) => (
                  <div key={`${o.a}-${o.b}`} className={styles.overlapCard}>
                    <div className={styles.overlapNames}>
                      {methodLabel(o.a)} vs {methodLabel(o.b)}
                    </div>
                    <div className={styles.overlapValue}>
                      {o.shared} <span>of {o.of} shared in the top {k}</span>
                    </div>
                  </div>
                ))}
              </div>

              <div className={styles.tableWrap}>
                <table className={styles.table}>
                  <thead>
                    <tr>
                      <th rowSpan={2} className={styles.left}>Cell line</th>
                      <th rowSpan={2} className={styles.left}>Disease</th>
                      {singleGene && <th rowSpan={2}>TPM<br /><span className={styles.thSub}>mean of sources</span></th>}
                      {RNA_METHODS.map((m) => (
                        <th key={m.value} colSpan={3} className={styles.group} style={{ color: m.color, background: m.bg }}>
                          {m.label}
                        </th>
                      ))}
                    </tr>
                    <tr>
                      {RNA_METHODS.map((m) => [
                        <th key={`${m.value}-r`} className={styles.sub}>Rank</th>,
                        <th key={`${m.value}-v`} className={styles.sub}>{m.valueLabel}</th>,
                        <th key={`${m.value}-s`} className={styles.sub}>Score</th>,
                      ])}
                    </tr>
                  </thead>
                  <tbody>
                    {comparison.rows.map((r) => (
                      <tr key={r.id}>
                        <td className={styles.left}>
                          <div className={styles.name}>{r.name}</div>
                          <div className={styles.achId}>{r.id}</div>
                        </td>
                        <td className={`${styles.left} ${styles.disease}`}>{r.disease || "-"}</td>
                        {singleGene && <td className={styles.num}>{formatTpm(r.tpm)}</td>}
                        {RNA_METHODS.map((m) => {
                          const rank = r.ranks[m.value];
                          const inTop = rank != null && rank <= k;
                          return [
                            <td key={`${m.value}-r`} className={styles.rankCell}>
                              <span
                                className={`${styles.rankPill} ${inTop ? "" : styles.outside}`}
                                style={inTop ? { color: m.color, background: m.bg } : undefined}
                              >
                                {rankText(rank, comparison.complete[m.value])}
                              </span>
                              {m.value !== "robust_z" && (
                                <RankChange rank={rank} baseline={r.ranks.robust_z} />
                              )}
                            </td>,
                            <td key={`${m.value}-v`} className={styles.num}>
                              {r.values[m.value].map((v) => (
                                <div key={v.gene}>
                                  {singleGene ? "" : <span className={styles.geneTag}>{v.gene}</span>}
                                  {formatRnaValue(v.value)}
                                </div>
                              ))}
                            </td>,
                            <td key={`${m.value}-s`} className={styles.num}>
                              {r.scores[m.value] == null ? "-" : r.scores[m.value].toFixed(3)}
                            </td>,
                          ];
                        })}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              <div className={styles.footer}>
                <div className={styles.notes}>
                  <p>
                    <b>What changes:</b> only how each RNA source is put on a common scale before DepMap and HPA are
                    averaged. Protein, weights and filters are the same in all three. The middle column of each
                    method is that method's RNA value (z, or log2 for the last one); Score is the 0 to 1 score made from it.
                  </p>
                  <ul>
                    {RNA_METHODS.map((m) => (
                      <li key={m.value}>
                        <b style={{ color: m.color }}>{m.label}</b>: {m.formula}
                      </li>
                    ))}
                  </ul>
                  <p>
                    Both z-scores take their reference (median and MAD, or mean and SD) from cell lines at or above 1
                    TPM only; lines below 1 TPM still get a score from their own value. log2(TPM + 1) has no reference.
                  </p>
                  <p>
                    Rank 1 = best for the requested direction. Grey = outside the top {k} of that method. "&gt; {LIST_LIMIT}"
                    = beyond the {LIST_LIMIT} cell lines compared. Arrows = rank change vs robust z.
                  </p>
                </div>
                <button className={styles.csvBtn} onClick={downloadCsv}>
                  Download CSV
                </button>
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
}
