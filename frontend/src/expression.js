// Expression categories from raw TPM (AZ priority 5, set in backend config.py).
// A cell line gets the highest band any source puts it in, so it only counts
// as "off" when DepMap and HPA both agree.

export const EXPRESSION_CATEGORIES = {
  below_detection: { label: "Below detection", short: "Off", color: "#c0392b", bg: "#fde8e8", rule: "< 0.1 TPM in every source" },
  negligible: { label: "Negligible", short: "Off", color: "#b7791f", bg: "#fff4dc", rule: "< 1 TPM in every source" },
  low_measurable: { label: "Low measurable", short: "Low", color: "#1d7a6c", bg: "#e3f4f1", rule: "1–3 TPM" },
  moderate: { label: "Moderate", short: "On", color: "#2b6cb0", bg: "#e6f0fa", rule: "> 3 TPM" },
  missing: { label: "No RNA data", short: "—", color: "#888", bg: "#f0efec", rule: "no RNA value" },
};

export const LOW_EXPRESSION_MODES = [
  { key: "off", label: "TPM below 1", hint: "Cell lines below 1 TPM in both DepMap and HPA first, tied on RNA" },
  { key: "measurable", label: "TPM greater than 1", hint: "Lowest cell lines from 1 TPM upwards first (DepMap and HPA agreeing first); below 1 TPM go last" },
];

export function categoryInfo(category) {
  return EXPRESSION_CATEGORIES[category] || EXPRESSION_CATEGORIES.missing;
}

// On-screen flag (requested by the AZ team): shown ONLY when the gene is off,
// i.e. below 1 TPM in DepMap and HPA. At 1 TPM or more only the TPM is shown.
// The full category (Below detection / Negligible / Low measurable / Moderate)
// is still used for ranking and is kept in the CSV export for analysis.
export const OFF_FLAG = { label: "< 1 TPM", color: "#b7791f", bg: "#fff4dc" };

export function isOff(info) {
  return !!info && (info.category === "below_detection" || info.category === "negligible");
}

export const DISAGREE_HINT = "DepMap and HPA disagree: one source is below 1 TPM, the other is not";

// "0.45 TPM" (mean of the sources) or "—" when there is no RNA value
export function formatTpm(info) {
  if (!info || info.tpm == null) return "—";
  return `${info.tpm < 10 ? info.tpm.toFixed(2) : info.tpm.toFixed(1)} TPM`;
}

// Tooltip text: the TPM of each source, plus why a flag or ⚠ is shown, e.g.
// "DepMap 0.6 · HPA 0.3 TPM\nBelow 1 TPM in DepMap and HPA"
export function expressionTooltip(info) {
  if (!info) return "";
  const names = { depmap: "DepMap", hpa: "HPA" };
  const sources = Object.entries(info.tpm_by_source || {})
    .map(([s, v]) => `${names[s] || s} ${v}`)
    .join(" · ");
  const lines = [sources ? `${sources} TPM` : "No RNA value"];
  if (isOff(info)) lines.push("Below 1 TPM in DepMap and HPA");
  if (info.sources_disagree) lines.push(DISAGREE_HINT);
  return lines.join("\n");
}
