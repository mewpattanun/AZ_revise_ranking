// Expression categories from raw TPM (AZ priority 5, set in backend config.py).
// A cell line gets the highest band any source puts it in, so it only counts
// as "off" when DepMap and HPA both agree.

export const EXPRESSION_CATEGORIES = {
  below_detection: { label: "Below detection", short: "Off", color: "#c0392b", bg: "#fde8e8", rule: "< 0.1 TPM in every source" },
  negligible: { label: "Negligible", short: "Off", color: "#b7791f", bg: "#fff4dc", rule: "< 1 TPM in every source" },
  low_relative: { label: "Low for this gene", short: "Off", color: "#9c4221", bg: "#feebdc", rule: "≥ 1 TPM but off under the chosen off cutoff" },
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

// Below 1 TPM in DepMap and HPA (always counts as off)
export function isOff(info) {
  return !!info && (info.category === "below_detection" || info.category === "negligible");
}

// 1 TPM or more, but off for THIS gene under a gene-relative off cutoff
// (only when the user picks a method other than "Fixed 1 TPM").
export const RELATIVE_OFF_FLAG = { color: "#9c4221", bg: "#feebdc" };
export function isRelativeOff(info) {
  return !!info && info.category === "low_relative";
}
export function relativeOffLabel(hugo) {
  return `Low for ${hugo}`;
}

// How a gene is judged "off" (backend config.OFF_METHODS). Below 1 TPM is
// always off; the other methods can also call a line off at >= 1 TPM when it
// is low for that gene (e.g. VIM 2.5 TPM in MCF7 vs ~1,000 in mesenchymal lines).
export const OFF_METHODS = [
  { key: "fixed", label: "Fixed 1 TPM (default)",
    hint: "Off = below 1 TPM in DepMap and HPA." },
  { key: "strict_mixture", label: "Strict mixture (recommended alternative)",
    hint: "Also off: lines in the low group of a gene with a clear on/off split (e.g. VIM). Genes without a clear split use 1 TPM." },
  { key: "distance_from_on", label: "Distance from on-group",
    hint: "Also off: lines more than 3 SD below the group of lines that express the gene." },
  { key: "mixture", label: "Mixture model (experimental)",
    hint: "Also off: lines in the low group of a 2-group split. Can call normal lines off (e.g. ERBB2)." },
  { key: "percent_of_reference", label: "% of reference (experimental)",
    hint: "Also off: lines below 10% of the gene's reference level (90th percentile of expressing lines). Calls many lines off." },
  { key: "custom", label: "Custom cutoff",
    hint: "Off = below the TPM you enter (at least 1), for every gene in the search." },
];

export function offMethodLabel(key) {
  return (OFF_METHODS.find((m) => m.key === key) || OFF_METHODS[0]).label.replace(/ \(.*\)$/, "");
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
  if (isRelativeOff(info) && info.off_reason) lines.push(`Off for this gene: ${info.off_reason}`);
  if (info.sources_disagree) lines.push(DISAGREE_HINT);
  return lines.join("\n");
}
