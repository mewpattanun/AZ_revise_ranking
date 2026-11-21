// Compare the three ways of putting RNA on a common scale (sidebar button in the main page).
// The backend ranks the same query once per rna_method; this file lines the three
// rankings up cell line by cell line. No React in here, so it is easy to check.

export const RNA_METHODS = [
  {
    value: "robust_z",
    label: "Robust z (MAD)",
    short: "Robust z",
    valueLabel: "Robust z",
    formula: "(x - median) / (1.4826 x MAD)",
    color: "#5b2c83",
    bg: "#efe7f7",
  },
  {
    value: "standard_z",
    label: "Standard z (SD)",
    short: "Standard z",
    valueLabel: "z (SD)",
    formula: "(x - mean) / SD",
    color: "#185fa5",
    bg: "#e6f1fb",
  },
  {
    value: "log2_tpm",
    label: "log2(TPM + 1)",
    short: "log2 TPM",
    valueLabel: "log2(TPM+1)",
    formula: "log2(TPM + 1), no z-score",
    color: "#2d6a4f",
    bg: "#e8f4e8",
  },
];

// Ask the backend for this many cell lines per method (its maximum), so a cell line that is in
// the top K of one method still has a rank in the others.
export const LIST_LIMIT = 200;

const PAIRS = [
  ["robust_z", "standard_z"],
  ["robust_z", "log2_tpm"],
  ["standard_z", "log2_tpm"],
];

export function methodLabel(value) {
  return RNA_METHODS.find((m) => m.value === value)?.label || value;
}

// resultsByMethod: { robust_z: [row, ...], standard_z: [...], log2_tpm: [...] }, rows as returned by
//   POST /api/rank (rank 1 = best for the requested direction).
// k: how many top cell lines to compare. geneNames: the searched genes; the TPM is shown for a
//   single gene, and each method's RNA value (z or log2) is kept per gene.
export function buildComparison(resultsByMethod, k, geneNames = []) {
  const firstGene = geneNames.length === 1 ? geneNames[0] : null;
  const methods = RNA_METHODS.map((m) => m.value);
  const byMethod = {};
  const meta = new Map();
  for (const m of methods) {
    const rows = resultsByMethod[m] || [];
    byMethod[m] = new Map(rows.map((r) => [r.ach_id, r]));
    for (const r of rows) if (!meta.has(r.ach_id)) meta.set(r.ach_id, r);
  }

  const top = {};
  for (const m of methods) {
    top[m] = new Set((resultsByMethod[m] || []).filter((r) => r.rank <= k).map((r) => r.ach_id));
  }

  const ids = new Set();
  for (const m of methods) for (const id of top[m]) ids.add(id);

  const rows = [...ids].map((id) => {
    const first = meta.get(id);
    const ranks = {};
    const scores = {};
    const values = {};   // RNA consensus of that method: [{ gene, value }, ...]
    for (const m of methods) {
      const r = byMethod[m].get(id);
      ranks[m] = r ? r.rank : null;
      scores[m] = r ? r.score : null;
      values[m] = geneNames.map((g) => ({ gene: g, value: r?.expression?.[g]?.rna_value ?? null }));
    }
    const info = firstGene ? first.expression?.[firstGene] : null;
    return {
      id,
      name: first.cell_line_name || id,
      disease: first.primary_disease || null,
      tpm: info && info.tpm != null ? info.tpm : null,
      ranks,
      scores,
      values,
    };
  });

  const key = (v) => (v == null ? Infinity : v);
  rows.sort((a, b) => {
    for (const m of methods) {
      const d = key(a.ranks[m]) - key(b.ranks[m]);
      if (d !== 0 && Number.isFinite(d)) return d;
      if (key(a.ranks[m]) !== key(b.ranks[m])) return key(a.ranks[m]) === Infinity ? 1 : -1;
    }
    return a.id < b.id ? -1 : 1;
  });

  const overlaps = PAIRS.map(([a, b]) => ({
    a,
    b,
    shared: [...top[a]].filter((id) => top[b].has(id)).length,
    of: Math.max(top[a].size, top[b].size),
  }));

  // list shorter than the limit = every ranked cell line is in it, so a missing one is not ranked
  const complete = {};
  for (const m of methods) complete[m] = (resultsByMethod[m] || []).length < LIST_LIMIT;

  return { rows, overlaps, k, complete };
}

export function formatRnaValue(v) {
  return v == null ? "-" : v.toFixed(2);
}

// "7" inside the list, "> 200" for a cell line beyond it, "-" when the method does not rank it
export function rankText(rank, complete) {
  if (rank != null) return String(rank);
  return complete ? "-" : `> ${LIST_LIMIT}`;
}

function csvCell(v) {
  if (v == null) return "";
  const s = String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

export function comparisonCsv(comparison, genes) {
  const methods = RNA_METHODS.map((m) => m.value);
  const header = [
    "ach_id",
    "cell_line_name",
    "primary_disease",
    ...(genes.length === 1 ? ["tpm_mean"] : []),
    ...methods.flatMap((m) => [
      `rank_${m}`,
      ...genes.map((g) => (genes.length === 1 ? `rna_value_${m}` : `rna_value_${m}_${g.hugo}`)),
      `score_${m}`,
    ]),
  ];
  const lines = [header.join(",")];
  for (const r of comparison.rows) {
    const cells = [
      r.id,
      r.name,
      r.disease,
      ...(genes.length === 1 ? [r.tpm] : []),
      ...methods.flatMap((m) => [r.ranks[m], ...r.values[m].map((v) => v.value), r.scores[m]]),
    ];
    lines.push(cells.map(csvCell).join(","));
  }
  return lines.join("\n");
}
