import { useState, useCallback, useMemo } from "react";
import SearchBar from "./components/SearchBar/SearchBar.jsx";
import FilterPanel from "./components/FilterPanel/FilterPanel.jsx";
import ResultsTable from "./components/ResultsTable/ResultsTable.jsx";
import DetailPanel from "./components/DetailPanel/DetailPanel.jsx";
import CompareView from "./components/CompareView/CompareView.jsx";
import ChatWidget from "./components/ChatWidget/ChatWidget.jsx";
import DocsModal from "./components/DocsModal/DocsModal.jsx";
import MethodCompare from "./components/MethodCompare/MethodCompare.jsx";
import { RNA_METHODS, methodLabel as rnaMethodLabel } from "./compareMethods";
import { rankCellLines } from "./api";
import { exportResultsCsv } from "./exportCsv";
import { VERSION_LABEL, GIT_COMMIT, BUILD_DATE } from "./version";
import styles from "./App.module.css";

const DEFAULT_FILTERS = {
  w_rna: 0.7,
  w_protein: 0.3,
  mutation_mode: "ignore",
  driver_only: false,
  fusion_mode: "ignore",
  disease_filter: null,
  lineage_filter: null,
  subtype_filter: null,
  filter_mode: "filter",  // "filter" = hard exclude, "boost" = Q6 soft boost
  core_only: false,
  require_protein: false,
  top_n: 20,
  sources: ["depmap", "hpa"],
  assay_type: null,
  msi_max: null,
  cin_max: null,
  exclude_metabolite: null,
  metabolite_threshold: null,
  exclude_mirna: null,
  mirna_threshold: null,
  // What "low" means for LOW genes (AZ priority 5): "off" | "measurable"
  low_expression_mode: "off",
  // How each RNA source is put on a common scale: "robust_z" (default) | "standard_z" | "log2_tpm"
  rna_method: "robust_z",
};

export default function App() {
  /* Search state */
  const [genes, setGenes] = useState([]);
  const [filters, setFilters] = useState(DEFAULT_FILTERS);

  /* Results state */
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  // Snapshot of filters at the time of the last successful search — used to
  // detect when the current filters differ, so we can prompt the user to
  // refresh the results.
  const [lastSearchedFilters, setLastSearchedFilters] = useState(null);
  // Genes at the time of the last successful search — same purpose.
  const [lastSearchedGenes, setLastSearchedGenes] = useState(null);
  // Backend scoring version from the last /api/rank response (shown in the
  // version badge tooltip; older backends do not send it).
  const [scoringVersion, setScoringVersion] = useState(null);
  // Plain-language notes from the backend about the last result, e.g. a LOW
  // gene that does not separate the cell lines shown.
  const [notes, setNotes] = useState([]);

  /* Selection & panels */
  const [selected, setSelected] = useState([]);
  const [detailId, setDetailId] = useState(null);
  const [showCompare, setShowCompare] = useState(false);
  const [methodCompareQuery, setMethodCompareQuery] = useState(null); // query shown in the RNA method comparison
  const [docsTab, setDocsTab] = useState(null); // null = closed, "docs" | "sources" | "methods" | "about"

  /* Gene management */
  const addGene = useCallback((gene) => {
    setGenes((prev) => [...prev, gene]);
  }, []);

  const removeGene = useCallback((hugo) => {
    setGenes((prev) => prev.filter((g) => g.hugo !== hugo));
  }, []);

  const handleFilterChange = useCallback((newFilters) => {
    setFilters(newFilters);
  }, []);

  /* Request body for /api/rank from the current genes and filters */
  function buildApiParams() {
    // Protein is always included in scoring (only the 3 RNA sources are
    // user-toggleable), so ensure "protein" is always in sources.
    const sourcesWithProtein = filters.sources.includes("protein")
      ? filters.sources
      : [...filters.sources, "protein"];

    // Route disease/lineage/subtype to filter or boost params based on mode
    const isBoost = filters.filter_mode === "boost";
    const apiParams = {
      genes: genes.map((g) => ({ hugo: g.hugo, direction: g.direction })),
      ...filters,
      sources: sourcesWithProtein,
      // In boost mode: clear hard filters, set Q6 targets
      disease_filter: isBoost ? null : filters.disease_filter,
      lineage_filter: isBoost ? null : filters.lineage_filter,
      subtype_filter: isBoost ? null : filters.subtype_filter,
      target_disease: isBoost ? filters.disease_filter : null,
      target_lineage: isBoost ? filters.lineage_filter : null,
      target_subtype: isBoost ? filters.subtype_filter : null,
      q6_boost_enabled: isBoost,
    };
    // frontend-only fields, don't send them
    delete apiParams.filter_mode;
    delete apiParams.driver_only;
    return apiParams;
  }

  /* Rank the same query with the three RNA methods and show them side by side */
  function openMethodCompare() {
    if (genes.length === 0) return;
    setMethodCompareQuery({ genes, params: buildApiParams() });
  }

  /* Run ranking */
  async function handleSearch() {
    if (genes.length === 0) return;
    setLoading(true);
    setError(null);
    setResults(null);
    setNotes([]);
    setSelected([]);

    try {
      const res = await rankCellLines(buildApiParams());
      setResults(res.results);
      setScoringVersion(res.query?.scoring_version ?? null);
      setNotes(res.notes ?? []);
      setLastSearchedFilters(filters);
      setLastSearchedGenes(genes);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  }

  /* Selection handlers */
  function handleToggleSelect(id) {
    if (id === "COMPARE") {
      setShowCompare(true);
      return;
    }
    if (id === "CLEAR") {
      setSelected([]);
      return;
    }
    setSelected((prev) =>
      prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]
    );
  }

  function handleSelectAll(selectAll) {
    if (!displayedResults) return;
    setSelected(selectAll ? displayedResults.map((r) => r.ach_id) : []);
  }

  // Apply the "Driver only" filter on the client after ranking. Only active
  // when Mutation Filter = Include and Driver only is checked. Ranks are
  // recomputed so the displayed list is contiguous (1, 2, 3, ...).
  const displayedResults = useMemo(() => {
    if (!results) return null;
    const shouldFilter =
      filters.mutation_mode === "include" && filters.driver_only;
    if (!shouldFilter) return results;
    const filtered = results.filter(
      (r) => Array.isArray(r.mutations) && r.mutations.some((m) => m.is_driver)
    );
    return filtered.map((r, i) => ({ ...r, rank: i + 1 }));
  }, [results, filters.mutation_mode, filters.driver_only]);

  /* Highest score in the current results — used so the detail panel can colour
     its score with the same ratio-based thresholds as the results table. */
  const maxScore =
    displayedResults && displayedResults.length > 0
      ? Math.max(...displayedResults.map((r) => r.score))
      : 0;

  /* Detect pending changes: filters (or gene list) differ from the last search.
     driver_only is excluded because it filters live on the frontend and doesn't
     need a re-search. */
  const hasPendingChanges = useMemo(() => {
    if (!results || !lastSearchedFilters) return false;
    const stripDriverOnly = (f) => {
      const { driver_only, ...rest } = f;
      return rest;
    };
    const filtersChanged =
      JSON.stringify(stripDriverOnly(filters)) !==
      JSON.stringify(stripDriverOnly(lastSearchedFilters));
    const genesChanged =
      JSON.stringify(genes) !== JSON.stringify(lastSearchedGenes);
    return filtersChanged || genesChanged;
  }, [filters, lastSearchedFilters, genes, lastSearchedGenes, results]);

  return (
    <div className={styles.app}>
      {/* ── Header ── */}
      <header className={styles.header}>
        <div className={styles.headerLeft}>
          <div className={styles.logo}>
            Cell<span className={styles.logoAccent}>Line</span>Finder
          </div>
          <span
            className={styles.headerTag}
            title={
              `Website ${VERSION_LABEL}\nCommit ${GIT_COMMIT}\nBuilt ${BUILD_DATE}` +
              (scoringVersion ? `\nScoring ${scoringVersion}` : "\nScoring version: run a search to see it")
            }
          >
            {/* Commit hash and build date are in the tooltip (hover) */}
            {VERSION_LABEL} · Bristol × AstraZeneca
          </span>
        </div>
        <div className={styles.headerRight}>
          <button className={styles.headerLink} onClick={() => setDocsTab("docs")}>Documentation</button>
          <button className={styles.headerLink} onClick={() => setDocsTab("sources")}>Data Sources</button>
          <button className={styles.headerLink} onClick={() => setDocsTab("methods")}>Methods</button>
          <button className={styles.headerLink} onClick={() => setDocsTab("about")}>About</button>
        </div>
      </header>

      <div className={styles.main}>
        {/* ── Left sidebar: Search Query ── */}
        <aside className={styles.sidebar}>
          <div className={styles.sidebarSection}>
            <h2 className={styles.sidebarTitle}>Search Query</h2>
            <SearchBar genes={genes} onAddGene={addGene} onRemoveGene={removeGene} />
          </div>

          <div className={styles.sidebarSection}>
            <h2 className={styles.sidebarTitle}>RNA Normalisation</h2>
            <div className={styles.methodToggle}>
              {RNA_METHODS.map((m) => (
                <button
                  key={m.value}
                  className={`${styles.methodBtn} ${filters.rna_method === m.value ? styles.methodActive : ""}`}
                  onClick={() => handleFilterChange({ ...filters, rna_method: m.value })}
                >
                  {m.short}
                </button>
              ))}
            </div>
            <p className={styles.methodHint}>
              How each RNA source (DepMap, HPA) is put on a common scale before the two are averaged.
              <b> Robust z</b> = (x - median) / (1.4826 x MAD), the default.
              <b> Standard z</b> = (x - mean) / SD.
              <b> log2 TPM</b> = log2(TPM + 1), no z-score.
              Both z-scores take their reference from cell lines at or above 1 TPM only.
            </p>
            <button
              className={styles.compareBtn}
              onClick={openMethodCompare}
              disabled={genes.length === 0}
            >
              Compare the 3 methods
            </button>
          </div>

          <div className={styles.sidebarSection}>
            <FilterPanel
              filters={filters}
              onChange={handleFilterChange}
              hasLowGene={genes.some((g) => g.direction === "low")}
            />
          </div>

          <div className={styles.sidebarSection}>
            <button
              className={styles.searchBtn}
              onClick={handleSearch}
              disabled={genes.length === 0 || loading}
            >
              {loading ? "Searching..." : "Find Cell Lines"}
            </button>
            {error && <div className={styles.error}>{error}</div>}
          </div>
        </aside>

        {/* ── Right content: Results ── */}
        <main className={styles.content}>
          {displayedResults && displayedResults.length > 0 && (
            <div className={styles.resultsHeader}>
              <h2 className={styles.resultsTitle}>Results</h2>
              <div style={{ display: "flex", alignItems: "center", gap: "12px", marginLeft: "auto" }}>
                {hasPendingChanges && (
                  <span style={{ fontSize: "13px", color: "#f0b429", fontWeight: 500 }}>
                    Filters changed — refresh to update
                  </span>
                )}
                <button
                  onClick={handleSearch}
                  disabled={loading || genes.length === 0}
                  style={{
                    padding: "6px 14px",
                    fontSize: "13px",
                    borderRadius: "6px",
                    border: hasPendingChanges ? "1px solid #f0b429" : "1px solid #ccc",
                    background: hasPendingChanges ? "#fff8e6" : "#fff",
                    color: hasPendingChanges ? "#b8860b" : "#555",
                    cursor: loading ? "wait" : "pointer",
                    fontWeight: hasPendingChanges ? 600 : 400,
                  }}
                >
                  {loading ? "Refreshing..." : "↻ Refresh"}
                </button>
                <button
                  onClick={() => exportResultsCsv({ results: displayedResults, genes, filters })}
                  style={{
                    padding: "6px 14px",
                    fontSize: "13px",
                    borderRadius: "6px",
                    border: "1px solid #ccc",
                    background: "#fff",
                    color: "#555",
                    cursor: "pointer",
                  }}
                  title="Download results as CSV"
                >
                  Export CSV
                </button>
                <span className={styles.resultsMeta}>
                  Showing <b>{displayedResults.length}</b> cell lines · Ranked by{" "}
                  <b>
                    {/* Backend defaults: scoring_method=robust_z, multi_gene_method=geometric */}
                    {rnaMethodLabel(lastSearchedFilters?.rna_method || "robust_z")}
                    {lastSearchedGenes?.length > 1 && " + Geometric Mean"}
                  </b>
                </span>
              </div>
            </div>
          )}

          {/* Empty-state: a search was run but no cell lines matched. */}
          {results && (!displayedResults || displayedResults.length === 0) && !loading && (
            <div className={styles.emptyState}>
              <div className={styles.emptyIcon}>🔍</div>
              <div className={styles.emptyTitle}>No cell lines match your search criteria</div>
              {hasPendingChanges && (
                <button
                  onClick={handleSearch}
                  disabled={loading || genes.length === 0}
                  style={{
                    marginTop: "12px",
                    padding: "8px 16px",
                    fontSize: "13px",
                    borderRadius: "6px",
                    border: "1px solid #f0b429",
                    background: "#fff8e6",
                    color: "#b8860b",
                    fontWeight: 600,
                    cursor: "pointer",
                  }}
                >
                  ↻ Refresh with current filters
                </button>
              )}
            </div>
          )}

          {!results && !loading && (
            <div className={styles.emptyState}>
              <div className={styles.emptyIcon}>🔬</div>
              <div className={styles.emptyTitle}>Add genes to begin</div>
              <p className={styles.emptyText}>
                Search for target genes on the left, set your filters, then click
                "Find Cell Lines" to rank cell lines by multi-omics expression.
              </p>
            </div>
          )}

          {/* Notes about the result (shown only when there are results) */}
          {!loading && displayedResults?.length > 0 && notes.length > 0 && (
            <div className={styles.resultNotes} role="note">
              {notes.map((n) => (
                <p key={n} className={styles.resultNote}>
                  <span className={styles.resultNoteIcon} aria-hidden="true">ℹ</span>
                  {n}
                </p>
              ))}
            </div>
          )}

          <ResultsTable
            results={displayedResults}
            loading={loading}
            selected={selected}
            onToggleSelect={handleToggleSelect}
            onSelectAll={handleSelectAll}
            onRowClick={setDetailId}
            mutationMode={filters.mutation_mode}
            fusionMode={filters.fusion_mode}
            assayType={filters.assay_type}
            boostMode={filters.filter_mode === "boost"}
          />
        </main>
      </div>

      {/* Detail slide-in */}
      {detailId && (
        <DetailPanel
          achId={detailId}
          genes={genes}
          resultRow={displayedResults?.find((r) => r.ach_id === detailId)}
          maxScore={maxScore}
          rnaMethod={lastSearchedFilters?.rna_method || "robust_z"}
          wRna={filters.w_rna}
          wProtein={filters.w_protein}
          mutationMode={filters.mutation_mode}
          fusionMode={filters.fusion_mode}
          onClose={() => setDetailId(null)}
        />
      )}

      {/* Compare modal */}
      {showCompare && (
        <CompareView
          achIds={selected}
          genes={genes}
          onClose={() => setShowCompare(false)}
        />
      )}

      {/* RNA method comparison modal */}
      {methodCompareQuery && (
        <MethodCompare
          genes={methodCompareQuery.genes}
          params={methodCompareQuery.params}
          onClose={() => setMethodCompareQuery(null)}
        />
      )}

      {/* Docs modal */}
      {docsTab && (
        <DocsModal initialTab={docsTab} onClose={() => setDocsTab(null)} />
      )}

      {/* Chat widget */}
      <ChatWidget />
    </div>
  );
}
