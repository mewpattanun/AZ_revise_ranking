import { useState } from "react";
import styles from "./DocsModal.module.css";
import { VERSION_LABEL } from "../../version";

const TABS = [
  { key: "docs", label: "Documentation" },
  { key: "sources", label: "Data Sources" },
  { key: "methods", label: "Methods" },
  { key: "about", label: "About" },
];

export default function DocsModal({ initialTab = "docs", onClose }) {
  const [tab, setTab] = useState(initialTab);

  return (
    <div className={styles.overlay} onClick={onClose}>
      <div className={styles.modal} onClick={(e) => e.stopPropagation()}>
        {/* Header */}
        <div className={styles.header}>
          <h2 className={styles.title}>CellLineFinder Documentation</h2>
          <button className={styles.closeBtn} onClick={onClose}>×</button>
        </div>

        {/* Tabs */}
        <div className={styles.tabs}>
          {TABS.map((t) => (
            <button
              key={t.key}
              className={`${styles.tab} ${tab === t.key ? styles.tabActive : ""}`}
              onClick={() => setTab(t.key)}
            >
              {t.label}
            </button>
          ))}
        </div>

        {/* Content */}
        <div className={styles.body}>
          {tab === "docs" && <DocsContent />}
          {tab === "sources" && <SourcesContent />}
          {tab === "methods" && <MethodsContent />}
          {tab === "about" && <AboutContent />}
        </div>
      </div>
    </div>
  );
}

/* ── Tab content components ─────────────────────────────── */

function DocsContent() {
  return (
    <div className={styles.content}>
      <h3>Getting Started</h3>
      <p>
        CellLineFinder is a multi-omics cell line recommendation platform designed
        to help researchers identify the most suitable cancer cell lines for their
        experiments based on gene expression, protein abundance, mutation status,
        and fusion data.
      </p>

      <h4>How to Search</h4>
      <ol>
        <li>
          <strong>Add target genes</strong> — Type a HUGO gene symbol (e.g. EGFR, TP53)
          into the search bar and click Add. Set each gene to HIGH (prefer high expression)
          or LOW (prefer low expression).
        </li>
        <li>
          <strong>Configure filters</strong> — Optionally filter by disease, lineage,
          subtype, mutation/fusion status, or adjust RNA/Protein weights.
        </li>
        <li>
          <strong>Find Cell Lines</strong> — Click the button to run the ranking pipeline.
          Results are sorted by score (best match first).
        </li>
        <li>
          <strong>Explore details</strong> — Click any row to see the full multi-omics
          evidence breakdown, including per-source scores, protein detection, and mutation details.
        </li>
        <li>
          <strong>Compare</strong> — Select 2 or more cell lines using the checkboxes,
          then click Compare to see a side-by-side view.
        </li>
      </ol>

      <h4>Multi-Gene Search</h4>
      <p>
        When searching for multiple genes, every cell line gets a 0 to 1 score for
        each gene (for a LOW gene the score is flipped, so higher is always better).
        The gene scores are combined with a weighted geometric mean, so one gene that
        fails pulls the overall score down. A cell line with no data for a gene is
        scored on the other genes and marked with ⚠. Per-gene ranks are shown in each
        results row so you can see which genes a cell line is strong or weak in.
      </p>

      <h4>Score Interpretation</h4>
      <table className={styles.infoTable}>
        <thead>
          <tr><th>Score Range</th><th>Interpretation</th></tr>
        </thead>
        <tbody>
          <tr><td className={styles.scoreHigh}>0.7 – 1.0</td><td>Strong match — high expression/protein for target genes</td></tr>
          <tr><td className={styles.scoreMid}>0.4 – 0.7</td><td>Moderate match — consider reviewing the evidence breakdown</td></tr>
          <tr><td className={styles.scoreLow}>0.0 – 0.4</td><td>Weak match — likely low expression for one or more genes</td></tr>
        </tbody>
      </table>

      <h4>Chatbot Assistant</h4>
      <p>
        Use the chat widget (bottom-right corner) to ask biology questions about
        your results, get explanations of cell line characteristics, or request
        help interpreting the ranking output.
      </p>
    </div>
  );
}

function SourcesContent() {
  return (
    <div className={styles.content}>
      <h3>Data Sources</h3>
      <p>
        CellLineFinder integrates data from multiple public databases, harmonised
        into a unified DuckDB warehouse. All data is mapped to a common cell line
        identifier (ACH ID from DepMap).
      </p>

      <h4>RNA Expression</h4>
      <table className={styles.infoTable}>
        <thead>
          <tr><th>Source</th><th>Description</th><th>Coverage</th></tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>DepMap</strong></td>
            <td>Cancer Cell Line Encyclopedia (CCLE) RNA-seq TPM values from the Broad Institute DepMap portal</td>
            <td>~1,400 cell lines</td>
          </tr>
          <tr>
            <td><strong>HPA</strong></td>
            <td>Human Protein Atlas RNA-seq expression data (TPM) across cancer cell lines</td>
            <td>~70 cell lines</td>
          </tr>
          <tr>
            <td><strong>GEO</strong></td>
            <td>Gene Expression Omnibus curated RNA-seq datasets, reprocessed to TPM. Not used in ranking (heterogeneous studies); shown as GEO support instead</td>
            <td>~250 cell lines</td>
          </tr>
        </tbody>
      </table>

      <h4>Protein Abundance</h4>
      <table className={styles.infoTable}>
        <thead>
          <tr><th>Source</th><th>Description</th><th>Coverage</th></tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>Gygi Proteomics</strong></td>
            <td>Mass spectrometry–based protein intensity measurements from the Gygi lab (via DepMap)</td>
            <td>~375 cell lines</td>
          </tr>
        </tbody>
      </table>

      <h4>Genomic Variants</h4>
      <table className={styles.infoTable}>
        <thead>
          <tr><th>Source</th><th>Description</th></tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>Mutations</strong></td>
            <td>Somatic mutations from DepMap, including driver/hotspot/damaging annotations and variant allele frequency</td>
          </tr>
          <tr>
            <td><strong>Fusions</strong></td>
            <td>Gene fusions from DepMap, including confidence level, reading frame, and fusion fragments per million (FFPM)</td>
          </tr>
        </tbody>
      </table>

      <h4>Cell Line Metadata</h4>
      <p>
        Cell line annotations (disease, lineage, subtype, sex, growth pattern)
        are sourced from the DepMap Model table and used for filtering and the
        tissue-match (Q6) boost.
      </p>
    </div>
  );
}

function MethodsContent() {
  return (
    <div className={styles.content}>
      <h3>Ranking Methodology</h3>

      <h4>Step 1: Per-Gene Expression Scoring</h4>
      <p>
        DepMap and HPA derive from the same underlying data, so they are not
        treated as independent evidence. Each source is standardised separately
        (a gene-wise robust z-score of TPM across cell lines, using the
        median and MAD) and the matched robust z-scores are averaged into one RNA
        consensus per cell line. The difference between the two sources is kept
        as a discordance QC value (shown in the detail panel); it does not change
        the score.
      </p>
      <p>
        The expression score comes from the <strong>size</strong> of the consensus robust z,
        not from its rank: z is compressed with asinh (so one extreme cell line does not
        squash all the others) and scaled to 0–1. A cell line far above the rest
        therefore scores clearly higher than one just above the median, and two cell
        lines with almost the same expression get almost the same score. The rank
        shown in the table is simply the sorted score. A cell line present in only one
        source is scored on that source alone, so missing one source is not penalised.
      </p>
      <p>
        The median and MAD of a gene are estimated after removing the cell lines below 1 TPM
        (the same "gene is off" cutoff as the expression categories), so tiny values such as
        0.01 in DepMap and the 0.1 steps in HPA do not set the scale: both sources give
        similar z-scores for the same cell line. The z-score itself still uses the real
        value of every cell line, so the order of cell lines is unchanged and the off lines
        get a z below the expressed ones. If MAD is 0 for a gene in a source (most cell
        lines at or above 1 TPM share the same value), the standard deviation is used
        instead and the cell line is flagged
        with ⚠ in the results table and the detail panel. The legacy rank-based methods
        (Z-Score / Percentile → RRF) can still be selected for comparison.
      </p>

      <h4>Step 1b: Expression Categories and Low Expression</h4>
      <p>
        Every cell line also gets an expression category from its raw TPM:{" "}
        <strong>Below detection</strong> (&lt; 0.1 TPM in every source),{" "}
        <strong>Negligible</strong> (&lt; 1 TPM in every source),{" "}
        <strong>Low measurable</strong> (1–3 TPM) or <strong>Moderate</strong> (&gt; 3 TPM).
        A cell line is only called off when DepMap and HPA both agree; if one source says
        off and the other says expressed, it is marked with ⚠.
      </p>
      <p>
        The Expression column shows the TPM (mean of DepMap and HPA; hover for each
        source). A <strong>&lt; 1 TPM</strong> flag is shown only when the gene is off,
        i.e. below 1 TPM in both sources; at 1 TPM or more only the number is shown.
        The full category is still used for ranking and is included in the CSV export.
      </p>
      <p>
        <strong>Off cutoff.</strong> A fixed 1 TPM cutoff cannot see that a gene can be
        off at more than 1 TPM: epithelial MCF7 has VIM 2.5 TPM, while mesenchymal lines
        have about 1,000. When a LOW gene is searched, the <strong>Off cutoff</strong>{" "}
        option chooses how "off" is decided. Below 1 TPM is always off; the other methods
        can also call a line off at 1 TPM or more, shown as <strong>Low for GENE</strong>{" "}
        (hover or open the detail panel for the reason). With "TPM below 1", lines below
        1 TPM come first (tied), then the Low for GENE lines ordered by their score
        (lowest RNA first; protein counts where it was measured):
      </p>
      <ul>
        <li><strong>Fixed 1 TPM</strong> (default): below 1 TPM only.</li>
        <li><strong>Strict mixture</strong>: for genes whose cell lines split clearly into
          an off and an on group (two Gaussian groups of log2(TPM + 1), separation
          Ashman's D &gt; 2), lines that are at least 90% likely to be in the low group.
          Other genes keep 1 TPM.</li>
        <li><strong>Distance from on-group</strong>: more than 3 SD below the group of
          lines that express the gene.</li>
        <li><strong>Mixture model</strong> (experimental): the low group of any clear
          enough 2-group split; can call normal lines off (e.g. ERBB2).</li>
        <li><strong>% of reference</strong> (experimental): below 10% of the gene's
          reference level (90th percentile of lines at 1 TPM or more); calls many lines off.</li>
        <li><strong>Custom cutoff</strong>: below a TPM you choose (at least 1), for example
          10 TPM for ESR1.</li>
      </ul>
      <p>
        On 29 cell lines with known biology, Fixed 1 TPM caught 8 of 13 known negatives,
        Strict mixture 10 of 13 with no false calls, and % of reference 12 of 13 but with 2
        false calls (see the threshold benchmark in the repository).
      </p>
      <p>
        For a LOW gene the category is applied first, then the score: the lowest number is
        not automatically the best low-expression candidate, because values below 1 TPM are
        mostly noise. With <strong>TPM below 1</strong> (default) the off cell lines come
        first and tie on RNA, so protein (where measured) orders them. With{" "}
        <strong>TPM greater than 1</strong> the 1–3 TPM cell lines come first (those where
        DepMap and HPA agree before those where one source says off) and off cell lines go
        last. The 1 TPM cutoff was checked against known negative cell lines and
        DepMap–HPA agreement (see the threshold benchmark in the repository).
      </p>
      <p>
        The same order applies in multi-gene searches: with <strong>TPM below 1</strong>,
        every cell line where all LOW genes are off comes before any line that expresses
        one of them, and the geometric mean (Step 4) only orders lines within each group.
        For example, in ERBB2 HIGH + ESR1 LOW, HCC1954 (very high ERBB2 but ESR1 3 TPM)
        comes after every ESR1-off line.
      </p>

      <h4>Step 2: Protein Scoring</h4>
      <p>
        Protein is scored the same way as RNA: a gene-wise robust z-score of the
        protein intensities across cell lines (median and MAD), compressed with asinh
        and scaled to 0–1, so a cell line with much higher protein than the rest scores
        clearly higher. If MAD is 0 the standard deviation is used and the cell line is
        flagged with ⚠. (The legacy methods still score protein with Z-Score +
        Percentile → RRF.)
      </p>

      <h4>Step 3: RNA + Protein Combination</h4>
      <p>
        RNA and protein scores are min-max normalised to [0, 1], then combined
        with user-configurable weights (default: 70% RNA, 30% Protein):
      </p>
      <code className={styles.formula}>
        combined_score = (w_rna × RNA_norm + w_protein × Protein_norm) / (w_rna + w_protein)
      </code>
      <p>
        Cell lines with only RNA or only protein data are scored on whichever
        layer is available, so missing data is not counted as zero. A layer with
        weight 0 is not used: with protein at 100%, cell lines without protein data
        are left out. The Coverage column shows how much of the weight a cell line
        has data for (1.00 = RNA and protein).
      </p>

      <h4>Step 4: Multi-Gene Combination</h4>
      <p>
        For multi-gene searches, the combined score of each gene (Step 3, between 0 and 1,
        flipped for LOW genes) is called D. The genes are combined with a
        <strong> weighted geometric mean</strong>:
      </p>
      <code className={styles.formula}>
        score = exp( Σ w_g × log(D_g + ε) / Σ w_g )
      </code>
      <p>
        w_g is the weight of each gene (equal by default) and ε is a tiny number (0.000001)
        that only avoids log(0). A single poor gene lowers the score a lot. Example with
        three genes: a cell line with 1.00, 0.98 and 0.05 has an arithmetic mean of 0.677 but
        a geometric mean of 0.366, while a balanced cell line with 0.82, 0.80 and 0.78 scores
        0.800 either way, so the balanced one moves above.
      </p>
      <p>
        If a cell line has no score for a gene, that gene is skipped for it (the other
        weights are rescaled) and the cell line is marked with ⚠. The older summed RRF rule is
        still available for comparison (multi_gene_method = summed_rrf).
      </p>

      <h4>Step 5: Hard Filters</h4>
      <p>
        Mutation and fusion filters (include/exclude/ignore) are applied before
        multi-gene fusion. Instability (MSI/CIN), metabolite, and miRNA filters
        are also available as hard exclusion criteria.
      </p>

      <h4>Step 6: Tissue-Match Soft Boost (Q6)</h4>
      <p>
        When a target disease, lineage, or subtype is specified, a Q6 score (0–1)
        rewards cell lines matching the target tissue. The boost formula multiplies
        the base score by (0.7 + 0.3 × q6_score), so matching cell lines receive up
        to a 30% score increase without completely overriding expression-based ranking.
      </p>

      <h4>Step 7: Assay Compatibility Warnings (Q7)</h4>
      <p>
        When an assay type is selected, cell lines are checked for compatibility
        based on growth pattern (adherent vs. suspension). Warnings are displayed
        but do not affect ranking order.
      </p>

      <h4>Step 8: GEO Support</h4>
      <p>
        GEO is not used in the score or the rank. It is only shown in the detail panel
        as extra evidence. GEO samples come from different studies, so each study is
        checked on its own: the cell line is compared with the other cell lines in the
        same study (robust z-score). When you search for HIGH, a study supports the cell
        line if its z is +1 or more, and contradicts it if z is -1 or less. For LOW it
        is the other way round.
      </p>
      <ul>
        <li><strong>Supportive</strong>: more studies support than contradict, and at least half support (e.g. 3/4)</li>
        <li><strong>Contradictory</strong>: more studies contradict than support, and at least half contradict</li>
        <li><strong>Mixed</strong>: anything else</li>
        <li><strong>Unavailable</strong>: no usable GEO study for this cell line</li>
      </ul>
      <p>
        x/y means x studies support out of y studies that include the cell line. The
        z-score of each study is listed in the detail panel.
      </p>

      <h4>Similarity Scoring (Similar Alternatives)</h4>
      <p>
        The detail panel offers a "Similar Alternatives" section that finds cell lines
        with the most similar molecular profile to the selected one. This helps researchers
        identify backup options or related models.
      </p>
      <p>
        The similarity pipeline works across multiple omics layers:
      </p>
      <ul>
        <li>
          <strong>Expression similarity</strong> — Built from the top 1,000 most variable
          genes (by variance across cell lines) in DepMap RNA-seq. Each cell line becomes a
          1,000-dimensional vector of z-scored log₂(TPM+1) values. Cosine similarity measures
          how closely two cell lines' expression profiles match.
        </li>
        <li>
          <strong>Protein similarity</strong> — Z-scored protein intensity vectors from Gygi
          proteomics, compared via cosine similarity. Combined with expression at a 70:30 weighting.
        </li>
        <li>
          <strong>Metabolomics &amp; miRNA</strong> — Reported as separate similarity figures
          alongside each result but not blended into the headline score, because testing showed
          they separate lineages far more weakly than expression.
        </li>
      </ul>
      <p>
        Cell lines sharing a lineage of derivation (parent–derivative relationships) are
        flagged with a "Derivative" badge so researchers can distinguish genuinely independent
        alternatives from related models.
      </p>

      <h4>Reciprocal Rank Fusion (RRF)</h4>
      <p>
        RRF is a well-established rank aggregation method from information retrieval
        (Cormack, Clarke & Butt, 2009). For each item, the RRF score is:
      </p>
      <code className={styles.formula}>
        score = Σ 1 / (k + rank_i)
      </code>
      <p>
        The parameter k=60 controls how much weight is given to top-ranked items.
        This non-linear formula ensures that being near the top of a ranking contributes
        far more than being in the middle or bottom.
      </p>
    </div>
  );
}

function AboutContent() {
  return (
    <div className={styles.content}>
      <h3>About CellLineFinder</h3>
      <p>
        CellLineFinder {VERSION_LABEL} is a multi-omics cell line recommendation platform
        developed as a collaborative project between the University of Bristol
        and AstraZeneca. The platform helps cancer researchers identify the most
        appropriate cell lines for their experimental needs by integrating and
        ranking data across multiple omics layers.
      </p>

      <h4>Team — Group 34</h4>
      <table className={styles.infoTable}>
        <thead>
          <tr><th>Member</th><th>Contributions</th></tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>Prab Wongsekleo</strong></td>
            <td>Full-stack development, deployment (EC2/Nginx/CI), LLM chatbot integration</td>
          </tr>
          <tr>
            <td><strong>Nattawas</strong></td>
            <td>LLM chatbot, biology questions drafting, data harmonisation</td>
          </tr>
          <tr>
            <td><strong>Victor</strong></td>
            <td>Data cleaning, similarity scoring, supervisor communication</td>
          </tr>
          <tr>
            <td><strong>Mew</strong></td>
            <td>Base ranking algorithm, project leadership, project management</td>
          </tr>
        </tbody>
      </table>

      <h4>Supervisor</h4>
      <p>Daniel D'Andrea</p>

      <h4>Technology Stack</h4>
      <ul>
        <li><strong>Backend:</strong> Python, FastAPI, DuckDB, Parquet</li>
        <li><strong>Frontend:</strong> React, Vite</li>
        <li><strong>Deployment:</strong> AWS EC2, Nginx, GitHub Actions CI/CD</li>
        <li><strong>AI:</strong> OpenAI GPT API for chatbot</li>
      </ul>

      <h4>Academic Context</h4>
      <p>
        MSc Data Science final project, University of Bristol, 2024–2025.
        Developed in collaboration with AstraZeneca to address the challenge
        of selecting appropriate cancer cell lines from thousands of candidates
        across multiple data modalities.
      </p>
    </div>
  );
}
