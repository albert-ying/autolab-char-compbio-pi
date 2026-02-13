---
name: scientific-writing
description: Write publication-ready scientific manuscripts. IMRAD structure, clear figures, proper citations, LaTeX formatting. Optimized for computational biology journals.
metadata:
    skill-author: Albert Ying
---

# Scientific writing

## When to use

- Drafting or revising manuscript sections
- Structuring results around figures
- Writing methods with reproducibility in mind
- Responding to reviewer comments

## Paper structure (comp bio)

1. **Introduction**: broad context → gap (pivot with "However") → "Here, we..." → key findings
2. **Results**: one figure per subsection. State the question, describe the approach in one sentence, present the result, interpret it.
3. **Discussion**: summary of findings → comparison to prior work → limitations → future directions
4. **Methods**: enough detail to reproduce. Include software versions, parameters, random seeds.

## Writing rules

- One idea per sentence. Short declarative sentences.
- "We" for collective work.
- Every citation does argumentative work. If it doesn't support a specific claim, cut it.
- No significance inflation. State facts; let the reader judge importance.
- Title: short, no colon, conceptual rather than descriptive.

## Figure conventions

- Panel labels: bold uppercase (A, B, C)
- All axes labeled with units
- Color-blind safe palettes (viridis, colorbrewer)
- Resolution: 300 DPI minimum for raster, vector preferred
- Self-contained: a reader should understand the figure from the caption alone

## LaTeX tips

```latex
\begin{figure}[!ht]
\centering
\includegraphics[width=\textwidth]{figures/fig1.pdf}
\caption{\textbf{Title.} (A) Description. (B) Description.}
\label{fig:main}
\end{figure}
```
