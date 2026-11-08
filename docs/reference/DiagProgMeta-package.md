# Meta-Analysis of Diagnostic and Prognostic Studies with Sensitivity Analysis for Publication Bias on the SROC and Time-Dependent SROC Curves

This package implements the meta-analysis of diagnostic test accuracy
(DTA) studies and prognostic studies by the summary receiver operating
characteristic (SROC) curve and the area under the SROC curve (SAUC),
together with the likelihood-based sensitivity analysis for publication
bias.

- **Diagnostic studies**: the bivariate random-effects model of the
  logit-transformed sensitivity and specificity (Reitsma et al., 2005)
  [doi:10.1016/j.jclinepi.2005.02.022](https://doi.org/10.1016/j.jclinepi.2005.02.022)
  , and the sensitivity analysis for publication bias on the SROC curve
  and SAUC (Zhou et al., 2023)
  [doi:10.1002/sim.9643](https://doi.org/10.1002/sim.9643) .

- **Prognostic studies**: the time-dependent SROC curve and SAUC, where
  the time-dependent sensitivity and specificity are reconstructed from
  the published survival probabilities (Hattori and Zhou, 2016)
  [doi:10.1002/sim.7029](https://doi.org/10.1002/sim.7029) and the
  sensitivity analysis for publication bias on the time-dependent SROC
  curve and SAUC (Zhou et al., 2025)
  [doi:10.1017/rsm.2025.14](https://doi.org/10.1017/rsm.2025.14) .

This package implements the likelihood-based sensitivity analysis for
publication bias on the summary receiver operating characteristic (SROC)
curve and the area under the SROC curve (SAUC) in meta-analysis of
diagnostic test accuracy, see Zhou et al. (2023)
[doi:10.1002/sim.9643](https://doi.org/10.1002/sim.9643) . It also
provides the sensitivity analysis on the time-dependent SROC curve in
meta-analysis of prognostic studies, where the time-dependent
sensitivity and specificity are reconstructed from the published
survival probabilities, see Combescure et al. (2016)
[doi:10.1177/0962280212464542](https://doi.org/10.1177/0962280212464542)
.

## Details

Main functions:

- [`srocsa.t`](reference/srocsa.t.md) Sensitivity analysis for
  publication bias on the SROC curve and SAUC.

- [`tsrocsa.lrt`](reference/tsrocsa.lrt.md) Sensitivity analysis for
  publication bias on the time-dependent SROC curve and SAUC.

- [`sroc`](reference/sroc.md) Bivariate random-effects model for the
  SROC curve without publication bias.

- [`tsroc`](reference/tsroc.md) Time-dependent SROC curve without
  publication bias.

- [`sauc`](reference/sauc.md) SAUC and its confidence interval.

Results output functions:

- [`print.srocsa`](reference/print.srocsa.md) Print the main results
  from 'srocsa' objects.

- [`summary.srocsa`](reference/summary.srocsa.md) Summarize the results
  from 'srocsa' objects.

- [`plot.srocsa`](reference/plot.srocsa.md) Plot the SROC curves from
  'srocsa' objects.

Example data:

- [`IVD`](reference/IVD.md), [`CD64`](reference/CD64.md) Meta-analyses
  of diagnostic studies.

## Fitting functions

|  |  |  |
|----|----|----|
|  | **Without publication bias** | **Sensitivity analysis for publication bias** |
| Diagnostic studies (SROC) | [`reitsma.sroc`](reference/reitsma.sroc.md) | [`srocsa.t`](reference/srocsa.t.md) |
| Prognostic studies (time-dependent SROC) | [`hattori.tsroc`](reference/hattori.tsroc.md) | [`tsrocsa.lrt`](reference/tsrocsa.lrt.md) |

[`sauc`](reference/sauc.md) calculates the SAUC and its confidence
interval.

## Methods for the results

- Objects of class `"sroc"` (from
  [`reitsma.sroc`](reference/reitsma.sroc.md) and
  [`hattori.tsroc`](reference/hattori.tsroc.md)):
  [`print.sroc`](reference/print.sroc.md),
  [`summary.sroc`](reference/summary.sroc.md) (printed by
  [`print.summary.sroc`](reference/print.summary.sroc.md)), and
  [`plot.sroc`](reference/plot.sroc.md).

- Objects of class `"srocsa"` (from [`srocsa.t`](reference/srocsa.t.md)
  and [`tsrocsa.lrt`](reference/tsrocsa.lrt.md)):
  [`print.srocsa`](reference/print.srocsa.md),
  [`summary.srocsa`](reference/summary.srocsa.md), and
  [`plot.srocsa`](reference/plot.srocsa.md).

For the time-dependent results, the plot methods draw one plot for each
time point (press \<Return\> for the next in an interactive session), or
all the time points in one panel (`overlay = TRUE`).

## Example data

- [`IVD`](reference/IVD.md), [`CD64`](reference/CD64.md): meta-analyses
  of diagnostic studies.

- [`Ki67`](reference/Ki67.md): meta-analysis of Ki-67 as a prognostic
  marker in early breast cancer.

## References

Hattori S, Zhou XH (2016). Time-dependent summary receiver operating
characteristics for meta-analysis of prognostic studies. *Statistics in
Medicine*, 35(26), 4746–4763.
[doi:10.1002/sim.7029](https://doi.org/10.1002/sim.7029)

Reitsma JB, Glas AS, Rutjes AW, Scholten RJ, Bossuyt PM, Zwinderman AH
(2005). Bivariate analysis of sensitivity and specificity produces
informative summary measures in diagnostic reviews. *Journal of Clinical
Epidemiology*, 58(10), 982–990.
[doi:10.1016/j.jclinepi.2005.02.022](https://doi.org/10.1016/j.jclinepi.2005.02.022)

Zhou Y, Huang A, Hattori S (2023). A likelihood-based sensitivity
analysis for publication bias on the summary receiver operating
characteristic in meta-analysis of diagnostic test accuracy. *Statistics
in Medicine*, 42(6), 781–798.
[doi:10.1002/sim.9643](https://doi.org/10.1002/sim.9643)

Zhou Y, Huang A, Hattori S (2025). Sensitivity analysis for reporting
bias on the time-dependent summary receiver operating characteristics
curve in meta-analysis of prognosis studies with time-to-event outcomes.
*Research Synthesis Methods*, 16(3), 528–549.
[doi:10.1017/rsm.2025.14](https://doi.org/10.1017/rsm.2025.14)

## See also

Useful links:

- <https://kusuri26.github.io/DiagProgMeta/>

- <https://github.com/kusuri26/DiagProgMeta/>

- Report bugs at <https://github.com/kusuri26/DiagProgMeta/issues>

## Author

**Maintainer**: Yi Zhou <zhouy@people.kobe-u.ac.jp>
([ORCID](https://orcid.org/0000-0001-9254-3245))

Authors:

- Yi Zhou <zhouy@people.kobe-u.ac.jp>
  ([ORCID](https://orcid.org/0000-0001-9254-3245))

- Ryota Kinugawa \[contributor\]

- Kota Nakaike \[contributor\]

- Yamato Hashimoto \[contributor\]
