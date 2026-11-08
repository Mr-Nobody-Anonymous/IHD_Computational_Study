"""External precedent for repeated cell dosing and engraftment measurement.

This module exists because a literature search found that the architecture ATLAS-3
called "invented" is largely **already precedented**, and that two of its modelling
assumptions were wrong. Recording the sources here lets the model be corrected
against evidence rather than against intuition.

What the search changed:

1. **The attempt limit was too low.** ATLAS-3 used a precedented maximum of 2 doses.
   CATO (NCT06145035) randomises patients to **four** intravenous doses two months
   apart, and HOPE-2 gave four doses three months apart. Four is the precedented
   number, so ATLAS-3's three attempts were never beyond precedent at all.

2. **The retry operator was wrong.** ATLAS-3 modelled retries as best-of-k (a
   `maximum`). Every preclinical study in this registry instead reports **cumulative**
   benefit: each dose adds an increment, and splitting a fixed total cell number into
   three doses beats giving it all at once. Best-of-k cannot represent that result at
   all, because it makes total dose irrelevant in the wrong way.

3. **The mechanism label was wrong.** These cells do not engraft, do not become
   cardiomyocytes, and do not regenerate myocardium. The benefit is paracrine,
   anti-inflammatory and anti-fibrotic. Modelling their LVEF gain as scar replacement
   overstates what is being claimed.

4. **Human repeat-dose efficacy is not established.** The one randomised blinded
   repeat-dose trial with a functional primary endpoint (Vrtovec/REMEDIUM) was
   **null** against single dosing, and found a ceiling effect: only 10 of 30 patients
   responded to the second dose, and responders had *lower* pre-dose LVEF.

Nothing here is a result of this repository. These are other people's findings, and
the species and design of each is carried so that animal work can never earn human
credit.
"""

from __future__ import annotations

from dataclasses import dataclass

from ihd_solutions.units import require_nonnegative, require_positive

#: Design labels, ordered from strongest to weakest for efficacy claims.
HUMAN_RCT = "human_rct"
HUMAN_OPEN = "human_open"
LARGE_ANIMAL = "large_animal"
RODENT = "rodent"
REVIEW = "review"

HUMAN_DESIGNS = (HUMAN_RCT, HUMAN_OPEN)

#: What each precedent bears on.
SUPPORTS_DOSES = "attempt_limit"
SUPPORTS_CUMULATIVE = "cumulative_operator"
SUPPORTS_FEEDBACK = "feedback"
SUPPORTS_MEASUREMENT = "measurement"
SUPPORTS_CEILING = "headroom_ceiling"
SUPPORTS_MECHANISM = "mechanism"


@dataclass(frozen=True)
class Precedent:
    """One external finding, with the design that limits what it can support."""

    key: str
    citation: str
    url: str
    design: str
    population: str
    doses_tested: int
    interval_months: float
    finding: str
    supports: tuple[str, ...]
    efficacy_direction: str

    def __post_init__(self) -> None:
        if self.design not in (HUMAN_RCT, HUMAN_OPEN, LARGE_ANIMAL, RODENT, REVIEW):
            raise ValueError(f"unknown design {self.design!r}")
        require_nonnegative("doses_tested", self.doses_tested)
        require_nonnegative("interval_months", self.interval_months)
        if self.efficacy_direction not in ("positive", "null", "ongoing", "not_applicable"):
            raise ValueError(f"unknown efficacy_direction {self.efficacy_direction!r}")
        if not self.supports:
            raise ValueError("a precedent must bear on at least one lever")
        if not self.url.startswith(("https://", "http://")):
            raise ValueError("a precedent must carry a source URL")

    @property
    def is_human(self) -> bool:
        return self.design in HUMAN_DESIGNS

    @property
    def earns_human_efficacy_credit(self) -> bool:
        """Only a completed human trial with a positive result earns efficacy credit.

        An ongoing trial supports feasibility and the dose count, never efficacy.
        """
        return self.is_human and self.efficacy_direction == "positive"


PRECEDENTS: tuple[Precedent, ...] = (
    Precedent(
        key="cato",
        citation=(
            "Bolli R et al. Design and rationale of CATO, a Phase IIA randomized, "
            "double-blind, placebo-controlled study of single or repeated intravenous "
            "UC-MSCs in ischemic cardiomyopathy. Am Heart J 294 (April 2026)."
        ),
        url="https://clinicaltrials.gov/study/NCT06145035",
        design=HUMAN_RCT,
        population="ischaemic cardiomyopathy, EF <=40% by MRI, infarct >=5% LV",
        doses_tested=4,
        interval_months=2.0,
        finding=(
            "60 patients randomised 1:1:1 to four placebo, one UC-MSC + three placebo, "
            "or four UC-MSC doses of 100M cells IV. Primary endpoint change in LVEF by "
            "cardiac MRI at 12 months. 47 enrolled as of 2025-11-30; follow-up to "
            "March 2027. No efficacy result yet."
        ),
        supports=(SUPPORTS_DOSES,),
        efficacy_direction="ongoing",
    ),
    Precedent(
        key="hope2",
        citation=(
            "McDonald CM et al. Repeated intravenous cardiosphere-derived cell therapy "
            "in late-stage Duchenne muscular dystrophy (HOPE-2). Lancet 2022;399:1049-58."
        ),
        url="https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00012-5/fulltext",
        design=HUMAN_RCT,
        population="Duchenne muscular dystrophy (NOT ischaemic heart disease)",
        doses_tested=4,
        interval_months=3.0,
        finding=(
            "Four IV infusions of 150M CDCs every 3 months, randomised double-blind "
            "placebo-controlled. First repeat sequential cell dosing for a cardiac "
            "indication. Upper-limb function benefit. Repeat dosing recruits adaptive "
            "immunity, untested before this trial."
        ),
        supports=(SUPPORTS_DOSES,),
        efficacy_direction="positive",
    ),
    Precedent(
        key="vrtovec_remedium",
        citation=(
            "Vrtovec B et al. Effects of Repetitive Transendocardial CD34+ Cell "
            "Transplantation in Patients With Nonischemic Dilated Cardiomyopathy. "
            "Circ Res 2018;123:389-396."
        ),
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC6400294/",
        design=HUMAN_RCT,
        population="non-ischaemic cardiomyopathy, LVEF <40%",
        doses_tested=2,
        interval_months=6.0,
        finding=(
            "60 patients randomised to two transendocardial CD34+ doses 6 months apart "
            "or one. **No difference versus single dose** on LVEF, volumes, exercise "
            "tolerance or NT-proBNP at 1 year. Only 10 of 30 responded to the second "
            "dose (>=5 pp LVEF); responders had lower pre-dose LVEF than non-responders "
            "(34.9% vs 40.0%, p=0.02), indicating a ceiling effect."
        ),
        supports=(SUPPORTS_DOSES, SUPPORTS_CEILING, SUPPORTS_FEEDBACK),
        efficacy_direction="null",
    ),
    Precedent(
        key="yao_2011",
        citation=(
            "Yao et al. Two intracoronary BM-MNC doses 3 months apart after first STEMI "
            "with LVEF 20-39%; n=39, MRI/SPECT endpoints."
        ),
        url="https://academic.oup.com/cardiovascres/article/118/4/951/6237881",
        design=HUMAN_RCT,
        population="first STEMI, LVEF 20-39%",
        doses_tested=2,
        interval_months=3.0,
        finding=(
            "Three arms: single dose, two doses 3 months apart, placebo. Two doses "
            "beat one dose beat control on LVEF, LV volumes and scar size."
        ),
        supports=(SUPPORTS_DOSES, SUPPORTS_CUMULATIVE),
        efficacy_direction="positive",
    ),
    Precedent(
        key="gu_2011",
        citation=(
            "Gu et al. Repeated intracoronary mononuclear cell infusion at 6 months in "
            "ischaemic heart failure; n=45."
        ),
        url="https://academic.oup.com/cardiovascres/article/118/4/951/6237881",
        design=HUMAN_OPEN,
        population="ischaemic heart failure",
        doses_tested=2,
        interval_months=6.0,
        finding=(
            "Two doses 6 months apart produced a more pronounced LVEF improvement than "
            "a single dose. Not placebo-controlled."
        ),
        supports=(SUPPORTS_DOSES, SUPPORTS_CUMULATIVE),
        efficacy_direction="positive",
    ),
    Precedent(
        key="bolli_rat_2016",
        citation=(
            "Guo Y, ... Bolli R. Repeated Administrations of Cardiac Progenitor Cells "
            "Are Markedly More Effective Than a Single Administration: A New Paradigm "
            "in Cell Therapy. Circ Res 2016;119:545-557."
        ),
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC4992406/",
        design=RODENT,
        population="rat, chronic ischaemic cardiomyopathy (old MI)",
        doses_tested=3,
        interval_months=1.17,
        finding=(
            "Three doses of 12M CPCs 35 days apart versus ONE dose of the same total "
            "36M cells. Cumulative LVEF gain +13.0 +/-0.8 pp versus +4.7 +/-0.8 pp "
            "(p<0.01), roughly triple. Increments per dose +3.4, +4.0, +5.5 pp. "
            "Formation of new cardiomyocytes was negligible in all groups: the effect "
            "cannot be ascribed to differentiation. The two extra injections could not "
            "be replaced by tripling the single dose."
        ),
        supports=(SUPPORTS_CUMULATIVE, SUPPORTS_MECHANISM),
        efficacy_direction="positive",
    ),
    Precedent(
        key="bolli_rat_2018",
        citation=(
            "Repeated Administrations of Cardiac Progenitor Cells Are Superior to a "
            "Single Administration of an Equivalent Cumulative Dose. "
            "J Am Heart Assoc 2018;7:e007400."
        ),
        url="https://www.ahajournals.org/doi/10.1161/JAHA.117.007400",
        design=RODENT,
        population="rat, post-infarction LV dysfunction",
        doses_tested=3,
        interval_months=1.17,
        finding=(
            "Three 12M doses beat one 36M dose: LVEF +13.1 +/-1.3 pp versus +0.8 "
            "+/-3.1 pp. Repeated doses reduced collagen and immune infiltration; the "
            "single dose did not. New cardiomyocyte formation negligible."
        ),
        supports=(SUPPORTS_CUMULATIVE, SUPPORTS_MECHANISM),
        efficacy_direction="positive",
    ),
    Precedent(
        key="bolli_pig_2024",
        citation=(
            "Tang XL et al. Intravenous infusions of mesenchymal stromal cells have "
            "cumulative beneficial effects in a porcine model of chronic ischaemic "
            "cardiomyopathy. Cardiovasc Res 2024;120:1939-1953."
        ),
        url="https://academic.oup.com/cardiovascres/article/120/15/1939/7737527",
        design=LARGE_ANIMAL,
        population="pig, chronic ischaemic cardiomyopathy, reperfused MI",
        doses_tested=3,
        interval_months=1.17,
        finding=(
            "Three IV MSC doses 35 days apart beat one dose on three independent "
            "modalities (echo, MRI, haemodynamics); reduced fibrosis and hypertrophy. "
            "Very few IV cells reach the heart and they do not engraft, so "
            "'myocardial regeneration is not plausible'; mechanism is systemic "
            "anti-inflammatory. Direct preclinical basis for CATO."
        ),
        supports=(SUPPORTS_CUMULATIVE, SUPPORTS_MECHANISM),
        efficacy_direction="positive",
    ),
    Precedent(
        key="bolli_strategy_2019",
        citation=(
            "Bolli R, Hare JM. Potential Strategies for Clinical Translation of "
            "Repeated Cell Therapy. Circ Res 2019;124:1088-1090."
        ),
        url="https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.118.314653",
        design=REVIEW,
        population="heart failure, ischaemic and non-ischaemic",
        doses_tested=0,
        interval_months=0.0,
        finding=(
            "Explicitly proposes replacing the fixed inter-dose interval with "
            "response-guided timing: monitor serially and re-dose when the effect of "
            "the previous dose ceases, rather than at an arbitrary fixed interval. "
            "Proposed, NOT tested. This is the open loop ATLAS-3 assumed."
        ),
        supports=(SUPPORTS_FEEDBACK,),
        efficacy_direction="not_applicable",
    ),
    Precedent(
        key="pet_reporter_gene_2012",
        citation=(
            "Early Stem Cell Engraftment Predicts Late Cardiac Functional Recovery: "
            "Preclinical Insights From Molecular Imaging. "
            "Circ Cardiovasc Imaging 2012;5:391-401."
        ),
        url="https://www.ahajournals.org/doi/10.1161/circimaging.111.969329",
        design=RODENT,
        population="mouse MI, human cardiac progenitor cells, n=34-37",
        doses_tested=1,
        interval_months=0.0,
        finding=(
            "Day-1 PET reporter-gene signal predicts later function: rho=0.76 with "
            "week-2 echocardiographic function and rho=0.80 with cardiac MRI LVEF "
            "(p<0.0001). Low-engraftment animals showed no significant functional "
            "gain. This is a measured, published feedback signal quality - but in "
            "mice, with a reporter gene the cells must be engineered to carry."
        ),
        supports=(SUPPORTS_MEASUREMENT, SUPPORTS_FEEDBACK),
        efficacy_direction="not_applicable",
    ),
    Precedent(
        key="iron_mri_not_viability_2009",
        citation=(
            "Combined Reporter Gene PET and Iron Oxide MRI for Monitoring Survival and "
            "Localization of Transplanted Cells in the Rat Heart. J Nucl Med "
            "2009;50:1088-1094."
        ),
        url="https://jnm.snmjournals.org/content/50/7/1088",
        design=RODENT,
        population="rat heart, human endothelial progenitor cells",
        doses_tested=1,
        interval_months=0.0,
        finding=(
            "Iron-oxide MRI signal was unchanged for 7 days while reporter-gene PET "
            "became undetectable by day 7 and histology showed extensive apoptosis. "
            "Iron is phagocytosed by macrophages after the cells die, so iron MRI "
            "localises but does NOT measure viability. Only reporter-gene PET tracked "
            "living cells."
        ),
        supports=(SUPPORTS_MEASUREMENT,),
        efficacy_direction="not_applicable",
    ),
    Precedent(
        key="retention_not_engraftment_2017",
        citation=(
            "Non-invasive in vivo imaging of cardiac stem/progenitor cell "
            "biodistribution and retention after intracoronary and intramyocardial "
            "delivery in a swine model. J Transl Med 2017;15:57."
        ),
        url="https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-017-1157-0",
        design=LARGE_ANIMAL,
        population="minipig, chronic ischaemia-reperfusion, n=6",
        doses_tested=1,
        interval_months=0.0,
        finding=(
            "Acute retention was statistically similar by intramyocardial (13.4 +/-3.4 "
            "%ID) and intracoronary (17.4 +/-4.1 %ID) delivery, yet engrafted cells "
            "were detected histologically only after intramyocardial injection. "
            "Acute retention does NOT correlate with engraftment."
        ),
        supports=(SUPPORTS_MEASUREMENT,),
        efficacy_direction="not_applicable",
    ),
)


def precedent_frame():
    """The registry as a table, for reporting."""
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "key": p.key,
                "design": p.design,
                "population": p.population,
                "doses_tested": p.doses_tested,
                "interval_months": p.interval_months,
                "efficacy_direction": p.efficacy_direction,
                "earns_human_efficacy_credit": p.earns_human_efficacy_credit,
                "supports": "|".join(p.supports),
                "finding": p.finding,
                "citation": p.citation,
                "url": p.url,
            }
            for p in PRECEDENTS
        ]
    )


def max_precedented_doses(*, human_only: bool = True) -> int:
    """Largest number of administrations actually given in a published design.

    Defaults to human designs only. An ongoing trial counts: it establishes that the
    dose count is considered feasible and permissible, which is what a precedent
    limit means. It says nothing about efficacy.
    """
    candidates = [
        p.doses_tested for p in PRECEDENTS if p.doses_tested > 0 and (p.is_human or not human_only)
    ]
    if not candidates:
        raise ValueError("no precedent carries a dose count")
    return max(candidates)


def measurement_quality_anchor() -> tuple[float, float]:
    """Published correlation between an early engraftment signal and later function.

    Returns the (low, high) range from the PET reporter-gene study: rho=0.76 against
    echocardiography and rho=0.80 against cardiac MRI. Rodent, engineered reporter.
    ATLAS-3 used 0.70 as an unattributed guess; this replaces the guess with a source,
    and the species is carried so it cannot pass as human measurement performance.
    """
    return (0.76, 0.80)


def retry_responder_fraction() -> float:
    """Fraction responding to a second dose in the one blinded human repeat trial.

    10 of 30 patients improved LVEF by at least 5 pp after the second dose.
    """
    return 10.0 / 30.0


def cumulative_increment_ratios() -> tuple[float, ...]:
    """Per-dose LVEF increments from the rat CPC study, normalised to the first.

    Reported increments were +3.4, +4.0, +5.5 pp. Normalising gives the shape of the
    cumulative curve without importing its absolute magnitude, which is rodent and
    would be a translation claim this repository is not entitled to make.
    """
    increments = (3.4, 4.0, 5.5)
    first = require_positive("first increment", increments[0])
    return tuple(value / first for value in increments)


def human_repeat_efficacy_verdict() -> str:
    """One line, honest about what the human repeat-dose evidence actually shows."""
    human = [p for p in PRECEDENTS if p.is_human and SUPPORTS_CUMULATIVE in p.supports]
    positive = [p for p in human if p.efficacy_direction == "positive"]
    blinded_null = [
        p
        for p in PRECEDENTS
        if p.is_human and p.design == HUMAN_RCT and p.efficacy_direction == "null"
    ]
    return (
        f"{len(positive)} positive human repeat-dose reports "
        f"({', '.join(p.key for p in positive)}), but the one randomised blinded trial "
        f"with a functional primary endpoint "
        f"({', '.join(p.key for p in blinded_null)}) was NULL versus single dosing. "
        "Repeat dosing is precedented and feasible; its efficacy in humans is not "
        "established."
    )
