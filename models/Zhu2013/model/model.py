"""
PS Module – Calvin Cycle (C3 photosynthesis)
=============================================
Converted from MATLAB ePhotosynthesis (cropsinsilico/ePhotosynthesis) to
Python/mxlpy.

Source files:
  PSInitial.m  →  ps_parameters() + ps_initial_conditions()
  PSRate.m     →  rate functions + derived quantities
  PSmb.m       →  stoichiometry dicts in add_reaction()

State variables (12 dynamic; CO2/O2/NADPH fixed as parameters in standalone):
  RuBP, PGA, DPGA, T3P, ADPG, FBP, E4P, S7P, SBP, ATP, HexP, PenP

Derived (algebraic) quantities:
  DHAP, GAP        – from T3P via triose-phosphate isomerase equilibrium KE4
  ADP              – conserved adenylate: ADP = PS_C_CA − ATP
  NADP             – conserved pyridine: NADP = PS_C_CN − NADPH (NADPH fixed)
  F6P, G6P, G1P    – from hexose-phosphate lump HexP via KE21/KE22
  Ru5P, Ri5P, Xu5P – from pentose-phosphate lump PenP via KE11/KE12
  Pi               – from total-phosphate conservation + quadratic solve
  TK_Den           – shared transketolase denominator (v7 and v10)
  ATPreg           – regulator for phosphate translocator export

Reactions (18):
  v1  RuBisCO           RuBP + CO2 → 2 PGA
  v2  PGA kinase         PGA + ATP ↔ DPGA
  v3  GAPDH              DPGA + NADPH ↔ T3P
  v5  Aldolase           2 T3P ↔ FBP
  v6  FBPase             FBP → HexP
  v7  Transketolase-1    HexP(F6P) + T3P(GAP) → E4P + PenP(Xu5P)
  v8  Aldolase-2         T3P(DHAP) + E4P → SBP
  v9  SBPase             SBP → S7P
  v10 Transketolase-2    S7P + T3P(GAP) → 2 PenP(Xu5P+Ri5P)
  v13 PRK                PenP(Ru5P) + ATP → RuBP
  v16 ATP synthase       ADP + Pi → ATP
  v23 ADPG PPase         HexP(G1P) + ATP → ADPG
  v24 Starch synthase    ADPG → starch
  v25 ATP overflow       ATP → HexP (futile/overflow reaction)
  v31 PT-DHAP export     T3P(DHAP) → cytosol
  v32 PT-PGA export      PGA → cytosol
  v33 PT-GAP export      T3P(GAP) → cytosol
"""

from __future__ import annotations

import math
import numpy as np
from mxlpy import Model, Simulator


# ─────────────────────────────────────────────────────────────────────────────
# 1. PARAMETERS  (from PSInitial.m, PSRate.m)
#    All PSRatio factors = 1 (default, unperturbed model).
# ─────────────────────────────────────────────────────────────────────────────

def ps_parameters() -> dict:
    """Return all kinetic constants and Vmax values for the PS module."""
    p = {}

    # ── Conserved pool totals ─────────────────────────────────────────────
    p["PS_C_CA"]  = 1.5    # mM  total adenylates (ATP + ADP)
    p["PS_C_CP"]  = 15.0   # mM  total phosphate
    p["PS_C_CN"]  = 1.0    # mM  total pyridine nucleotide (NADP + NADPH)
    p["PS_PEXT"]  = 0.5    # mM  cytosolic Pi (external phosphate)

    # ── Fixed environmental inputs (standalone PS; dynamics switched off) ─
    p["CO2"]   = 0.012     # mM  stromal CO2
    p["O2"]    = 0.2646    # mM  stromal O2  (0.21 * 1.26)
    p["NADPH"] = 0.21      # mM  NADPH (kept constant in standalone mode)

    # ── Michaelis / inhibition / equilibrium constants ────────────────────
    # Reaction 1: RuBisCO  RuBP + CO2 → 2 PGA
    p["KM11"]  = 0.0115    # CO2
    p["KM12"]  = 0.222     # O2
    p["KM13"]  = 0.02      # RuBP
    p["KI11"]  = 0.84      # PGA   (competitive inhibition)
    p["KI12"]  = 0.04      # FBP
    p["KI13"]  = 0.075     # SBP
    p["KI14"]  = 0.9       # Pi
    p["KI15"]  = 0.07      # NADPH

    # Reaction 2: PGA kinase  PGA + ATP ↔ ADP + DPGA
    p["KM21"]  = 0.240     # PGA
    p["KM22"]  = 0.390     # ATP
    p["KM23"]  = 0.23      # ADP

    # Reaction 3: GAPDH  DPGA + NADPH ↔ T3P
    p["KM31a"] = 0.004     # DPGA
    p["KM32b"] = 0.1       # NADPH

    # Reaction 4: Triose-phosphate isomerase  DHAP ↔ GAP  (equilibrium only)
    p["KE4"]   = 1.0 / 0.05  # = 20  (→ DHAP dominant)

    # Reaction 5: Aldolase-1  GAP + DHAP ↔ FBP
    p["KM51"]  = 0.3
    p["KM52"]  = 0.4
    p["KM53"]  = 0.02
    p["KE5"]   = 7.1

    # Reaction 6: FBPase  FBP → F6P + Pi
    p["KM61"]  = 0.033
    p["KI61"]  = 0.7       # F6P  (product inhibition)
    p["KI62"]  = 12.0      # Pi
    p["KE6"]   = 6.66e5

    # Reaction 7 & 10: Transketolase (shared enzyme)
    p["KM71"]   = 0.100    # Xu5P
    p["KM72"]   = 0.100    # E4P
    p["KM73"]   = 0.1      # F6P
    p["KM74"]   = 0.100    # GAP
    p["KE7"]    = 0.1
    # Extended transketolase denominator constants (PSRate.m local, base values):
    p["KE57"]   = 1.005 * 0.1    # equilibrium for TK1 (F6P + GAP)
    p["Km8p5p"] = 0.118          # pentose-phosphate Km
    p["Km5p5p"] = 0.616
    p["KE810"]  = 0.8446         # equilibrium for TK2 (S7P + GAP)
    p["Km5gap"] = 0.2727
    p["Km8f6p"] = 0.5443
    p["Km8s7p"] = 0.01576
    p["Km8gap"] = 0.09

    # Reaction 8: Aldolase-2  DHAP + E4P ↔ SBP
    p["KM8"]   = 0.02
    p["KM81"]  = 0.4       # DHAP
    p["KM82"]  = 0.2       # E4P
    p["KE8"]   = 1.017

    # Reaction 9: SBPase  SBP → S7P + Pi
    p["KM9"]   = 0.05
    p["KI9"]   = 12.0      # Pi
    p["KE9"]   = 6.66e5

    # Reaction 10: see Reaction 7 (shared TK constants above)

    # Reaction 13: PRK  Ru5P + ATP → RuBP + ADP
    p["KM131"] = 0.05      # Ru5P
    p["KM132"] = 0.059     # ATP
    p["KI131"] = 2.0       # PGA
    p["KI132"] = 0.7       # RuBP
    p["KI133"] = 4.0       # Pi
    p["KI134"] = 2.5       # ADP
    p["KI135"] = 0.4       # ADP (second inhibition term)
    p["KE13"]  = 6.846e3

    # Reaction 16: ATP synthase  ADP + Pi → ATP
    p["KM161"] = 0.014     # ADP
    p["KM162"] = 0.3       # Pi
    p["KM163"] = 0.3       # ATP
    p["KE16"]  = 5.734

    # Reactions 21/22: hexose-phosphate pool equilibria
    p["KE21"]  = 2.3       # F6P ↔ G6P
    p["KE22"]  = 0.058     # G6P ↔ G1P

    # Reaction 23: ADPG pyrophosphorylase  G1P + ATP → ADPG + PPi
    p["KM231"] = 0.031     # G1P
    p["KM232"] = 0.045     # ATP
    p["KM233"] = 0.14      # ADPG
    p["KM234"] = 0.8       # PPi
    p["KA231"] = 0.23      # PGA activator
    p["KI231"] = 0.9       # Pi inhibitor
    p["KVmo"]  = 0.007     # minimum Vmax
    p["KE23"]  = 7.6e-3

    # Reaction 24: Starch synthase  ADPG → starch
    p["KM241"] = 0.2       # ADPG
    p["KM242"] = 0.6       # ADP
    p["KE24"]  = 7.4e5

    # Reaction 25 (ATP overflow / RuBP stabilisation)
    p["KE25"]     = 1.2e7   # used in Pi quadratic solve AND v25
    p["V25max"]   = 0.5 / 100.0 / 5.0   # 0.001
    p["MaxCoeff"] = 5.0

    # Phosphate translocator constants (reactions 31, 32, 33)
    p["KM311"] = 0.077     # DHAP
    p["KM312"] = 0.63      # Pi (stromal)
    p["KM313"] = 0.74      # Pi (external)
    p["KM32"]  = 0.25      # PGA
    p["KM33"]  = 0.075     # GAP

    # ── Vmax values ──────────────────────────────────────────────────────
    # Scaling coefficients from PSInitial.m (defaults)
    SC   = 1.0
    SC1  = 1.0
    STOM1, STOM2 = 1.0, 1.0

    p["V1"]  = 2.93  * SC1 / STOM1   # RuBisCO
    p["V2"]  = 30.15 * SC  * STOM2   # PGA kinase
    p["V3"]  = 4.04  * SC  * STOM2   # GAPDH
    p["V5"]  = 1.22  * SC            # Aldolase-1
    p["V6"]  = 0.734 * SC  / STOM1   # FBPase
    p["V7"]  = 3.12  * SC  * 4       # Transketolase (both v7 and v10)
    p["V8"]  = 1.22  * SC            # Aldolase-2
    p["V9"]  = 0.32  * 3             # SBPase
    p["V13"] = 10.81 * SC1           # PRK
    p["V16"] = 5.47                  # ATP synthase
    p["V23"] = 2.0                   # ADPG PPase
    p["V24"] = 2.0                   # Starch synthase
    p["V31"] = 1.0  * 20             # PT-DHAP
    p["V32"] = 1.0                   # PT-PGA
    p["V33"] = 1.0  * 20             # PT-GAP

    # Photorespiration coupling (zero in standalone mode)
    p["PR2PS_Pgca"] = 0.0

    return p


# ─────────────────────────────────────────────────────────────────────────────
# 2. INITIAL CONDITIONS  (from PSInitial.m, PSs vector)
# ─────────────────────────────────────────────────────────────────────────────

def ps_initial_conditions() -> dict:
    """Return initial metabolite concentrations (mM) for dynamic variables."""
    return {
        "RuBP":  2.000,   # PSs(1)
        "PGA":   2.400,   # PSs(2)
        "DPGA":  0.0011,  # PSs(3)
        "T3P":   0.5,     # PSs(4)  — lumped DHAP+GAP via KE4
        "ADPG":  0.005,   # PSs(5)
        "FBP":   0.670,   # PSs(6)
        "E4P":   0.050,   # PSs(7)
        "S7P":   2.000,   # PSs(8)
        "SBP":   0.300,   # PSs(9)
        "ATP":   0.68,    # PSs(10)
        # NADPH = 0.21  → parameter (dNADPH/dt = 0 in standalone)
        # CO2   = 0.012 → parameter
        # O2    = 0.2646 → parameter
        "HexP":  2.2,     # PSs(14)  — lumped F6P+G6P+G1P
        "PenP":  0.25,    # PSs(15)  — lumped Ru5P+Ri5P+Xu5P
    }


# ─────────────────────────────────────────────────────────────────────────────
# 3. DERIVED QUANTITY FUNCTIONS
#    These map to add_derived() calls in mxlpy.
#    Convention: all arguments are explicitly named and match add_derived args=[].
# ─────────────────────────────────────────────────────────────────────────────

def calc_DHAP(T3P, KE4):
    """DHAP from T3P lump: DHAP = T3P / (1 + KE4)"""
    return T3P / (1.0 + KE4)

def calc_GAP(T3P, KE4):
    """GAP from T3P lump: GAP = KE4 * T3P / (1 + KE4)"""
    return KE4 * T3P / (1.0 + KE4)

def calc_ADP(ATP, PS_C_CA):
    """ADP from conserved adenylate pool."""
    return PS_C_CA - ATP

def calc_NADP(NADPH, PS_C_CN):
    """NADP from conserved pyridine nucleotide pool."""
    return PS_C_CN - NADPH

# ── Hexose-phosphate sub-pool fractions ──────────────────────────────────────
def calc_F6P(HexP, KE21, KE22):
    """Fructose-6-phosphate from hexose-phosphate lump."""
    denom = 1.0 + 1.0/KE21 + KE22
    return (HexP / KE21) / denom

def calc_G6P(HexP, KE21, KE22):
    """Glucose-6-phosphate from hexose-phosphate lump."""
    denom = 1.0 + 1.0/KE21 + KE22
    return HexP / denom

def calc_G1P(HexP, KE21, KE22):
    """Glucose-1-phosphate from hexose-phosphate lump."""
    denom = 1.0 + 1.0/KE21 + KE22
    return (HexP * KE22) / denom

# ── Pentose-phosphate sub-pool fractions ─────────────────────────────────────
def calc_Ru5P(PenP, KE11, KE12):
    """Ribulose-5-phosphate from pentose-phosphate lump."""
    denom = 1.0 + 1.0/KE11 + 1.0/KE12
    return PenP / denom

def calc_Ri5P(PenP, KE11, KE12):
    """Ribose-5-phosphate from pentose-phosphate lump."""
    denom = 1.0 + 1.0/KE11 + 1.0/KE12
    return (PenP / KE11) / denom

def calc_Xu5P(PenP, KE11, KE12):
    """Xylulose-5-phosphate from pentose-phosphate lump."""
    denom = 1.0 + 1.0/KE11 + 1.0/KE12
    return (PenP / KE12) / denom

# ── Pi from total-phosphate conservation (quadratic) ─────────────────────────
def calc_Pi(
    RuBP, PGA, DPGA, T3P, FBP, E4P, SBP, S7P, ATP,
    HexP, PenP,
    KE4, KE21, KE22, KE11, KE12, KE25,
    PS_C_CP, PR2PS_Pgca,
):
    """
    Free Pi via total-phosphate balance + quadratic PPi equilibrium.

    Pit = PS_C_CP − (sum of all phosphorylated species)
    Pi  = 0.5 * (−KE25 + sqrt(KE25² + 4·Pit·KE25))
    """
    DHAP = T3P / (1.0 + KE4)
    GAP  = KE4 * T3P / (1.0 + KE4)
    denom_hex = 1.0 + 1.0/KE21 + KE22
    F6P = (HexP / KE21) / denom_hex
    G6P = HexP / denom_hex
    G1P = (HexP * KE22) / denom_hex
    denom_pen = 1.0 + 1.0/KE11 + 1.0/KE12
    Ru5P = PenP / denom_pen
    Ri5P = (PenP / KE11) / denom_pen
    Xu5P = (PenP / KE12) / denom_pen

    Pit = (PS_C_CP
           - PGA - 2.0*DPGA - GAP - DHAP
           - 2.0*FBP - F6P - E4P
           - 2.0*SBP - S7P
           - Xu5P - Ri5P - Ru5P
           - 2.0*RuBP - G6P - G1P
           - ATP - PR2PS_Pgca)

    Pi = 0.5 * (-KE25 + math.sqrt(KE25*KE25 + 4.0*Pit*KE25))
    return max(Pi, 1e-9)   # guard against negative values

# ── Transketolase shared denominator ─────────────────────────────────────────
def calc_TK_Den(GAP, F6P, S7P, Xu5P, E4P, Ri5P,
                Km5gap, Km8f6p, Km8s7p, Km8gap, Km8p5p, Km5p5p):
    """
    Shared denominator for transketolase reactions v7 and v10.
    Den = 1 + (1 + GAP/Km5gap)*(F6P/Km8f6p + S7P/Km8s7p)
            + GAP/Km8gap
            + (1/Km8p5p)*(Xu5P*(1 + E4P*Ri5P/Km5p5p) + E4P + Ri5P)
    """
    den = (1.0
           + (1.0 + GAP/Km5gap) * (F6P/Km8f6p + S7P/Km8s7p)
           + GAP/Km8gap
           + (1.0/Km8p5p) * (Xu5P * (1.0 + E4P*Ri5P/Km5p5p) + E4P + Ri5P))
    return den

# ── Regulation factor for phosphate translocators ────────────────────────────
def calc_ATPreg(PGA):
    """
    ATPreg = min(PGA / 3, 1.0)
    Scales phosphate translocator export rates.
    (RedoxReg_RA_com == 0 branch from PSRate.m)
    """
    return min(PGA / 3.0, 1.0)


# ─────────────────────────────────────────────────────────────────────────────
# 4. RATE FUNCTIONS  (from PSRate.m)
# ─────────────────────────────────────────────────────────────────────────────

# ── v1: RuBisCO (RUBISCOMETHOD=2, the default) ───────────────────────────────
def v1_rubisco(RuBP, CO2, O2, PGA, FBP, SBP, Pi, NADPH,
               V1, KM11, KM12, KM13, KI11, KI12, KI13, KI14, KI15):
    V1Reg = 1.0 + PGA/KI11 + FBP/KI12 + SBP/KI13 + Pi/KI14 + NADPH/KI15
    sat   = RuBP / (RuBP + KM13 * V1Reg)
    v1    = V1 * sat * CO2 / (CO2 + KM11 * (1.0 + O2/KM12))
    threshold = V1 / 2.5
    if RuBP < threshold:
        v1 = v1 * RuBP / threshold
    return v1

# ── v2: PGA kinase ────────────────────────────────────────────────────────────
def v2_pga_kinase(PGA, ATP, ADP, V2, KM21, KM22, KM23):
    return V2 * PGA * ATP / ((PGA + KM21) * (ATP + KM22 * (1.0 + ADP/KM23)))

# ── v3: GAPDH ─────────────────────────────────────────────────────────────────
def v3_gapdh(DPGA, NADPH, V3, KM31a, KM32b):
    return V3 * DPGA * NADPH / ((DPGA + KM31a) * (NADPH + KM32b))

# ── v5: Aldolase-1  GAP + DHAP ↔ FBP ────────────────────────────────────────
def v5_aldolase1(GAP, DHAP, FBP, V5, KM51, KM52, KM53, KE5):
    num = GAP * DHAP - FBP / KE5
    den = (KM51 * KM52) * (1.0 + GAP/KM51 + DHAP/KM52 + FBP/KM53
                            + GAP*DHAP/(KM51*KM52))
    return V5 * num / den

# ── v6: FBPase  FBP → F6P + Pi ──────────────────────────────────────────────
def v6_fbpase(FBP, F6P, Pi, V6, KM61, KI61, KI62, KE6):
    num = FBP - F6P * Pi / KE6
    den = FBP + KM61 * (1.0 + F6P/KI61 + Pi/KI62)
    return V6 * num / den

# ── v7: Transketolase-1  F6P + GAP → E4P + Xu5P ─────────────────────────────
def v7_tk1(F6P, GAP, E4P, Xu5P, TK_Den, V7, KE57, Km8p5p, Km5p5p):
    num = F6P * GAP * KE57 - E4P * Xu5P
    return V7 * num / (Km8p5p * Km5p5p * TK_Den)

# ── v8: Aldolase-2  DHAP + E4P → SBP ────────────────────────────────────────
def v8_aldolase2(DHAP, E4P, SBP, V8, KM81, KM82, KE8):
    num = DHAP * E4P - SBP / KE8
    den = (E4P + KM82) * (DHAP + KM81)
    return V8 * num / den

# ── v9: SBPase  SBP → S7P + Pi ──────────────────────────────────────────────
def v9_sbpase(SBP, S7P, Pi, V9, KM9, KI9, KE9):
    num = SBP - Pi * S7P / KE9
    den = SBP + KM9 * (1.0 + Pi/KI9)
    return V9 * num / den

# ── v10: Transketolase-2  S7P + GAP → Ri5P + Xu5P ───────────────────────────
def v10_tk2(S7P, GAP, Xu5P, Ri5P, TK_Den, V7, KE810, Km8p5p, Km5p5p):
    num = S7P * GAP * KE810 - Xu5P * Ri5P
    return V7 * num / (Km8p5p * Km5p5p * TK_Den)

# ── v13: PRK  Ru5P + ATP → RuBP + ADP ───────────────────────────────────────
def v13_prk(Ru5P, ATP, ADP, RuBP, Pi, PGA,
            V13, KM131, KM132, KI131, KI132, KI133, KI134, KI135, KE13):
    num = ATP * Ru5P - ADP * RuBP / KE13
    den_atp = ATP * (1.0 + ADP/KI134) + KM132 * (1.0 + ADP/KI135)
    den_ru5 = Ru5P + KM131 * (1.0 + PGA/KI131 + RuBP/KI132 + Pi/KI133)
    return V13 * num / (den_atp * den_ru5)

# ── v16: ATP synthase  ADP + Pi → ATP ────────────────────────────────────────
def v16_atp_synthase(ADP, Pi, ATP, V16, KM161, KM162, KM163, KE16):
    num = ADP * Pi - ATP / KE16
    den = (KM161 * KM162
           * (1.0 + ADP/KM161 + Pi/KM162 + ATP/KM163
              + ADP*Pi/(KM161*KM162)))
    return V16 * num / den

# ── v23: ADPG pyrophosphorylase  G1P + ATP → ADPG + PPi ─────────────────────
def v23_adpg_ppase(G1P, ATP, ADPG, Pi, PGA,
                   V23, KVmo, KA231, KI231, KM231, KM232, KM233, KM234, KE23):
    # Assume PPi ≈ Pi (Pit − Pi); use Pi as proxy for OPOP
    OPOP = max(Pi * 0.1, 1e-9)   # simplified; full model tracks PPi separately
    Va = KVmo + V23 * (PGA / (KA231 * (1.0 + PGA/KA231)))
    num = Va * (ATP * G1P - ADPG * OPOP / KE23)
    den = ((1.0 + Pi/KI231)
           * KM231 * KM232
           * (1.0 + ATP/KM232 + G1P/KM231 + ATP*G1P/(KM231*KM232)
              + ADPG/KM233 + OPOP/KM234 + ADPG*OPOP/(KM233*KM234)))
    return num / den

# ── v24: Starch synthase  ADPG → starch ──────────────────────────────────────
def v24_starch(ADPG, V24, KM241):
    return V24 * ADPG / (KM241 * (1.0 + ADPG/KM241))

# ── v25: ATP overflow / RuBP stabilisation ───────────────────────────────────
def v25_overflow(RuBP, ATP, V25max, MaxCoeff):
    return V25max * (1.0 - RuBP/MaxCoeff) * ATP / (ATP + 1.0)

# ── v31: Phosphate translocator DHAP export ──────────────────────────────────
def v31_pt_dhap(DHAP, ATPreg, V31, KM311, KM313, PS_PEXT):
    return V31 * DHAP / (DHAP + KM311) * PS_PEXT / (PS_PEXT + KM313) * ATPreg

# ── v32: Phosphate translocator PGA export ───────────────────────────────────
def v32_pt_pga(PGA, ATPreg, V32, KM32, KM313, PS_PEXT):
    return V32 * PGA / (PGA + KM32) * PS_PEXT / (PS_PEXT + KM313) * ATPreg

# ── v33: Phosphate translocator GAP export ───────────────────────────────────
def v33_pt_gap(GAP, ATPreg, V33, KM33, KM313, PS_PEXT):
    return V33 * GAP / (GAP + KM33) * PS_PEXT / (PS_PEXT + KM313) * ATPreg


# ─────────────────────────────────────────────────────────────────────────────
# 5. MODEL ASSEMBLY
# ─────────────────────────────────────────────────────────────────────────────

def build_ps_model() -> Model:
    """
    Assemble and return the standalone PS (Calvin cycle) mxlpy Model.

    Stoichiometry is derived from PSmb.m:
      dRuBP  =  v13 − v1
      dPGA   =  2·v1 − v2 − v32
      dDPGA  =  v2 − v3
      dT3P   =  v3 − 2·v5 − v7 − v8 − v10 − v31 − v33
      dADPG  =  v23 − v24
      dFBP   =  v5 − v6
      dE4P   =  v7 − v8
      dS7P   =  v9 − v10
      dSBP   =  v8 − v9
      dATP   =  v16 − v2 − v23 − v13 − v25
      dHexP  =  v6 − v7 − v23 + v25
      dPenP  =  v7 + 2·v10 − v13
    """
    params = ps_parameters()
    inits  = ps_initial_conditions()

    model = Model()

    # ── Parameters ────────────────────────────────────────────────────────
    model.add_parameters(params)

    # Extra equilibrium constants needed in PenP / HexP derived quantities
    # (not in PSInitial.m; using standard values)
    model.add_parameter("KE11", 0.4)    # Ri5P ↔ Ru5P
    model.add_parameter("KE12", 0.67)   # Xu5P ↔ Ru5P

    # ── Dynamic state variables ───────────────────────────────────────────
    model.add_variables(inits)

    # ── Derived (algebraic) quantities ───────────────────────────────────
    model.add_derived("DHAP",    fn=calc_DHAP,    args=["T3P", "KE4"])
    model.add_derived("GAP",     fn=calc_GAP,     args=["T3P", "KE4"])
    model.add_derived("ADP",     fn=calc_ADP,     args=["ATP", "PS_C_CA"])
    model.add_derived("NADP",    fn=calc_NADP,    args=["NADPH", "PS_C_CN"])
    model.add_derived("F6P",     fn=calc_F6P,     args=["HexP", "KE21", "KE22"])
    model.add_derived("G6P",     fn=calc_G6P,     args=["HexP", "KE21", "KE22"])
    model.add_derived("G1P",     fn=calc_G1P,     args=["HexP", "KE21", "KE22"])
    model.add_derived("Ru5P",    fn=calc_Ru5P,    args=["PenP", "KE11", "KE12"])
    model.add_derived("Ri5P",    fn=calc_Ri5P,    args=["PenP", "KE11", "KE12"])
    model.add_derived("Xu5P",    fn=calc_Xu5P,    args=["PenP", "KE11", "KE12"])

    model.add_derived("Pi", fn=calc_Pi, args=[
        "RuBP", "PGA", "DPGA", "T3P", "FBP", "E4P", "SBP", "S7P", "ATP",
        "HexP", "PenP",
        "KE4", "KE21", "KE22", "KE11", "KE12", "KE25",
        "PS_C_CP", "PR2PS_Pgca",
    ])

    model.add_derived("TK_Den", fn=calc_TK_Den, args=[
        "GAP", "F6P", "S7P", "Xu5P", "E4P", "Ri5P",
        "Km5gap", "Km8f6p", "Km8s7p", "Km8gap", "Km8p5p", "Km5p5p",
    ])

    model.add_derived("ATPreg", fn=calc_ATPreg, args=["PGA"])

    # ── Reactions ─────────────────────────────────────────────────────────

    # v1: RuBisCO  RuBP + CO2 → 2 PGA
    model.add_reaction("v1", fn=v1_rubisco,
        args=["RuBP", "CO2", "O2", "PGA", "FBP", "SBP", "Pi", "NADPH",
              "V1", "KM11", "KM12", "KM13", "KI11", "KI12", "KI13", "KI14", "KI15"],
        stoichiometry={"RuBP": -1, "PGA": 2},
    )

    # v2: PGA kinase  PGA + ATP → DPGA
    model.add_reaction("v2", fn=v2_pga_kinase,
        args=["PGA", "ATP", "ADP", "V2", "KM21", "KM22", "KM23"],
        stoichiometry={"PGA": -1, "ATP": -1, "DPGA": 1},
    )

    # v3: GAPDH  DPGA → T3P
    model.add_reaction("v3", fn=v3_gapdh,
        args=["DPGA", "NADPH", "V3", "KM31a", "KM32b"],
        stoichiometry={"DPGA": -1, "T3P": 1},
    )

    # v5: Aldolase-1  2 T3P → FBP
    model.add_reaction("v5", fn=v5_aldolase1,
        args=["GAP", "DHAP", "FBP", "V5", "KM51", "KM52", "KM53", "KE5"],
        stoichiometry={"T3P": -2, "FBP": 1},
    )

    # v6: FBPase  FBP → HexP (+Pi released, tracked via Pi derived)
    model.add_reaction("v6", fn=v6_fbpase,
        args=["FBP", "F6P", "Pi", "V6", "KM61", "KI61", "KI62", "KE6"],
        stoichiometry={"FBP": -1, "HexP": 1},
    )

    # v7: Transketolase-1  HexP(F6P) + T3P(GAP) → E4P + PenP(Xu5P)
    model.add_reaction("v7", fn=v7_tk1,
        args=["F6P", "GAP", "E4P", "Xu5P", "TK_Den", "V7", "KE57", "Km8p5p", "Km5p5p"],
        stoichiometry={"HexP": -1, "T3P": -1, "E4P": 1, "PenP": 1},
    )

    # v8: Aldolase-2  T3P(DHAP) + E4P → SBP
    model.add_reaction("v8", fn=v8_aldolase2,
        args=["DHAP", "E4P", "SBP", "V8", "KM81", "KM82", "KE8"],
        stoichiometry={"T3P": -1, "E4P": -1, "SBP": 1},
    )

    # v9: SBPase  SBP → S7P (+Pi released)
    model.add_reaction("v9", fn=v9_sbpase,
        args=["SBP", "S7P", "Pi", "V9", "KM9", "KI9", "KE9"],
        stoichiometry={"SBP": -1, "S7P": 1},
    )

    # v10: Transketolase-2  S7P + T3P(GAP) → 2 PenP
    model.add_reaction("v10", fn=v10_tk2,
        args=["S7P", "GAP", "Xu5P", "Ri5P", "TK_Den", "V7", "KE810", "Km8p5p", "Km5p5p"],
        stoichiometry={"S7P": -1, "T3P": -1, "PenP": 2},
    )

    # v13: PRK  PenP(Ru5P) + ATP → RuBP + ADP
    model.add_reaction("v13", fn=v13_prk,
        args=["Ru5P", "ATP", "ADP", "RuBP", "Pi", "PGA",
              "V13", "KM131", "KM132", "KI131", "KI132", "KI133", "KI134", "KI135", "KE13"],
        stoichiometry={"PenP": -1, "ATP": -1, "RuBP": 1},
    )

    # v16: ATP synthase  ADP + Pi → ATP
    model.add_reaction("v16", fn=v16_atp_synthase,
        args=["ADP", "Pi", "ATP", "V16", "KM161", "KM162", "KM163", "KE16"],
        stoichiometry={"ATP": 1},
    )

    # v23: ADPG PPase  G1P + ATP → ADPG + PPi
    model.add_reaction("v23", fn=v23_adpg_ppase,
        args=["G1P", "ATP", "ADPG", "Pi", "PGA",
              "V23", "KVmo", "KA231", "KI231",
              "KM231", "KM232", "KM233", "KM234", "KE23"],
        stoichiometry={"HexP": -1, "ATP": -1, "ADPG": 1},
    )

    # v24: Starch synthase  ADPG → starch (sink)
    model.add_reaction("v24", fn=v24_starch,
        args=["ADPG", "V24", "KM241"],
        stoichiometry={"ADPG": -1},
    )

    # v25: ATP overflow  ATP → HexP
    model.add_reaction("v25", fn=v25_overflow,
        args=["RuBP", "ATP", "V25max", "MaxCoeff"],
        stoichiometry={"ATP": -1, "HexP": 1},
    )

    # v31: PT-DHAP export  T3P(DHAP) → cytosol
    model.add_reaction("v31", fn=v31_pt_dhap,
        args=["DHAP", "ATPreg", "V31", "KM311", "KM313", "PS_PEXT"],
        stoichiometry={"T3P": -1},
    )

    # v32: PT-PGA export  PGA → cytosol
    model.add_reaction("v32", fn=v32_pt_pga,
        args=["PGA", "ATPreg", "V32", "KM32", "KM313", "PS_PEXT"],
        stoichiometry={"PGA": -1},
    )

    # v33: PT-GAP export  T3P(GAP) → cytosol
    model.add_reaction("v33", fn=v33_pt_gap,
        args=["GAP", "ATPreg", "V33", "KM33", "KM313", "PS_PEXT"],
        stoichiometry={"T3P": -1},
    )

    return model


# ─────────────────────────────────────────────────────────────────────────────
# 6. QUICK-RUN HELPER
# ─────────────────────────────────────────────────────────────────────────────

def simulate_ps(t_end: float = 200.0, n_points: int = 500):
    """Run PS standalone simulation and return (variables, fluxes) DataFrames."""
    model = build_ps_model()
    time_points = np.linspace(0, t_end, n_points)
    sim = Simulator(model).simulate_time_course(time_points)
    result = sim.get_result().value
    return result.variables, result.fluxes


if __name__ == "__main__":
    vars_, fluxes = simulate_ps(t_end=300.0, n_points=600)
    print("=== Final concentrations (t=300 s) ===")
    print(vars_.iloc[-1].to_string())
    print("\n=== Final fluxes (t=300 s) ===")
    print(fluxes.iloc[-1].to_string())
