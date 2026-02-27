# Implementation Plan: Check 3 + Anchor 8 Redesign + 5 Reality Checks

## Overview
Three targeted changes to `validate_biology.py` and supporting C++/XML:

1. **Check 3** — three-way vismodegib rank (virtual, reads from cache)
2. **Anchor 8 redesign** — HA_DEGRADE / COL_DEGRADE component-specific ECM degradation
3. **5 Reality Checks** — integration gate (all virtual, EA cannot start until all 5 pass)

---

## Part 1 — C++ + XML: ha_degrade_strength / col_degrade_strength

### Why
Anchor 8 needs independent degradation of the HA (diffusion) and collagen (mechanical) ECM fractions.
The existing knob partition blocks gene-network interventions (HAS2/COL1A1). Instead, two scalar
user-parameters drive a weighted decomposition of the existing barrier formulas — no new BioFVM
substrate required.

### 1a. `Stroma_world/PhysiCell/config/PhysiCell_settings.xml`
Add two user_parameters (default 0.0):
```xml
<ha_degrade_strength type="double" units="dimensionless">0.0</ha_degrade_strength>
<col_degrade_strength type="double" units="dimensionless">0.0</col_degrade_strength>
```

### 1b. `custom_modules/tumor_calibration_knobs.h`
Add two fields to `TumorCalibrationKnobs` (environmental, not EA-targetable):
```cpp
// --- ECM component degradation (set by external physical intervention, not EA) ---
double ha_degrade_strength  = 0.0;   // 0=no HA depletion → 1=full HA depletion
double col_degrade_strength = 0.0;   // 0=no collagen depletion → 1=full collagen depletion
```

### 1c. `custom_modules/tumor_calibration_knobs.cpp`
Add a free function `load_ecm_degradation_params(TumorCalibrationKnobs&)` that reads from
PhysiCell's parameter store:
```cpp
void load_ecm_degradation_params(TumorCalibrationKnobs& knobs)
{
    if (parameters.doubles.find_index("ha_degrade_strength") >= 0)
        knobs.ha_degrade_strength  = parameters.doubles("ha_degrade_strength");
    if (parameters.doubles.find_index("col_degrade_strength") >= 0)
        knobs.col_degrade_strength = parameters.doubles("col_degrade_strength");
}
```
Declare the function in the .h file.

### 1d. `custom_modules/custom.cpp`
Call `load_ecm_degradation_params(g_active_knobs)` after the calibration profile is loaded in
`setup_tissue()`.

### 1e. `custom_modules/tumor_cell.cpp`
Replace the two barrier formula lines with weight-based decomposition.
HA weight in diffusion barrier = 0.70; collagen weight = 0.30.
Collagen weight in mechanical stiffness = 0.80; HA weight = 0.20.

**Drug delivery barrier** (current: `1.0 - 0.5 * ecm`):
```cpp
const double ha_frac   = 1.0 - 0.70 * g_active_knobs.ha_degrade_strength;
const double col_frac  = 1.0 - 0.30 * g_active_knobs.col_degrade_strength;
ecm_delivery_barrier = 1.0 - 0.5 * local_ecm_val * (0.70 * ha_frac + 0.30 * col_frac);
```
Wait — this doesn't simplify cleanly. Correct formulation:
```cpp
// HA contributes 70% of diffusion barrier; collagen contributes 30%
const double ecm_barrier_ha  = 0.70 * (1.0 - g_active_knobs.ha_degrade_strength);
const double ecm_barrier_col = 0.30 * (1.0 - g_active_knobs.col_degrade_strength);
const double ecm_barrier_total = ecm_barrier_ha + ecm_barrier_col;   // in [0,1]
ecm_delivery_barrier = 1.0 - 0.5 * local_ecm_val * ecm_barrier_total;
```
Default (both=0): `ecm_barrier_total = 1.0` → same as current. ✓
Full HA depletion (ha=1): `total = 0.30` → barrier = `1 - 0.15*ecm` (70% reduction). ✓
Full COL depletion (col=1): `total = 0.70` → barrier = `1 - 0.35*ecm` (30% reduction). ✓

**Solid stress** (current: `simple_pressure * (1 + ecm)`):
```cpp
// Collagen contributes 80% of mechanical stiffness; HA contributes 20%
const double mech_col = 0.80 * (1.0 - g_active_knobs.col_degrade_strength);
const double mech_ha  = 0.20 * (1.0 - g_active_knobs.ha_degrade_strength);
solid_stress = simple_pressure * (1.0 + local_ecm_val * (mech_col + mech_ha));
```
Default: `mech_col + mech_ha = 1.0` → same as current. ✓
Full COL depletion: `= 0.20` → 80% reduction in mechanical amplification. ✓
Full HA depletion: `= 0.80` → 20% reduction in mechanical amplification. ✓

Independence is guaranteed by construction: HA depleted improves penetration more than COL depleted;
COL depleted releases confinement more than HA depleted.

---

## Part 2 — Python: validate_biology.py

### 2a. ScenarioSpec dataclass — add `virtual` field
```python
@dataclass
class ScenarioSpec:
    ...
    virtual: bool = False   # True = no simulation; evaluator reads from baseline_cache only
```

### 2b. run_all() — route virtual scenarios differently
```python
if spec.virtual:
    result = self._run_virtual_check(spec)
else:
    result = self._run_scenario(spec)
```

### 2c. New `_run_virtual_check(spec)` method
Calls `spec.evaluator({}, {}, self._baseline_cache)` with empty metrics/extras.
Returns a `ScenarioResult` with `exit_code=0` (or -1 on evaluator exception).
Wall time = 0.0.

### 2d. New `_cached_extras(baseline, name)` static method
```python
@staticmethod
def _cached_extras(baseline, name):
    e = baseline.get(name, {}).get("extras", {})
    return e if isinstance(e, dict) else {}
```

### 2e. ValidationSummary — add check fields
```python
@dataclass
class ValidationSummary:
    ...
    checks_total: int
    checks_passed: int
    checks_failed: int
```
`run_all()` populates these by inspecting `extras.get("check_tests", [])` (separate key from `anchor_tests`).

### 2f. Anchor 8 restructure
Remove `ANCHOR_8_COLLAGEN_DEPLETED_PROXY` and `ANCHOR_8_TWO_COMPONENT_BARRIER`.
Add three new specs (all with `anchor_id=8`):

**ANCHOR_8_HA_DEPLETED_DRUG** — HA depletion + drug, data collection only (0 anchor_tests):
- `user_parameter_overrides = {**cytotoxic_user, "ha_degrade_strength": 0.9}`
- `variable_overrides = cytotoxic_var`
- `dependencies = ["ANCHOR_7_ECM_DEGRADATION_LIMITS"]`
- Evaluator stores proxy marker, returns True

**ANCHOR_8_COL_DEPLETED_DRUG** — collagen depletion + drug, data collection only (0 anchor_tests):
- `user_parameter_overrides = {**cytotoxic_user, "col_degrade_strength": 0.9}`
- `variable_overrides = cytotoxic_var`
- `dependencies = ["ANCHOR_7_ECM_DEGRADATION_LIMITS"]`
- Evaluator stores proxy marker, returns True

**ANCHOR_8_BOTH_DEPLETED_DRUG** — both + drug, independence test (4 anchor_tests P8a–P8d):
- `user_parameter_overrides = {**cytotoxic_user, "ha_degrade_strength": 0.9, "col_degrade_strength": 0.9}`
- `variable_overrides = cytotoxic_var`
- `dependencies = ["ANCHOR_8_HA_DEPLETED_DRUG", "ANCHOR_8_COL_DEPLETED_DRUG"]`
- Evaluator reads HA arm and COL arm from cache, runs 4 independence tests:
  - P8a: `pen_ha > pen_col` (HA depletion improves penetration more than collagen depletion)
  - P8b: `extent_col > extent_ha` (collagen depletion releases confinement more than HA depletion)
  - P8c: `pen_both > pen_ha` (combined depletion better than HA alone)
  - P8d: `pen_both > pen_col` (combined depletion better than COL alone)

Anchor test count stays at 34 (0+0+4 = same as old 0+4).

ANCHOR_10_SPATIAL_SANCTUARY dependencies: replace `"ANCHOR_8_TWO_COMPONENT_BARRIER"` with
`"ANCHOR_8_BOTH_DEPLETED_DRUG"`.

### 2g. 5 Reality Checks — all virtual, appended after ANCHOR_10

**CHECK_1_NATURAL_HISTORY** (virtual, anchor_id=11, depends on A1+A5+A6):
Evaluator reads extras from cache. Single check_test:
- Verifies: CAF activated (A1: activated_caf_fraction > 0), ECM gradient (A1: ecm_peri>ecm_boundary),
  central hypoxia (A5: hif1a_core > hif1a_periphery), peripheral EMT (A6: emt_peripheral > emt_core).
Returns single PASS/FAIL stored in `extras["check_tests"]`.

**CHECK_2_CHEMO_FAILURE** (virtual, anchor_id=12, depends on A2+A7):
Evaluator reads from cache. Single check_test:
- Drug penetration declined with barrier maturity (A2: P2a/P2b passed), residual survivors persist
  (A7: viability > 0), barrier rebuilds after degradation (A7: ecm_delta > 0).

**CHECK_3_VISMODEGIB_PARADOX** (virtual, anchor_id=13, depends on A1+A3+A7):
Evaluator reads tumor_count from cached metrics for A1, A3, A7. Single check_test:
- `total_tumor_cells(A7) < total_tumor_cells(A1) < total_tumor_cells(A3)`
  i.e. SHH+drug kills most, control in middle, SHH-only worst.

**CHECK_4_FITNESS_RANKING** (virtual, anchor_id=14, depends on A2+A7+A8_BOTH):
Evaluator reads fitness scores and viabilities from cache. Single check_test:
- `viability(A8_BOTH: combination) < viability(A2: drug-only) < 1.0`
  Combination strategy outperforms single-agent; neither produces complete eradication.

**CHECK_5_SANCTUARY_REGROWTH** (virtual, anchor_id=15,
 depends on CHECK_1+CHECK_2+CHECK_3+CHECK_4+A10):
Gate check: only runs when Checks 1-4 pass.
Reads A10 extras from cache. Single check_test:
- Survivor localization (live_mean_ecm > dead_mean_ecm), distance shielding, regrowth after nadir,
  barrier re-emergence (ecm_late > ecm_at_nadir). All 4 must hold simultaneously.

### 2h. Update 34-test warning
Change `tests_total != 34` check to accommodate stable anchor count but not flag on check tests.
The check `if skipped_count == 0 and tests_total != 34` stays; Reality Check tests go in
`extras["check_tests"]` not `anchor_tests`, so they don't affect this count.

---

## Part 3 — Python: tests/test_validation_anchors.py

### 3a. test_anchor_suite_has_expected_order_and_dependencies
- `len(specs) == 17` (11 - 2 + 3 + 5)
- `set(spec.anchor_id for spec in specs) == set(range(1, 11)) | {11, 12, 13, 14, 15}`
- `specs[-1].name == "CHECK_5_SANCTUARY_REGROWTH"`
- `by_name["ANCHOR_8_BOTH_DEPLETED_DRUG"].dependencies == ["ANCHOR_8_HA_DEPLETED_DRUG", "ANCHOR_8_COL_DEPLETED_DRUG"]`
- `"ANCHOR_8_BOTH_DEPLETED_DRUG" in by_name["ANCHOR_10_SPATIAL_SANCTUARY"].dependencies`

### 3b. test_dependency_gate_skips_all_downstream_when_anchor1_fails
- Monkeypatch both `_run_scenario` AND `_run_virtual_check`
- `summary.total == 17`, `summary.passed == 0`, `summary.failed == 17`

### 3c. test_dependency_gate_skips_only_affected_branches
- `"ANCHOR_8_BOTH_DEPLETED_DRUG" not in calls` (was TWO_COMPONENT_BARRIER)
- `"ANCHOR_8_HA_DEPLETED_DRUG" not in calls` (also downstream of A7)
- `"ANCHOR_8_COL_DEPLETED_DRUG" not in calls` (also downstream of A7)
- Reality checks referencing A3 or A7 also skipped (CHECK_2, CHECK_3, CHECK_4, CHECK_5)

---

## Execution Order
1. C++ + XML changes (5 files)
2. validate_biology.py (2f Anchor 8 first, then 2a–2e infrastructure, then 2g Reality Checks)
3. test_validation_anchors.py
4. Build: `PHYSICELL_DIR=.../Stroma_world/PhysiCell make -j$(nproc) -C /home/a0hirs04/PROJECT-NORTHSTAR`
