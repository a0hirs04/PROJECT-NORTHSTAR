from __future__ import annotations

from pathlib import Path

from python.validation.validate_biology import BiologyValidator, ScenarioResult, ScenarioSpec


def _make_result(spec: ScenarioSpec, success: bool, details: str = "ok") -> ScenarioResult:
    return ScenarioResult(
        anchor_id=spec.anchor_id,
        scenario=spec.name,
        success=success,
        criterion=spec.criterion,
        details=details,
        dependencies=list(spec.dependencies),
        run_dir="/tmp/fake",
        exit_code=0 if success else 1,
        wall_time=0.01,
        fitness=0.5 if success else 0.0,
        metrics={"total_tumor_cells": 50},
        extras={},
    )


def _make_validator(ea_test_paths: dict, tmp_path: Path) -> BiologyValidator:
    return BiologyValidator(
        binary_path=ea_test_paths["binary"],
        config_path=ea_test_paths["config"],
        output_dir=tmp_path / "validation_runs",
        timeout_seconds=5,
        sim_max_time=10.0,
    )


def test_anchor_suite_has_expected_order_and_dependencies(ea_test_paths, tmp_path):
    validator = _make_validator(ea_test_paths, tmp_path)
    specs = validator._scenario_specs()  # noqa: SLF001

    assert len(specs) == 10
    assert sorted(spec.anchor_id for spec in specs) == list(range(1, 11))
    assert specs[0].name == "ANCHOR_1_SELF_ASSEMBLY"
    assert specs[1].name == "ANCHOR_3_SHH_PARADOX"
    assert specs[-1].name == "ANCHOR_10_SPATIAL_SANCTUARY"

    by_name = {spec.name: spec for spec in specs}
    assert by_name["ANCHOR_6_PERIPHERAL_EMT"].dependencies == ["ANCHOR_5_CENTRAL_HYPOXIA"]
    assert by_name["ANCHOR_8_TWO_COMPONENT_BARRIER"].dependencies == ["ANCHOR_7_ECM_DEGRADATION_LIMITS"]
    assert by_name["ANCHOR_10_SPATIAL_SANCTUARY"].dependencies == ["ANCHOR_9_BARRIER_DENSITY_PROGNOSTIC"]


def test_dependency_gate_skips_all_downstream_when_anchor1_fails(ea_test_paths, tmp_path, monkeypatch):
    validator = _make_validator(ea_test_paths, tmp_path)
    calls: list[str] = []

    def fake_run(spec: ScenarioSpec) -> ScenarioResult:
        calls.append(spec.name)
        if spec.name == "ANCHOR_1_SELF_ASSEMBLY":
            return _make_result(spec, success=False, details="forced fail")
        return _make_result(spec, success=True)

    monkeypatch.setattr(validator, "_run_scenario", fake_run)

    summary = validator.run_all()

    assert calls == ["ANCHOR_1_SELF_ASSEMBLY"]
    assert summary.total == 10
    assert summary.passed == 0
    assert summary.failed == 10
    for result in summary.scenario_results[1:]:
        assert result.exit_code == -2
        assert "failed dependencies" in result.details


def test_dependency_gate_skips_only_affected_branches(ea_test_paths, tmp_path, monkeypatch):
    validator = _make_validator(ea_test_paths, tmp_path)
    calls: list[str] = []

    def fake_run(spec: ScenarioSpec) -> ScenarioResult:
        calls.append(spec.name)
        if spec.name == "ANCHOR_3_SHH_PARADOX":
            return _make_result(spec, success=False, details="forced fail")
        return _make_result(spec, success=True)

    monkeypatch.setattr(validator, "_run_scenario", fake_run)

    summary = validator.run_all()
    by_name = {r.scenario: r for r in summary.scenario_results}

    assert "ANCHOR_7_ECM_DEGRADATION_LIMITS" not in calls
    assert "ANCHOR_8_TWO_COMPONENT_BARRIER" not in calls

    assert "ANCHOR_2_DRUG_PENETRATION_MATURITY" in calls
    assert "ANCHOR_9_BARRIER_DENSITY_PROGNOSTIC" in calls
    assert "ANCHOR_10_SPATIAL_SANCTUARY" in calls

    assert by_name["ANCHOR_7_ECM_DEGRADATION_LIMITS"].exit_code == -2
    assert by_name["ANCHOR_8_TWO_COMPONENT_BARRIER"].exit_code == -2
    assert by_name["ANCHOR_10_SPATIAL_SANCTUARY"].success is True
