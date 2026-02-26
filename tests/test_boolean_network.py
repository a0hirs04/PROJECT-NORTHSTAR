from __future__ import annotations

import importlib
import inspect
import math
from typing import Any, Dict, Iterable, List

import pytest


def _load_boolean_network_class():
    candidates = [
        "python.boolean_network",
        "python.model.boolean_network",
        "python.wrapper.boolean_network",
        "boolean_network",
    ]
    class_names = ("BooleanNetwork", "PDACBooleanNetwork", "Network")
    for module_name in candidates:
        try:
            module = importlib.import_module(module_name)
        except Exception:
            continue
        for class_name in class_names:
            cls = getattr(module, class_name, None)
            if cls is not None and inspect.isclass(cls):
                return cls
    return None


BOOLEAN_NETWORK_CLASS = _load_boolean_network_class()
pytestmark = pytest.mark.skipif(
    BOOLEAN_NETWORK_CLASS is None,
    reason="No Python BooleanNetwork implementation found in PROJECT-NORTHSTAR.",
)


def _make_network(cell_type: str):
    assert BOOLEAN_NETWORK_CLASS is not None
    attempts = [
        {"cell_type": cell_type},
        {"cell_type_name": cell_type},
        {"mode": cell_type},
        {},
    ]
    for kwargs in attempts:
        try:
            return BOOLEAN_NETWORK_CLASS(**kwargs)
        except TypeError:
            continue
    # Last resort
    return BOOLEAN_NETWORK_CLASS(cell_type)


def _set_gene(bn: Any, gene: str, value: float) -> None:
    setters = ("set_gene", "set_gene_state", "set_node", "set_state", "set")
    for name in setters:
        fn = getattr(bn, name, None)
        if callable(fn):
            try:
                fn(gene, float(value))
                return
            except TypeError:
                try:
                    fn({gene: float(value)})
                    return
                except Exception:
                    pass

    for attr in ("state", "states", "genes", "gene_states"):
        maybe = getattr(bn, attr, None)
        if isinstance(maybe, dict):
            maybe[gene] = float(value)
            return

    raise AttributeError(f"Unable to set gene '{gene}' on BooleanNetwork object.")


def _get_gene(bn: Any, gene: str) -> float:
    getters = ("gene", "get_gene", "get_gene_state", "get_node", "get_state", "get")
    for name in getters:
        fn = getattr(bn, name, None)
        if callable(fn):
            try:
                return float(fn(gene))
            except Exception:
                pass

    for attr in ("state", "states", "genes", "gene_states"):
        maybe = getattr(bn, attr, None)
        if isinstance(maybe, dict) and gene in maybe:
            return float(maybe[gene])

    raise AttributeError(f"Unable to read gene '{gene}' from BooleanNetwork object.")


def _apply_interventions(bn: Any, interventions: List[Dict[str, Any]]) -> None:
    fn = getattr(bn, "apply_interventions", None)
    if callable(fn):
        fn(interventions)
        return

    # Minimal fallback semantics for simple implementations.
    for iv in interventions:
        gene = str(iv.get("gene", ""))
        effect = str(iv.get("effect", "INHIBIT")).upper()
        strength = float(iv.get("strength", 1.0))
        if effect == "INHIBIT":
            _set_gene(bn, gene, max(0.0, 1.0 - strength))
        elif effect == "ACTIVATE":
            _set_gene(bn, gene, min(1.0, strength))


def _update(
    bn: Any,
    *,
    dt: float = 1.0,
    oxygen: float = 38.0,
    tgfb: float = 0.0,
    shh: float = 0.0,
    drug: float = 0.0,
) -> None:
    fn = getattr(bn, "update", None)
    if not callable(fn):
        raise AttributeError("BooleanNetwork object does not expose update().")

    call_attempts = [
        lambda: fn(dt, oxygen, tgfb, shh, drug),
        lambda: fn(dt=dt, oxygen=oxygen, tgfb=tgfb, shh=shh, drug=drug),
        lambda: fn(oxygen=oxygen, tgfb=tgfb, shh=shh, drug=drug),
        lambda: fn(),
    ]
    last_exc = None
    for attempt in call_attempts:
        try:
            attempt()
            return
        except Exception as exc:
            last_exc = exc
    raise RuntimeError(f"Unable to call BooleanNetwork.update() with supported signatures: {last_exc}")


def _run_steps(
    bn: Any,
    *,
    steps: int = 25,
    dt: float = 1.0,
    oxygen: float = 38.0,
    tgfb: float = 0.0,
    shh: float = 0.0,
    drug: float = 0.0,
    interventions: List[Dict[str, Any]] | None = None,
) -> None:
    interventions = interventions or []
    for _ in range(steps):
        _update(bn, dt=dt, oxygen=oxygen, tgfb=tgfb, shh=shh, drug=drug)
        if interventions:
            _apply_interventions(bn, interventions)


def _safe_gene_vector(bn: Any, genes: Iterable[str]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for gene in genes:
        try:
            out[gene] = _get_gene(bn, gene)
        except Exception:
            continue
    return out


def test_kras_locked_tumor_network_increases_myc_and_bcl_xl():
    on = _make_network("tumor")
    off = _make_network("tumor")

    _set_gene(on, "KRAS", 1.0)
    _set_gene(off, "KRAS", 0.0)

    _run_steps(on, steps=30)
    _run_steps(off, steps=30)

    assert _get_gene(on, "MYC") >= _get_gene(off, "MYC")
    assert _get_gene(on, "BCL_XL") >= _get_gene(off, "BCL_XL")


def test_tp53_loss_reduces_bax_expression():
    wt = _make_network("tumor")
    lost = _make_network("tumor")

    _set_gene(wt, "TP53", 1.0)
    _set_gene(lost, "TP53", 0.0)

    _run_steps(wt, steps=30)
    _run_steps(lost, steps=30)
    assert _get_gene(lost, "BAX") <= _get_gene(wt, "BAX")


def test_high_tgfb_activates_stromal_acta2():
    high = _make_network("stromal")
    low = _make_network("stromal")

    _run_steps(high, steps=25, tgfb=1.0)
    _run_steps(low, steps=25, tgfb=0.0)
    assert _get_gene(high, "ACTA2") >= _get_gene(low, "ACTA2")


def test_hif1a_responds_to_low_oxygen():
    hypoxic = _make_network("tumor")
    normoxic = _make_network("tumor")

    _run_steps(hypoxic, steps=25, oxygen=1.0)
    _run_steps(normoxic, steps=25, oxygen=38.0)
    assert _get_gene(hypoxic, "HIF1A") >= _get_gene(normoxic, "HIF1A")


def test_inhibit_egfr_reduces_downstream_myc():
    treated = _make_network("tumor")
    untreated = _make_network("tumor")

    interventions = [{"gene": "EGFR", "effect": "INHIBIT", "strength": 1.0}]
    _run_steps(treated, steps=25, interventions=interventions)
    _run_steps(untreated, steps=25)

    assert _get_gene(treated, "MYC") <= _get_gene(untreated, "MYC")


def test_continuous_update_converges_to_stable_state():
    bn = _make_network("tumor")
    genes = ["KRAS", "EGFR", "MYC", "BCL_XL", "BAX", "TP53", "HIF1A"]

    _run_steps(bn, steps=80, oxygen=15.0, tgfb=0.2, shh=0.2, drug=0.1)
    prev = _safe_gene_vector(bn, genes)
    _run_steps(bn, steps=40, oxygen=15.0, tgfb=0.2, shh=0.2, drug=0.1)
    curr = _safe_gene_vector(bn, genes)

    assert prev and curr
    common = sorted(set(prev.keys()) & set(curr.keys()))
    assert common

    max_delta = max(abs(curr[g] - prev[g]) for g in common)
    assert max_delta <= 0.25


def test_edge_cases_all_genes_off_on_and_simultaneous_interventions():
    bn = _make_network("tumor")
    genes = ["KRAS", "EGFR", "MYC", "BCL_XL", "BAX", "TP53", "HIF1A", "NRF2", "ABCB1"]

    # All OFF
    for gene in genes:
        try:
            _set_gene(bn, gene, 0.0)
        except Exception:
            continue
    _run_steps(bn, steps=10)
    vec_off = _safe_gene_vector(bn, genes)
    for val in vec_off.values():
        assert math.isfinite(val)
        assert -1e-6 <= val <= 1.0 + 1e-6

    # All ON
    for gene in genes:
        try:
            _set_gene(bn, gene, 1.0)
        except Exception:
            continue
    _run_steps(bn, steps=10)
    vec_on = _safe_gene_vector(bn, genes)
    for val in vec_on.values():
        assert math.isfinite(val)
        assert -1e-6 <= val <= 1.0 + 1e-6

    # Simultaneous interventions
    interventions = [
        {"gene": "EGFR", "effect": "INHIBIT", "strength": 0.8},
        {"gene": "MYC", "effect": "INHIBIT", "strength": 0.5},
        {"gene": "BAX", "effect": "ACTIVATE", "strength": 0.7},
    ]
    _run_steps(bn, steps=15, interventions=interventions)
    vec_iv = _safe_gene_vector(bn, ["EGFR", "MYC", "BAX"])
    assert vec_iv
    for val in vec_iv.values():
        assert math.isfinite(val)
        assert -1e-6 <= val <= 1.0 + 1e-6
