from __future__ import annotations

import numpy as np
import pytest

from python.wrapper.output_parser import OutputParser


def test_parser_extracts_correct_cell_counts(mock_physicell_output):
    parser = OutputParser(mock_physicell_output)
    metrics = parser.parse_final_state()

    assert metrics.total_tumor_cells == 2
    assert metrics.total_stromal_cells == 2
    assert metrics.live_tumor_cells == 1
    assert metrics.activated_cafs == 2
    assert metrics.mesenchymal_tumor_cells == 1
    assert metrics.mean_tumor_drug_sensitivity == pytest.approx(0.45, rel=1e-6)
    assert metrics.tumor_extent == pytest.approx(100.0, rel=1e-6)


def test_stroma_barrier_score_computation(mock_physicell_output):
    parser = OutputParser(mock_physicell_output)

    tumor_positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [100.0, 0.0, 0.0],
        ]
    )
    coords = np.array(
        [
            [0.0, 0.0, 0.0],    # radius 50 from centroid: excluded
            [150.0, 0.0, 0.0],  # radius 100: included
            [250.0, 0.0, 0.0],  # radius 200: included
            [-100.0, 0.0, 0.0], # radius 150: included
        ]
    )
    ecm_values = np.array([0.1, 0.5, 0.6, 0.7])

    score = parser._compute_stroma_barrier_score(  # noqa: SLF001
        tumor_positions, {"coords": coords, "values": ecm_values}
    )
    assert score == pytest.approx((0.5 + 0.6 + 0.7) / 3.0, rel=1e-6)


def test_parse_timeseries(mock_physicell_output):
    parser = OutputParser(mock_physicell_output)
    df = parser.parse_timeseries()

    assert len(df) == 2
    assert list(df["time"]) == [0.0, 60.0]
    assert "tumor_count" in df.columns
    assert "mean_ecm" in df.columns
    assert "drug_penetration" in df.columns

    final = df.iloc[-1]
    assert int(final["tumor_count"]) == 2
    assert int(final["stroma_count"]) == 2
    assert int(final["live_tumor_cells"]) == 1
    assert float(final["stroma_barrier_score"]) == pytest.approx(0.6, rel=1e-6)
    assert float(final["drug_penetration"]) == pytest.approx(0.3, rel=1e-6)
