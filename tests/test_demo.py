import json

from aav9_sma.demo import run_demo


def test_demo_is_self_contained(tmp_path) -> None:
    summary = run_demo(tmp_path)

    assert summary["mode"] == "synthetic_smoke_test"
    assert summary["rows"] == 64
    assert summary["passed_packaging_gate"] > 0
    for path in summary["outputs"].values():
        assert (tmp_path / path.split("/")[-1]).exists()
    assert json.loads((tmp_path / "demo_summary.json").read_text()) == summary
