import sys
from pathlib import Path
import pandas as pd

SESSION = Path(__file__).parent / "fixtures" / "session"


def test_pipeline_runs_end_to_end(monkeypatch):
    import pipeline
    monkeypatch.setattr(pipeline, "PROJECT_ROOT", SESSION)
    monkeypatch.setattr(
        sys, "argv",
        ["pipeline.py", "--shape", "S-Test", "--ink", "SF5.5", "--bg", "Black",
         "--no-sam2", "--workers", "1"]
    )
    assert pipeline.main() == 0
    csv = SESSION / "Analises/Python/Printability/results/S-Test/SF5.5/Black/master_printability.csv"
    assert csv.exists()
    df = pd.read_csv(csv)
    assert len(df) == 2
    assert "Pr_consensus_mean" in df.columns
    assert "agreement_score" in df.columns


def test_pipeline_uses_cache_on_rerun(monkeypatch, capsys):
    import pipeline
    monkeypatch.setattr(pipeline, "PROJECT_ROOT", SESSION)
    monkeypatch.setattr(
        sys, "argv",
        ["pipeline.py", "--shape", "S-Test", "--ink", "SF5.5", "--bg", "Black",
         "--no-sam2", "--workers", "1"]
    )
    assert pipeline.main() == 0
    captured = capsys.readouterr()
    assert "SKIP" in captured.out
