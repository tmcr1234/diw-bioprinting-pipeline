from aggregator import write_row, read_csv, COLUMNS


def _make_row(Vp, Flow, Pr_mean):
    row = {c: None for c in COLUMNS}
    row.update({
        "ink": "SF5.5", "bg": "Black", "shape": "S-Test",
        "Vp": Vp, "Flow": Flow,
        "Pr_consensus_mean": Pr_mean,
    })
    return row


def test_write_row_creates_file(tmp_path):
    csv = tmp_path / "out.csv"
    write_row(csv, _make_row(10, 1.57, 1.05))
    df = read_csv(csv)
    assert len(df) == 1
    assert df.iloc[0]["Pr_consensus_mean"] == 1.05


def test_write_row_overwrites_same_key(tmp_path):
    csv = tmp_path / "out.csv"
    write_row(csv, _make_row(10, 1.57, 1.05))
    write_row(csv, _make_row(10, 1.57, 1.20))
    df = read_csv(csv)
    assert len(df) == 1
    assert df.iloc[0]["Pr_consensus_mean"] == 1.20


def test_write_row_appends_distinct_keys(tmp_path):
    csv = tmp_path / "out.csv"
    write_row(csv, _make_row(10, 1.57, 1.05))
    write_row(csv, _make_row(15, 2.0, 0.95))
    df = read_csv(csv).sort_values(["Vp", "Flow"]).reset_index(drop=True)
    assert len(df) == 2
