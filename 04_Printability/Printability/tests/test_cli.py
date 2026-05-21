from cli import parse_args


def test_required_args_parsed():
    ns = parse_args(["--shape", "S-Test", "--ink", "SF5.5", "--bg", "Black"])
    assert ns.shape == "S-Test"
    assert ns.ink == "SF5.5"
    assert ns.bg == "Black"
    assert ns.use_sam2 is True


def test_no_sam2_flag():
    ns = parse_args(["--shape", "S-Test", "--ink", "SF5.5", "--bg", "Black", "--no-sam2"])
    assert ns.use_sam2 is False


def test_limit_flag():
    ns = parse_args(["--shape", "S-Test", "--ink", "SF5.5", "--bg", "Black", "--limit", "3"])
    assert ns.limit == 3
