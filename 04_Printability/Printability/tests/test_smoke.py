def test_package_importable():
    import importlib
    # Modules will be added as tasks complete; this baseline checks the folder is a package.
    spec = importlib.util.find_spec("segmentation")
    assert spec is not None, "segmentation/ should be a package"
