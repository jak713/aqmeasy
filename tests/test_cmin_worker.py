import os
import pytest
from pathlib import Path

# Lazy import to avoid heavy dependency loading during test collection
# (scipy, ASE, torch are pulled in by CMINWorker and cause SIGABRT on some systems)
def _import_cmin_worker():
    from aqmeasy.controllers.CMIN_worker import CMINWorker
    return CMINWorker


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")


def _create_sdf_with_elements(file_path, element_symbols):
    """Create a minimal SDF file with specified elements using V2000 format."""
    atom_count = len(element_symbols)
    sdf_content = f"""Simple molecule
  Test  01011212003D 0   0   0   0   0   0            999 V2000
M  CHG  0
{atom_count:>3}  0  0  0  0  0  0  0  0999 V2000
"""
    # Add atom block - each element at (0, 0, 0)
    for i, elem in enumerate(element_symbols):
        sdf_content += f"    0.0000    0.0000    0.0000 {elem:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n"
    
    sdf_content += "M  END\n$$$$\n"
    
    file_path.parent.mkdir(parents=True, exist_ok=True)
    file_path.write_text(sdf_content, encoding="utf-8")


def test_parse_sdf_summary_extracts_counts_and_energies(tmp_path):
        CMINWorker = _import_cmin_worker()
        sdf_path = tmp_path / "summary.sdf"
        sdf_path.write_text(
                """Mol 1
    AQME  01011212003D

    0  0  0  0  0  0            999 V2000
M  END
>  <Energy>
1.5

$$$$
Mol 2
    AQME  01011212003D

    0  0  0  0  0  0            999 V2000
M  END
>  <Energy>
0.5

$$$$
""",
                encoding="utf-8",
        )

        worker = CMINWorker([], parameters={})
        count, energies, structures = worker._parse_sdf_summary(str(sdf_path))

        assert count == 2
        assert energies == [1.5, 0.5]
        assert len(structures) == 2
        assert structures[0]["energy"] == 0.5
        assert structures[0]["rank"] == 1


def test_collect_results_reports_file_and_conformer_elimination_from_all_confs(tmp_path, monkeypatch):
    CMINWorker = _import_cmin_worker()
    input_a = tmp_path / "mol_a.sdf"
    input_b = tmp_path / "mol_b.sdf"
    _touch(input_a)
    _touch(input_b)

    cmin_dir = tmp_path / "CMIN"
    all_confs_dir = cmin_dir / "All_confs"
    filtered_a = cmin_dir / "mol_a_xtb.sdf"
    filtered_b = cmin_dir / "mol_b_xtb.sdf"
    all_a = all_confs_dir / "mol_a_xtb_all_confs.sdf"
    all_b = all_confs_dir / "mol_b_xtb_all_confs.sdf"
    _touch(filtered_a)
    _touch(filtered_b)
    _touch(all_a)
    _touch(all_b)

    fake_counts = {
        str(filtered_a): (2, [1.0, 1.5]),
        str(filtered_b): (0, []),
        str(all_a): (4, [1.0, 1.2, 1.5, 2.0]),
        str(all_b): (3, [0.5, 0.8, 1.1]),
        str(input_a): (4, []),
        str(input_b): (3, []),
    }

    def fake_parse(file_path):
        return fake_counts.get(str(file_path), (0, []))

    monkeypatch.setattr(CMINWorker, "_parse_sdf_summary", lambda self, p: (*fake_parse(p), []))

    worker = CMINWorker([str(input_a), str(input_b)], parameters={"program": "xtb"})
    results = worker._collect_results(output_dir="", run_dir=str(tmp_path))

    assert results["input_count"] == 2
    assert results["output_count"] == 2
    assert results["file_level_eliminated"] == 1
    assert results["conformer_level_eliminated"] == 5
    assert results["input_conformer_count"] == 7
    assert results["output_conformer_count"] == 2

    per_file = {row["input_name"]: row for row in results["per_file_data"]}
    assert per_file["mol_a.sdf"]["status"] == "partial"
    assert per_file["mol_a.sdf"]["eliminated_conformers"] == 2
    assert per_file["mol_b.sdf"]["status"] == "eliminated"
    assert per_file["mol_b.sdf"]["eliminated_conformers"] == 3


def test_collect_results_excludes_unmatched_outputs_from_counts(tmp_path, monkeypatch):
    CMINWorker = _import_cmin_worker()
    input_a = tmp_path / "starter.sdf"
    _touch(input_a)

    cmin_dir = tmp_path / "CMIN"
    all_confs_dir = cmin_dir / "All_confs"
    valid_filtered = cmin_dir / "starter_ani.sdf"
    valid_all = all_confs_dir / "starter_ani_all_confs.sdf"
    unrelated = cmin_dir / "random_results.sdf"
    _touch(valid_filtered)
    _touch(valid_all)
    _touch(unrelated)

    fake_counts = {
        str(valid_filtered): (1, [1.5]),
        str(valid_all): (2, [1.2, 1.5]),
        str(unrelated): (9, [0.1]),
        str(input_a): (2, []),
    }

    monkeypatch.setattr(
        CMINWorker,
        "_parse_sdf_summary",
        lambda self, p: (*fake_counts.get(str(p), (0, [])), []),
    )

    worker = CMINWorker([str(input_a)], parameters={"program": "ani"})
    results = worker._collect_results(output_dir="", run_dir=str(tmp_path))

    assert results["output_count"] == 1
    assert results["file_level_eliminated"] == 0
    assert results["conformer_level_eliminated"] == 1
    assert len(results["unmatched_output_files"]) == 1
    assert "did not match" in results["warnings"][0]


def test_collect_results_separates_failed_files_from_filtered_zero(tmp_path, monkeypatch):
    CMINWorker = _import_cmin_worker()
    input_a = tmp_path / "alpha.sdf"
    input_b = tmp_path / "beta.sdf"
    _touch(input_a)
    _touch(input_b)

    cmin_dir = tmp_path / "CMIN"
    filtered_a = cmin_dir / "alpha_xtb.sdf"
    all_a = cmin_dir / "alpha_xtb_all_confs.sdf"
    _touch(filtered_a)
    _touch(all_a)

    fake_counts = {
        str(filtered_a): (2, [0.1, 0.3]),
        str(all_a): (3, [0.1, 0.2, 0.3]),
        str(input_a): (3, []),
        str(input_b): (2, []),
    }

    monkeypatch.setattr(
        CMINWorker,
        "_parse_sdf_summary",
        lambda self, p: (*fake_counts.get(str(p), (0, [])), []),
    )

    worker = CMINWorker([str(input_a), str(input_b)], parameters={"program": "xtb"})
    results = worker._collect_results(output_dir="", run_dir=str(tmp_path))

    assert results["failed_file_count"] == 1
    assert results["file_level_eliminated"] == 0

    per_file = {row["input_name"]: row for row in results["per_file_data"]}
    assert per_file["beta.sdf"]["execution_status"] == "failed"
    assert per_file["beta.sdf"]["filter_outcome"] == "execution_failed"
    assert per_file["beta.sdf"]["status"] == "failed"
    assert per_file["beta.sdf"]["output_conformers"] is None


def test_prepare_aqme_paths_uses_common_root_relative_paths(tmp_path, monkeypatch):
    CMINWorker = _import_cmin_worker()
    file_a = tmp_path / "fromsource" / "CSEARCH" / "crest_xyz" / "mol_1.xyz"
    file_b = tmp_path / "aqme-test1_aqme" / "3030_rdkit.sdf"
    _touch(file_a)
    _touch(file_b)

    def fake_convert(self, xyz_path):
        sdf_path = os.path.splitext(xyz_path)[0] + ".sdf"
        with open(sdf_path, "w", encoding="utf-8") as handle:
            handle.write("")

    monkeypatch.setattr(CMINWorker, "_convert_xyz_to_sdf", fake_convert)

    worker = CMINWorker([], parameters={})
    aqme_files, run_dir, staging_dir = worker._prepare_aqme_paths([str(file_a), str(file_b)])

    assert run_dir == staging_dir
    assert os.path.isdir(staging_dir)
    assert aqme_files == [
        "mol_1.sdf",
        "3030_rdkit.sdf",
    ]
    assert os.path.isfile(os.path.join(staging_dir, "mol_1.sdf"))
    assert os.path.isfile(os.path.join(staging_dir, "3030_rdkit.sdf"))


def test_prepare_aqme_paths_single_file_returns_basename(tmp_path):
    CMINWorker = _import_cmin_worker()
    file_a = tmp_path / "one" / "molecule.sdf"
    _touch(file_a)

    worker = CMINWorker([], parameters={})
    aqme_files, run_dir, staging_dir = worker._prepare_aqme_paths([str(file_a)])

    assert run_dir == staging_dir
    assert aqme_files == ["molecule.sdf"]


def test_prepare_aqme_paths_rejects_duplicate_basenames(tmp_path):
    CMINWorker = _import_cmin_worker()
    file_a = tmp_path / "set_a" / "same.sdf"
    file_b = tmp_path / "set_b" / "same.sdf"
    _touch(file_a)
    _touch(file_b)

    worker = CMINWorker([], parameters={})

    with pytest.raises(ValueError, match="duplicate names"):
        worker._prepare_aqme_paths([str(file_a), str(file_b)])


def test_prepare_aqme_paths_rejects_duplicate_generated_sdf_names(tmp_path, monkeypatch):
    CMINWorker = _import_cmin_worker()
    file_a = tmp_path / "set_a" / "same.xyz"
    file_b = tmp_path / "set_b" / "same.sdf"
    _touch(file_a)
    _touch(file_b)

    monkeypatch.setattr(CMINWorker, "_convert_xyz_to_sdf", lambda self, path: None)

    worker = CMINWorker([], parameters={})

    with pytest.raises(ValueError, match="duplicate staged SDF names"):
        worker._prepare_aqme_paths([str(file_a), str(file_b)])


def test_validate_ani_input_elements_rejects_unsupported_elements(monkeypatch):
    CMINWorker = _import_cmin_worker()
    worker = CMINWorker([], parameters={})

    monkeypatch.setattr(
        CMINWorker,
        "_collect_input_atomic_numbers",
        lambda self, files: {1, 6, 8, 17},
    )

    with pytest.raises(ValueError, match="ANI1ccx") as exc:
        worker._validate_ani_input_elements("ani", ["dummy.sdf"], "ANI1ccx")

    message = str(exc.value)
    assert "Cl(Z=17)" in message
    assert "Supported elements for ANI1ccx" in message


def test_validate_ani_input_elements_accepts_supported_elements(monkeypatch):
    CMINWorker = _import_cmin_worker()
    worker = CMINWorker([], parameters={})

    monkeypatch.setattr(
        CMINWorker,
        "_collect_input_atomic_numbers",
        lambda self, files: {1, 6, 8, 17},
    )

    worker._validate_ani_input_elements("ani", ["dummy.sdf"], "ANI2x")


def test_validate_ani_input_elements_reports_when_no_ani_method_can_handle_input(monkeypatch):
    CMINWorker = _import_cmin_worker()
    worker = CMINWorker([], parameters={})

    monkeypatch.setattr(
        CMINWorker,
        "_collect_input_atomic_numbers",
        lambda self, files: {1, 15, 35, 46, 77},
    )

    with pytest.raises(ValueError) as exc:
        worker._validate_ani_input_elements("ani", ["dummy.sdf"], "ANI2x")

    message = str(exc.value)
    assert "No ANI method in AQME supports this element set" in message
    assert "Choose ANI2x" not in message


def test_validate_ani_input_elements_ignores_non_ani_program(monkeypatch):
    CMINWorker = _import_cmin_worker()
    worker = CMINWorker([], parameters={})

    monkeypatch.setattr(
        CMINWorker,
        "_collect_input_atomic_numbers",
        lambda self, files: {17},
    )

    worker._validate_ani_input_elements("xtb", ["dummy.sdf"], "ANI1ccx")
