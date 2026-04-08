import os

from HCRProbeDesign import genomeMask as gm
from HCRProbeDesign import _datadir


def test_genomemask_uses_absolute_index(monkeypatch, tmp_path):
    captured = {}

    def fake_call(cmd):
        captured["cmd"] = cmd
        return 0

    monkeypatch.setattr(gm.subprocess, "call", fake_call)
    monkeypatch.chdir(tmp_path)

    gm.genomemask(">read1\nACGT\n", handleName="test", index="/abs/index")

    idx = captured["cmd"].index("-x") + 1
    assert captured["cmd"][idx] == "/abs/index"


def test_genomemask_uses_relative_index(monkeypatch, tmp_path):
    captured = {}

    def fake_call(cmd):
        captured["cmd"] = cmd
        return 0

    monkeypatch.setattr(gm.subprocess, "call", fake_call)
    data_dir = str(tmp_path / "datadir")
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", data_dir)
    monkeypatch.chdir(tmp_path)

    gm.genomemask(">read1\nACGT\n", handleName="test", index="indices/mm10/mm10")

    idx = captured["cmd"].index("-x") + 1
    assert captured["cmd"][idx] == os.path.join(data_dir, "indices/mm10/mm10")
