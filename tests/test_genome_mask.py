import os

from HCRProbeDesign import genomeMask as gm


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
    monkeypatch.setattr(gm, "package_directory", "/pkg")
    monkeypatch.chdir(tmp_path)

    gm.genomemask(">read1\nACGT\n", handleName="test", index="indices/mm10/mm10")

    idx = captured["cmd"].index("-x") + 1
    assert captured["cmd"][idx] == os.path.join("/pkg", "indices/mm10/mm10")
