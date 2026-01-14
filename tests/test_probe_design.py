import sys

from HCRProbeDesign import HCR
from HCRProbeDesign import probeDesign
from HCRProbeDesign import tiles


class _FakeHairpin:
    def __init__(self, tm=0.0, structure_found=False):
        self.tm = tm
        self.structure_found = structure_found


def test_fasta_input_creates_hcr_probes(monkeypatch, tmp_path, capsys):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(">target1\nACGTACGTAC\n")

    def fake_calc_hairpin(_seq):
        return _FakeHairpin()

    def fake_calc_tm(_seq):
        return 50.0

    monkeypatch.setattr(probeDesign.primer3, "calc_hairpin", fake_calc_hairpin)
    monkeypatch.setattr(probeDesign.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(tiles.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(probeDesign, "outputRunParams", lambda _args: None)

    argv = [
        "probeDesign",
        str(fasta_path),
        "-g",
        "--tileSize",
        "10",
        "--minGC",
        "0",
        "--maxGC",
        "100",
        "--minGibbs",
        "-1000",
        "--maxGibbs",
        "1000",
        "--targetGibbs",
        "0",
        "--maxRunLength",
        "999",
        "--maxProbes",
        "1",
        "--channel",
        "B1",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    probeDesign.main()

    out = capsys.readouterr().out.strip().splitlines()
    assert out[0].startswith("name\tprobe\tstart\tlength\tP1\tP2\tchannel")
    assert len(out) == 2

    row = out[1].split("\t")
    assert row[0].startswith("target1:")
    assert row[1] == "gtacgtacgt"
    assert row[4].startswith(HCR.initiators["B1"]["odd"])
    assert row[5].endswith(HCR.initiators["B1"]["even"])
    assert row[6] == "B1"
