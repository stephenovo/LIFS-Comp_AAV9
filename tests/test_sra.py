import xml.etree.ElementTree as ET

from aav9_sma.data.sra import _package_to_row, summarize_sra_manifest


def test_parses_sra_package_and_summarizes() -> None:
    package = ET.fromstring(
        """
        <EXPERIMENT_PACKAGE>
          <EXPERIMENT accession="SRX1" alias="hammerhead_Brain_a1_r1">
            <TITLE>Brain</TITLE>
          </EXPERIMENT>
          <SAMPLE accession="SRS1"><SAMPLE_ATTRIBUTES>
            <SAMPLE_ATTRIBUTE><TAG>library</TAG><VALUE>Hammerhead</VALUE></SAMPLE_ATTRIBUTE>
            <SAMPLE_ATTRIBUTE><TAG>type</TAG><VALUE>Biodistribution</VALUE></SAMPLE_ATTRIBUTE>
            <SAMPLE_ATTRIBUTE><TAG>sample</TAG><VALUE>Brain</VALUE></SAMPLE_ATTRIBUTE>
          </SAMPLE_ATTRIBUTES></SAMPLE>
          <RUN_SET><RUN accession="SRR1" size="10" total_spots="2" /></RUN_SET>
        </EXPERIMENT_PACKAGE>
        """
    )

    row = _package_to_row(package)
    summary = summarize_sra_manifest([row])

    assert row["run_accession"] == "SRR1"
    assert row["library"] == "Hammerhead"
    assert summary["run_count"] == 1
    assert summary["total_bytes"] == 10
