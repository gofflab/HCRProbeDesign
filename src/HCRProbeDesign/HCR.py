
# Split-initiator sequences from
# Sequences retrieved from https://dev.biologists.org/content/develop/suppl/2018/06/21/145.12.dev165753.DC1/DEV165753supp.pdf
# Sequences below include 2nt spacers on 3' end of odd and 5' end of even to allow for bending to present initiator.
initiators = {
    "B1": {
        "odd": "gAggAgggCAgCAAACggAA",
        "even": "TAgAAgAgTCTTCCTTTACg"
    },
    "B2": {
        "odd": "CCTCgTAAATCCTCATCAAA",
        "even": "AAATCATCCAgTAAACCgCC"
    },
    "B3": {
        "odd": "gTCCCTgCCTCTATATCTTT",
        "even": "TTCCACTCAACTTTAACCCg"
    },
    "B4": {
        "odd": "CCTCAACCTACCTCCAACAA",
        "even": "ATTCTCACCATATTCgCTTC"
    },
    "B5": {
        "odd": "CTCACTCCCAATCTCTATAA",
        "even": "AACTACCCTACAAATCCAAT"
    },
}

def addInitiator(tile,initiator="B1"):
    pass
