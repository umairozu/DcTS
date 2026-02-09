import re
from dataclasses import dataclass

@dataclass
class OligoSequence:
    """
    (from paper fig.4a)
    length (nt) for attributes below:
    adapter_L    -> 20 nt
    rs_code      -> 20 nt
    payload      -> 64 nt
    seed         -> 16 nt
    cutting_site -> 4 nt
    adapter_R    -> 20 nt
    """
    adapter_L: str
    rs_code: str
    data_payload: str
    seed: str
    cutting_site: str
    adapter_R: str
    def __post_init__(self):
        if not (re.fullmatch(r'[ACGT]{20}', self.adapter_L) and re.fullmatch(r'[ACGT]{20}', self.rs_code)
                and re.fullmatch(r'[ACGT]{64}', self.data_payload) and re.fullmatch(r'[ACGT]{16}', self.seed)
                and re.fullmatch(r'[ACGT]{4}', self.cutting_site) and re.fullmatch(r'[ACGT]{20}', self.adapter_R)):
            raise RuntimeError("Invalid sequence, allowed bases are [A/C/G/T] for specific length only.")


    def sequence(self) -> str:
        return self.adapter_L + self.rs_code + self.data_payload + self.seed + self.cutting_site + self.adapter_R

    def validation(self):
        assert (len(self.adapter_L)    == 20)
        assert (len(self.rs_code)      == 20)
        assert (len(self.data_payload) == 64)
        assert (len(self.seed)         == 16)
        assert (len(self.cutting_site) == 4)
        assert (len(self.adapter_R)    == 20)

        assert (len(self.sequence())   == 144)


