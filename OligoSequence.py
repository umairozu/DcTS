
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
    payload: str
    seed: str
    cutting_site: str
    adapter_R: str

    def sequence(self) -> str:
        return self.adapter_L + self.rs_code + self.payload + self.seed + self.cutting_site + self.adapter_R

    def validation(self):
        assert (len(self.adapter_L)    == 20)
        assert (len(self.rs_code)      == 20)
        assert (len(self.payload)      == 64)
        assert (len(self.seed)         == 16)
        assert (len(self.cutting_site) == 4)
        assert (len(self.adapter_R)    == 20)

        assert (len(self.sequence())   == 144)


