from dataclasses import dataclass, field
from typing import Optional, List
from DNA_Payload import DNA_Payload

@dataclass
class Partition:

    index: int
    # Ensuring every Partition get separate payload memory not shared as in
    # for example:
    # payload: Optional[DNA_Payload] = DNA_Payload()
    payload : DNA_Payload = field(default_factory = DNA_Payload)

    """Supplementary text 1"""
    @staticmethod
    def partitions_for_label(label: str) -> int:
        n = len(label)
        return 9 + 3 * n

    def isEmpty(self) -> bool:
        for item in self.payload.oligos:
            if item[0] is None:
                return True
        return False
