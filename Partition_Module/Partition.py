from dataclasses import dataclass
from typing import Optional
from DNA_Payload import DNA_Payload

@dataclass
class Partition:

    index: int
    payload: Optional[DNA_Payload] = None # see DNA_payload & OligoSequence class for info

    """Supplementary text 1"""
    @staticmethod
    def partitions_for_label(label: str) -> int:
        n = len(label)
        return 9 + 3 * n

    def isEmpty(self) -> bool:
        return self.payload is None
