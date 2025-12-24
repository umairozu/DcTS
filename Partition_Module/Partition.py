
from typing import Optional

class Partition:

    index: int
    payload: Optional[str] = None # rn a DNA string, may add encoded chunks

    """Supplementary text 1"""
    @staticmethod
    def partitions_for_label(label: str) -> int:
        n = len(label)
        return 9 + 3 * n

    def isEmpty(self) -> bool:
        return self.payload is None
