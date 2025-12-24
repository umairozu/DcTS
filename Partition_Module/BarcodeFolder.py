from typing import Optional, List
from Partition import Partition


class BarcodeFolder:

    label: str
    slot: List[Partition]

    """
    
    [      (index 1, payload 1), (index2, payload 2), ...   ]
    ^^^         ^^^^
    slot      partition
    
    """

    def creating_slots(self):
        cap = Partition.partitions_for_label(self.label)
        self.slot = [Partition() for _ in range(cap)]

    def capacity(self) -> int:
        return len(self.slot)

    """first available space """
    def first_empty(self) -> Optional[Partition]:
        for p in self.slot:
            if p.isEmpty():
                return p
        return None

    def get_slot(self, index: int)-> Partition:
        if index < 1 or index > self.capacity():
            raise IndexError(f"desired partition {index} out of range, 1..{self.capacity()} for folder {self.label}")
        return self.slot[index]







