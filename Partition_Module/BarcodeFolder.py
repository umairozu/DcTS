import dataclasses
from dataclasses import dataclass
from typing import Optional, List
from Partition_Module.Partition import Partition


@dataclass
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
        self.slot = [Partition(index= i) for i in range(cap)]

    def capacity(self) -> int:
        return len(self.slot)

    """first available space """
    def first_empty(self) -> Optional[Partition]:
        for p in self.slot:
            if p.isEmpty():
                return p
        return None

    def get_slot(self, index: int)-> Partition:
        if index < 0 or index > self.capacity():
            raise IndexError(f"desired partition {index} out of range, 0..{self.capacity() - 1} for folder {self.label}")
        return self.slot[index]







