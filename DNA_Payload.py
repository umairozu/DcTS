
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from OligoSequence import OligoSequence


TOTAL_SPACE = 4

"""Payload "STORED" in one partition"""
@dataclass
class DNA_Payload:
    # this could be folder's label or maybe chunks identifier
    #source_id: str

    # oligo's as list of [oligo sequence , Co_copies, encapsulated]
    oligos: List[Tuple[Optional[OligoSequence.sequence], int, bool]] = field(default_factory=list)


    #Remove this if not used
    """def __post_init__(self):
        if not self.oligos:
            self.oligos = [(None,0,False)] * TOTAL_SPACE"""

    def isEmpty(self):
        for item in self.oligos:
            if item[0] is None or item[0] == "":
                return True
        return False
    def get_copies(self):
        return self.oligos[0][1]



