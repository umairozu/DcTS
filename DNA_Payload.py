
from dataclasses import dataclass
from typing import List, Optional
from OligoSequence import OligoSequence


"""Payload "STORED" in one partition"""
@dataclass
class DNA_Payload:
    source_id: str # this could be folder's label or maybe chunks identifier
    encapsulated: bool
    oligos: List[OligoSequence]
    Co_copies: Optional[float] = None


