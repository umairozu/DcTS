from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Tuple, List
import re
from Partition_Module.BarcodeFolder import BarcodeFolder

@dataclass
class TapeFS:

    """
    Acceptable/ used addressing format is absolute.
    format is <label>_<index>
    e,g: JK Li_5, JK Li_22
    """

    folders: Dict[str, BarcodeFolder] = field(default_factory= dict)

    def createFolder(self,label: str) -> BarcodeFolder:
        if label in self.folders:
            return self.folders[label]
        folder = BarcodeFolder(label=label, slot= [])
        folder.creating_slots()
        self.folders[label] = folder
        return folder


    def address_info(self, address: str) ->  int:
        address = re.match(r"^(.*)_(\d+)$",address)
        if address is None:
            raise ValueError(f"invalid address provided. format is <label>_<index>") #ValueError is for invalid DATA boii
        index = int(address.group(2)) # 'group()' is a part of the text that was captured by parentheses () in the pattern
        return index


    """payload could/will be changed to encoded chunks instead of DNA string later"""
    def deposit(self,label: str, payload: str):
        folder = self.createFolder(label)
        p = folder.first_empty()
        if p is None:
            raise RuntimeError(f"folder '{label}' is full (capacity = {folder.capacity()})")
            # what if I want to deposit at a specific folder/location instead of creating a new one?
            # answer: createFolder() method automatically does that. check that method!!
            # it first checks if a folder corresponding to that label exists, otherwise creates a new folder
            # if the specified folder exists but is full, Error is generated
        p.payload = payload
        print(f"deposited at {label}_{p.index}")
        return f"{label}_{p.index}"

    """
    - checking address(label) at first : valid or not valid
    - valid label returns label, index 
    - next, payload is checked, empty or not empty
    - if empty raise error
    """
    # returning payload as string
    def retrieve(self, label: str) -> str:
        index = self.address_info(label)
        label = label.rsplit("_",1)[-2]
        if label not in self.folders:
            raise KeyError(f"folder '{label}' does not exist")
        retrieval_folder = self.folders[label]
        payload = retrieval_folder.get_slot(index).payload
        if payload is None:
            raise KeyError(f"no payload at {label}")
        return payload

    def removal(self, label: str) -> None:
        index = self.address_info(label)
        if label not in self.folders:
            raise KeyError(f"folder '{label}' does not exist")
        removal_folder = self.folders[label]
        payload = removal_folder.get_slot(index)
        if payload is None:
            raise KeyError(f"payload at {label}_{index} is empty, nothing for removal")
        payload = None

    """for that particular folder, display all partitions (index 1, false) info from that slot"""
    def list_folder(self, label:str) -> List[Tuple[int,bool]]:
        folder = self.folders[label]
        if folder is None:
            raise KeyError(f"Unknown folder '{label}'")
        p = folder.slot
        return [(p.index, not p.isEmpty()) for p in folder.slot]







