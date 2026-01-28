import os
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Tuple, List, Optional
import re

from DNA_Payload import DNA_Payload
from OligoSequence import OligoSequence
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


    @staticmethod
    def address_info(address: str) ->  int:
        address = re.match(r"^(.*)_(\d+)$",address)
        if address is None:
            raise ValueError(f"invalid address provided. format is <label>_<index>") #ValueError is for invalid DATA boii
        index = int(address.group(2)) # 'group()' is a part of the text that was captured by parentheses () in the pattern
        return index


    """
    def deposit(self, label: str, sequence: OligoSequence, num_copies: int, encapsulated: bool):
        folder = self.createFolder(label)
        partition = folder.first_empty()
        if partition is None:
            raise RuntimeError(f"folder '{label}' is full (capacity = {folder.capacity()})")
            # what if I want to deposit at a specific folder/location instead of creating a new one?
            # answer: createFolder() method automatically does that. check that method!!
            # it first checks if a folder corresponding to that label exists, otherwise creates a new folder
            # if the specified folder exists but is full, Error is generated

        partition.payload.oligos.append((sequence,num_copies,encapsulated)) # rn , no limit on appending amount of oligos inside the partition, appending as many oligos as wished in a single partition
        print(f"deposited at {label}_{partition.index}")
        return f"{label}_{partition.index}"
    """

    def deposit(self, label, sequence, num_copies, encapsulated):
        index = self.address_info(label)
        label = label.rsplit("_", 1)[-2]
        folder = self.createFolder(label)
        partition = folder.slot[index]
        # woking on : if there is a space in payload[] insert there instead of appending at the end
        found_space = False
        for i, item in enumerate(partition.payload.oligos):
            if item[0] is None:
                partition.payload.oligos[i] = (sequence,num_copies,encapsulated)
                found_space = True
                print(f"1.deposited at {label}_{index}")
                break
        if not found_space:
            #partition.payload.oligos.append((sequence,num_copies,encapsulated))
            print(f"2. Can't find space at {label}_{index}")
        return f"{label}_{index}"


    """
    - checking address(label) at first : valid or not valid
    - valid label returns label, index 
    - next, payload is checked, empty or not empty
    - if empty raise error
    """
    # returning payload as DNA_Payload for the required DNA oligo and specified index
    def retrieve(self, label: str, adapter_L: str, adapter_R: str) -> OligoSequence:
        index = self.address_info(label)
        label = label.rsplit("_",1)[-2]
        if label not in self.folders:
            raise KeyError(f"folder '{label}' does not exist")
        retrieval_folder = self.folders[label]

        p = retrieval_folder.slot[index]
        sequences = [item[0] for item in p.payload.oligos]

        payload: Optional[OligoSequence] = None
        for seq in sequences:
            if seq is not None:
                if seq.adapter_L == adapter_L and seq.adapter_R == adapter_R:
                    payload = seq
                    print("RETRIEVED SUCCESSFULLY!")
                    print(f"Retrieved Oligo: {seq.sequence()}")
                    if payload is None or payload == "":
                        raise KeyError(f"no payload at {label}_{index}")
                    break
                else:
                    print(f"Adapters provided for retrieval didn't match with stored Oligo at location [{label}_{index}]")
        if payload is None:
            print("<--- PAYLOAD is NONE here, DON'T PRINT --->")
        return payload



    def removal(self, label: str, adapter_L: str, adapter_R: str) -> None:
        index = self.address_info(label)
        label = label.rsplit("_", 1)[-2]
        if label not in self.folders:
            raise KeyError(f"folder '{label}' does not exist")
        removal_folder = self.folders[label]

        # payload = removal_folder.get_slot(index)
        p = removal_folder.slot[index]
        if p is None:
            raise KeyError(f"payload at {label}_{index} is empty, nothing for removal")

        for i, item in enumerate(p.payload.oligos):
            seq = item[0]
            if seq is not None:
                if seq.adapter_L == adapter_L and seq.adapter_R == adapter_R:
                    p.payload.oligos[i] = (None,0,False)
                    print(f"PAYLOAD REMOVED SUCCESSFULLY! @ {label}_{index}")




    """for that particular folder, display all partitions (index 1, false) info from that slot"""
    def list_folder(self, label:str) -> List[Tuple[int,bool]]:
        folder = self.folders[label]
        if folder is None:
            raise KeyError(f"Unknown folder '{label}'")
        return [(p.index, p.payload.isEmpty()) for p in folder.slot]

    """
    - For insertion operation below
    - first creating a temporary '.txt' file and loading existing payload data from the partition into the text file
    - Opening created temporary file in notepad via subprocess tool
    - then user edit text
    - after edit, save and close file
    - Open the save file and read edited text from it
    - return edited/ updated text
    """
    @staticmethod
    def edit_in_notepad(text):
        with tempfile.NamedTemporaryFile(suffix = ".txt", delete = False, mode = 'w') as tf:
            tf.write(text)
            temp_path = tf.name
        """1. Launch notepad, 2. open specified file inside it"""
        subprocess.run(['notepad.exe', temp_path])

        with open(temp_path, 'r') as f:
            updated_text = f.read()

        os.remove(temp_path)
        return updated_text

    """Insertion operation on file data"""
    """
    - adapter_l & adapter_R are used as oligo identifiers!!
    """
    def insertion(self, label: str, adapter_L: str, adapter_R: str):
        partition_index = self.address_info(label)
        label = label.rsplit("_", 1)[-2]
        folder = self.folders[label]
        if folder is None:
            raise KeyError(f"Unknown folder '{label}'")
        p = folder.slot[partition_index]

        # As per paper, a payload can have multiple different oligos in it, find the correct oligo and get payload
        sequences = [item[0] for item in p.payload.oligos]
        for seq in sequences:
            if seq.adapter_L == adapter_L and seq.adapter_R == adapter_R:
                payload = seq.data_payload
                break

        updated_text = self.edit_in_notepad(payload)

        for seq in sequences:
            if seq is not None:
                if seq.adapter_L == adapter_L and seq.adapter_R == adapter_R:
                    seq.data_payload = updated_text
                    print(seq.data_payload)

                ###################################################









