import os
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Tuple, List, Optional
import re

from natsort import natsorted

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

    """Supplementary text 1"""
    @staticmethod
    def partitions_for_label(label: str) -> int:
        n = len(label)
        return 9 + 3 * n

    @staticmethod
    def address_info(address: str) ->  int:
        address = re.match(r"^(.*)_(\d+)$",address)
        if address is None:
            raise ValueError(f"invalid address provided. format is <label>_<index>") #ValueError is for invalid DATA boii
        index = int(address.group(2)) # 'group()' is a part of the text that was captured by parentheses () in the pattern
        return index


    def createFolder(self,label: str):
        new_path = fr'{os.getcwd()}\{label}'
        num_partitions = self.partitions_for_label(label)
        if not os.path.exists(new_path):
            os.mkdir(new_path)
            for i in range(num_partitions):
                path_0 = f"Partition{i}.txt"
                path_1 = fr'{os.getcwd()}\{label}\{path_0}'
                open(path_1,'w')
        else:
            print(f"folder '{label}' already exists!")

    def deposit(self, label, DNA_Payload: DNA_Payload):
        index = self.address_info(label)
        label = label.rsplit("_", 1)[-2]
        folder = fr'{os.getcwd()}\{label}'

        if os.path.exists(folder):
            deposit_file = fr'{folder}\Partition{index}.txt'

            with open(deposit_file, "a+") as file:
                file.seek(0)
                num_lines = len(file.readlines())
                #print(num_lines)
                deposition_amount = DNA_Payload.__getcopies__()
                if deposition_amount and deposition_amount <= (10 - num_lines) :
                        for _ in range(deposition_amount):
                            file.writelines(str(DNA_Payload.oligos[0][0] +","+ DNA_Payload.oligos[0][2].__str__())+"\n") # if needed we can rewrite num_copies as well with data_payload & encapsulated
                            print(f"Deposition Successful at {label}_{index}")
                else:
                    print(f"Partition_{index} can not accommodate this deposition, choose a different location or update number of copies \n"
                          f"deposition amount = {deposition_amount} \n"
                          f"Available Space = {10 - num_lines}")


    def empty_partition(self, label):
        index = self.address_info(label)
        label = label.rsplit("_", 1)[-2]
        folder = fr'{os.getcwd()}\{label}'
        if os.path.exists(folder):
            file = fr'{folder}\Partition{index}.txt'
            os.truncate(file,0)
            print(f"Partition_{index} Emptied Successfully!")


    def retrieval(self, label, Adapter_L, Adapter_R):
        index = self.address_info(label)
        label = label.rsplit("_", 1)[-2]
        folder = fr'{os.getcwd()}\{label}'
        if os.path.exists(folder):
            file = fr'{folder}\Partition{index}.txt'

            #creating a separate .txt file for retrieval, if its empty, deleting it
            counter = 1
            retrieval_file = fr'{folder}\retrieval_from_Partition{index}({counter}).txt'

            while os.path.exists(retrieval_file):
                counter += 1
                retrieval_file = fr'{folder}\retrieval_from_Partition{index}({counter}).txt'

            with open(file) as read_file:

                open(retrieval_file, "a+")

                retrieval_count = 0
                print("Retrieved Oligos:")
                for line in read_file:
                    if line[:20] == Adapter_L and line[124:144] == Adapter_R:
                        with open(retrieval_file, "a+") as new_file:
                            new_file.writelines(line[:144] + "\n")
                        print(f"{line[:144]}")
                        retrieval_count += 1
                print(f"Retrieved count = {retrieval_count}")

                if retrieval_count == 0:
                    print(f"No Oligo/s found with the Adapters provided at Partition_{index}")
                    os.remove(retrieval_file)
        else:
            print(f"Specified Folder doesn't exist [{label}_{index}]")

    def removal(self, label, Adapter_L, Adapter_R):
        index = self.address_info(label)
        label = label.rsplit("_", 1)[-2]
        folder = fr'{os.getcwd()}\{label}'

        if os.path.exists(folder):
            file = fr'{folder}\Partition{index}.txt'
            with open(file) as f:
                lines = f.readlines()

            new_lines = []
            removal_count = 0
            for line in lines:
                if line[:20] == Adapter_L and line[124:144] == Adapter_R: # Avoid the specified sequences and then append the same file with remaining oligos
                    removal_count += 1
                else:
                    new_lines.append(line)
            if removal_count == 0:
                print("Nothing to remove!")

            with open(file, "w") as file:
                file.writelines(new_lines)

    @staticmethod
    def list_folder(label):
        folder = fr'{os.getcwd()}\{label}'
        if os.path.exists(folder):
            list_f = os.listdir(folder)
            for i, partition in enumerate(natsorted(list_f)):
                partition = fr'{folder}\{partition}'
                with open(partition) as file:
                    lines = file.readlines()
                    print(list_f[i], f"\nRemaining space in partition: {10- len(lines)}\n")
        else:
            print(f"folder '{label}' doesn't exist in the current directory [{os.getcwd()}] ")


    def edit_Oligo(self, label, Adapter_L, Adapter_R, new_payload):
        index = self.address_info(label)
        label = label.rsplit("_", 1)[-2]
        folder = fr'{os.getcwd()}\{label}'

        if os.path.exists(folder):
            file = fr'{folder}\Partition{index}.txt'
            edit_count = 0
            with open(file) as f:
                lines = f.readlines()

                new_lines = []
                for line in lines:
                    if line[:20] == Adapter_L and line[124:144] == Adapter_R:
                        print(f"previous Oligo payload: {line[40:104]}")
                        line = line[:40] + new_payload + line[104:]
                        print(f"new Oligo payload: {new_payload}\n")
                        edit_count += 1
                        new_lines.append(line)
                    else:
                        new_lines.append(line)
            with open(file, "w+") as file:
                file.writelines(new_lines)

            if edit_count == 0:
                print("Nothing to edit at the specified location!")
        else:
            print(f"folder '{label}' doesn't exist in the current directory [{os.getcwd()}] ")


    ################################################################

    """    
    folders: Dict[str, BarcodeFolder] = field(default_factory= dict)

    def createFolder(self,label: str) -> BarcodeFolder:
        if label in self.folders:
            return self.folders[label]
        folder = BarcodeFolder(label=label, slot= [])
        folder.creating_slots()
        self.folders[label] = folder
        return folder
    """

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





    """
    - checking address(label) at first : valid or not valid
    - valid label returns label, index 
    - next, payload is checked, empty or not empty
    - if empty raise error
    """
    # returning payload as DNA_Payload for the required DNA oligo and specified index
    """
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
        
        """




    """
    
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
                    

    """




    """for that particular folder, display all partitions (index 1, false) info from that slot"""
    """
    
    def list_folder(self, label:str) -> List[Tuple[int,bool]]:
        folder = self.folders[label]
        if folder is None:
            raise KeyError(f"Unknown folder '{label}'")
        return [(p.index, p.payload.isEmpty()) for p in folder.slot]
        
    """



    """
    - For insertion operation below
    - first creating a temporary '.txt' file and loading existing payload data from the partition into the text file
    - Opening created temporary file in notepad via subprocess tool
    - then user edit text
    - after edit, save and close file
    - Open the save file and read edited text from it
    - return edited/ updated text
    """


    """
    
    @staticmethod
    def edit_in_notepad(text):
        with tempfile.NamedTemporaryFile(suffix = ".txt", delete = False, mode = 'w') as tf:
            tf.write(text)
            temp_path = tf.name
        #1. Launch notepad, 2. open specified file inside it
        subprocess.run(['notepad.exe', temp_path])

        with open(temp_path, 'r') as f:
            updated_text = f.read()

        os.remove(temp_path)
        return updated_text

    """

    """
    Insertion operation on file data
    - adapter_l & adapter_R are used as oligo identifiers!!
    """
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


        """

                        ###################################################









