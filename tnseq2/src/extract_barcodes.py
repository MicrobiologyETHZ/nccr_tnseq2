from tnseq2.src.sequence import stream_fa, FastA

from typing import Tuple, List

def extract_barcode_host(r1: FastA, tp2: str, bc2tp2:int=13, bcLen:int=17, before:bool=True) -> Tuple[str, str]:
    '''
   Extract barcode and host sequence from read with tp2.
   Return (None, None) if the barcode sequence is not complete (17bp)
    :param r1:
    :param tp2:
    :return:

    -----|BARCODE|----------|TN end sequence (tp2)|---Host------
    -----|-bcLen-|--bc2tp2--|---------tp2---------|-------------
    -----|-17bp--|---13bp---|---------15bp--------|----?--------
    ---(-30)---(-13)-------(0)---------------------------------
    '''

    splits: List[str] = r1.sequence.split(tp2) # check tha tp2 in sequence?
    if before:
        bcStart = -(bcLen+bc2tp2)
        bcEnd = -bc2tp2
        bcSeq = splits[0]
        hostSeq = splits[1]
    else:
        bcStart = bc2tp2
        bcEnd = bc2tp2+bcLen
        bcSeq = splits[1]
        hostSeq = splits[0]

    if len(bcSeq) >= bc2tp2 + bcLen:
        barcode: str = bcSeq[bcStart:bcEnd]
        host_sequence: str = hostSeq
        return barcode, host_sequence
    else:
        return None, None
