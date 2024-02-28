import sys
import numpy as np 
from pyspcm import *
from spcm_tools import *




'''
Various functions to help out with programming the Spectrum Instruments DAC card
Unless indicated otherwise, these functions have been tested out on our M4i.6621-x8 card
They're supposed to be compatible with other SI cards, but proceed with caution

File created by Justin on 2023/11/2
'''














def to_AWG(array):
    amp = np.max(np.abs(array)) # array maximum 
    arr_to_AWG = (array / amp * 32767).astype(np.int16) # convert to 16 bit array
    return int(amp*1000), arr_to_AWG 
    


def performance_checks(hCard):
    # read type, function and sn and check for D/A card
    lCardType = int32(0)
    spcm_dwGetParam_i32(hCard, SPC_PCITYP, byref(lCardType))
    lSerialNumber = int32(0)
    spcm_dwGetParam_i32(hCard, SPC_PCISERIALNO, byref(lSerialNumber))
    lFncType = int32(0)
    spcm_dwGetParam_i32(hCard, SPC_FNCTYPE, byref(lFncType))

    sCardName = szTypeToName(lCardType.value)
    if lFncType.value == SPCM_TYPE_AO:
        sys.stdout.write("Found: {0} sn {1:05d}\n".format(sCardName, lSerialNumber.value))
    else:
        sys.stdout.write("This code does not support Card: {0} sn {1:05d}\n".format(sCardName, lSerialNumber.value))
        spcm_vClose(hCard)
        exit(1)



def single_mode_init(hCard, loops):
    '''
    Data generation from on-board memory after one trigger event.
    ***Will only detect ONE single trigger event.*** 
    Repeats the complete programmed memory *loops* times.
    If loops = 0, it loops continuously after one single trigger event.
    '''
    spcm_dwSetParam_i32(hCard, SPC_CARDMODE, SPC_REP_STD_SINGLE)
    spcm_dwSetParam_i64(hCard, SPC_LOOPS, int64(loops))
    return

def single_restart_mode_init(hCard, loops):
    '''
    Data generation from on-board memory after a trigger event.
    ***Use this mode if you have multiple trigger events***
    If loops = 0, the programmed memory replays once on every 
    trigger event
    If loops != 0, the complete programmed memory will be replayed
    once per trigger events until *loops* number of replays/trigger events
    have occured.
    '''
    spcm_dwSetParam_i32(hCard, SPC_CARDMODE, SPC_REP_STD_SINGLERESTART)
    spcm_dwSetParam_i64(hCard, SPC_LOOPS, int64(loops))
    return

def multi_mode_init(hCard, loops):
    '''
    In "MULTI" mode the on-board memory is divided into equal size segments determined by *SegmentSize*,
    *ReplayNum* is the number of segments to be replayed (ReplayNum*SegmentSize = total memsize). 
    At each trigger event the segments will be replayed (ie: if ReplayNum = 3, then we have segments 1, 2, 3;
    the first trigger event will play segment 1, the second trigger event will play segment 2, and so on).
    If loops = 0 segments will loop continuously at each trigger event,
    if loops = 1 segments will stop playing after ReplayNum of segments has been reached.

    *** Keep in mind available memsize is determined by the number of output channels ***
        1 Channel: available memsize =  2GSample, max SegmentSize = memsize/2
        2 Channels: available memsize =  1GSample, max SegmentSize = memsize/2
    '''
    spcm_dwSetParam_i32(hCard, SPC_CARDMODE, SPC_REP_STD_MULTI)
    spcm_dwSetParam_i64(hCard, SPC_LOOPS, int64(loops))

    # '''
    # Filler for when we figure out how this works
    # '''
    # print("This function isn't written yet")
    return



def sequence_mode_init(hCard, MaxSegments):
    '''
    Divides on-board memory into selected number of *MaxSegments* and plays sequence of segments in
    order defined by user.
    If using two channels, segments must be written in the following form: CH1 Seg0, CH2 Seg0, CH1 Seg1,
    CH2 Seg1, CH1 Seg2, CH2 Seg2, etc.
    '''
    spcm_dwSetParam_i32(hCard, SPC_CARDMODE, SPC_REP_STD_SEQUENCE)
    lBytesPerSample = int32(0)
    spcm_dwGetParam_i32(hCard, SPC_MIINST_BYTESPERSAMPLE,  byref(lBytesPerSample))
    spcm_dwSetParam_i32 (hCard, SPC_SEQMODE_MAXSEGMENTS, MaxSegments)
    print("If using 2 channels upload sequence in following form:\n CH1 Seg0, CH2 Seg0, CH1 Seg1, CH2 Seg1, etc.")
    return

def double_mode_init(hCard, loops):
    '''
    Outputs the same waveform on both channels.
    Test to see if both channels can be enabled or only channel 1.
    '''
    spcm_dwSetParam_i32(hCard, SPC_CARDMODE, SPC_DOUBLEOUT0)
    spcm_dwSetParam_i64(hCard, SPC_LOOPS, int64(loops))
    return    



def configure_external_trigger_single(hcard, source, level, slope, level2 = 1):
 
    # Turn off the software trigger
    spcm_dwSetParam_i32(hcard, SPC_TRIG_ORMASK, SPC_TMASK_NONE)
    
    # Level: trigger level in V.     
    if level < -10 or level > 10:
        print("Trigger level out of range")
        return
    cLevel = int32(int(level*1000))

    if source.upper() == 'EXT0':
        Tmask, Mcommand, Lcommand = SPC_TMASK_EXT0, SPC_TRIG_EXT0_MODE, SPC_TRIG_EXT0_LEVEL0
    else:
        Tmask, Mcommand, Lcommand = SPC_TMASK_EXT1, SPC_TRIG_EXT1_MODE,  SPC_TRIG_EXT1_LEVEL0
    
    if slope.upper() == 'POS':
        SlopeCommand = SPC_TM_POS
    elif slope.upper() == 'NEG':
        SlopeCommand = SPC_TM_NEG
    elif slope.upper() == 'HIGH':
        SlopeCommand = SPC_TM_HIGH
    elif slope.upper() == 'LOW':
        SlopeCommand = SPC_TM_LOW
    elif slope.upper() == 'BOTH':
        SlopeCommand = SPC_TM_BOTH
    elif slope.upper() == 'ARMNEG':
    #     #Re-arms the trigger on the negative edge at user defined LEVEL0, and triggers at LEVEL1
        cLevel2 = int32(int(level2*1000))
        SlopeCommand =  (SPC_TM_NEG | SPC_TM_REARM)
        spcm_dwSetParam_i32(hcard, Lcommand, cLevel2)
        spcm_dwSetParam_i32(hcard, SPC_TRIG_EXT0_LEVEL1, cLevel)
        spcm_dwSetParam_i32(hcard, Mcommand, SlopeCommand)
        spcm_dwSetParam_i32(hcard, SPC_TRIG_ORMASK, Tmask)
        return
    elif slope.upper() == 'ARMPOS':
    #     #Re-arms the trigger on the negative edge at user defined LEVEL0, and triggers at LEVEL1
        cLevel2 = int32(int(level2*1000))
        SlopeCommand =  (SPC_TM_NEG | SPC_TM_REARM)
        spcm_dwSetParam_i32(hcard, Lcommand, cLevel)
        spcm_dwSetParam_i32(hcard, SPC_TRIG_EXT0_LEVEL1, cLevel2)
        spcm_dwSetParam_i32(hcard, Mcommand, SlopeCommand)
        spcm_dwSetParam_i32(hcard, SPC_TRIG_ORMASK, Tmask)
        return
    elif slope.upper() == 'WINENTER':
         SlopeCommand =  SPC_TM_WINENTER
         cLevel2 = int32(int(level2*1000))
         spcm_dwSetParam_i32(hcard, Lcommand, cLevel2)
         spcm_dwSetParam_i32(hcard, SPC_TRIG_EXT0_LEVEL1, cLevel)
         spcm_dwSetParam_i32(hcard, Mcommand, SlopeCommand)
         spcm_dwSetParam_i32(hcard, SPC_TRIG_ORMASK, Tmask)
         return
    else:
        print("Invalid slope")
        return
    
    spcm_dwSetParam_i32(hcard, Lcommand, cLevel)
    spcm_dwSetParam_i32(hcard, Mcommand, SlopeCommand)
    spcm_dwSetParam_i32(hcard, SPC_TRIG_ORMASK, Tmask)

    return



def configure_external_trigger_double(hcard, level, slope):
    # UNTESTED

    
    # Turn off the software trigger
    spcm_dwSetParam_i32(hcard, SPC_TRIG_ORMASK, SPC_TMASK_NONE)
    
    # Level: trigger level in V.     
    if level < -10 or level > 10:
        print("Trigger level out of range")
        return
    cLevel = int32(int(level*1000))

    if slope.upper() == 'POS':
        SlopeCommand = SPC_TM_POS
    elif slope.upper() == 'NEG':
        SlopeCommand == SPC_TM_NEG
    else:
        print("Invalid slope")
        return

    spcm_dwSetParam_i32(hcard, SPC_TRIG_EXT0_LEVEL0, cLevel)
    spcm_dwSetParam_i32(hcard, SPC_TRIG_EXT1_LEVEL0, cLevel)
    spcm_dwSetParam_i32(hcard, SPC_TRIG_EXT0_MODE, SlopeCommand)
    spcm_dwSetParam_i32(hcard, SPC_TRIG_EXT1_MODE, SlopeCommand)
    spcm_dwSetParam_i32(hcard, SPC_TRIG_ORMASK, SPC_TMASK_EXT1 | SPC_TMASK_EXT0)
    return


def set_amplitude(hcard, channel, amplitude):
    if amplitude < 80 or amplitude > 2500:
        print("Choose a valid amplitude")
        return
    if channel == 1:
        spcm_dwSetParam_i32(hcard, SPC_AMP0 + 0 * (SPC_AMP1 - SPC_AMP0), int32(amplitude))
    elif channel == 2:
        spcm_dwSetParam_i32(hcard, SPC_AMP0 + 1 * (SPC_AMP1 - SPC_AMP0), int32(amplitude))
    elif type(channel) == type('string') and channel.upper() == 'BOTH':
        spcm_dwSetParam_i32(hcard, SPC_AMP0 + 0 * (SPC_AMP1 - SPC_AMP0), int32(amplitude))
        spcm_dwSetParam_i32(hcard, SPC_AMP0 + 1 * (SPC_AMP1 - SPC_AMP0), int32(amplitude))
    else:
        print('Choose valid channel')
    return


def get_amplitude(hcard, channel):
    amplitude = int32(0)
    if channel == 1:
        spcm_dwGetParam_i32(hcard, SPC_AMP0 + 0 * (SPC_AMP1 - SPC_AMP0), byref(amplitude))
    elif channel == 2:
        spcm_dwGetParam_i32(hcard, SPC_AMP0 + 1 * (SPC_AMP1 - SPC_AMP0), byref(amplitude))
    elif type(channel) == type('string') and channel.upper() == 'BOTH':
        spcm_dwGetParam_i32(hcard, SPC_AMP0 + 0 * (SPC_AMP1 - SPC_AMP0), byref(amplitude))
        spcm_dwGetParam_i32(hcard, SPC_AMP0 + 1 * (SPC_AMP1 - SPC_AMP0), byref(amplitude))
    else:
        print('Choose valid channel')
    return amplitude.value

def set_ext_clock(hCard, frequency, sample_rate):
    # sets clock reference mode to external
    spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE, SPC_CM_EXTREFCLOCK)
    spcm_dwSetParam_i32 (hCard, SPC_REFERENCECLOCK, frequency)
    spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, sample_rate)



def SpcMCalcSignal(hCard, dwLenInSamples, eShape, dwLoops, frequency, amp, dwGainP = 1, dwPhase_degree = 0):
    llMinFS, llMaxFS, llValue = int64(0), int64(0), int64(0)
    i = uint32(0)
    lResolution = int32(0)
    dwLenInSamples = int(dwLenInSamples)
    dScale = 0.0
    pnData = []
    lSampRate = int64(0)
    spcm_dwGetParam_i64(hCard, SPC_SAMPLERATE, byref(lSampRate))
    lMaxADCValue = int32(0)
    spcm_dwGetParam_i32(hCard, SPC_MIINST_MAXADCVALUE, byref(lMaxADCValue))

    # examine the resolution, bytewidth and min/max values
    dwByteWidth = 2
    llMinFS = -lMaxADCValue.value-1
    llMaxFS = lMaxADCValue.value
    dScale = (float(llMaxFS) * dwGainP)

    # calculation of different signal shapes
    dSineXScale = 2.0 * np.pi / dwLenInSamples * dwLoops
    dwBlockLen = dwLenInSamples // dwLoops
    dwBlockHalf = dwBlockLen // 2
    dwPosInBlock = 0
    dSpan = float(llMaxFS - llMinFS)
    dPhase_rad = 2.0 * np.pi * dwPhase_degree / 360.0
    dwPhase_samples = int((float(dwBlockLen) * dwPhase_degree) / 360.0)
    # if (eShape == "eSine" or eShape =="eInvertedSine") and (frequency == None):
    #     dwLenInSamples = int(dwLenInSamples)
    # elif (eShape == "eSine" or eShape =="eInvertedSine") and (frequency != None):
    # dwLenInSamples = adjust_sample_length(lSampRate.value, dwLenInSamples, frequency)
    if frequency != None and (eShape == "eSine" or eShape == "eInvertedSine"):
        n = 32*np.arange(1,100000)
        nsamp = np.argmin(residual_phase(frequency, lSampRate.value, n))
        dwLenInSamples = n[nsamp]

    for i in range(dwLenInSamples):
        dwPosInBlock = i % dwBlockLen

        # calculation of value
        if eShape == "eDCZero":
            llValue = 0
        elif eShape == "eDCPlusFS":
            llValue = llMaxFS
        elif eShape == "eDCMinusFS":
            llValue = llMinFS
        elif eShape == "eSine" or eShape == "eInvertedSine":
            llValue = int(dScale * np.sin(2*np.pi/(lSampRate.value/frequency) * i + dPhase_rad))
        elif eShape == "eRectangle" or eShape == "eInvertedRectangle":
            if (dwPosInBlock + dwPhase_samples) % dwBlockLen < dwBlockHalf:
                llValue = llMinFS
            else:
                llValue = llMaxFS
        elif eShape == "eTriangle" or eShape == "eInvertedTriangle":
            if (dwPosInBlock - dwPhase_samples) % dwBlockLen < dwBlockHalf:
                llValue = int(llMinFS + ((dwPosInBlock - dwPhase_samples) % dwBlockLen) * dSpan / dwBlockHalf)
            else:
                llValue = int(llMaxFS - (((dwPosInBlock - dwPhase_samples) % dwBlockLen) - dwBlockHalf) * dSpan / dwBlockHalf)
        elif eShape == "eSawtooth" or eShape == "eInvertedSawtooth":
            llValue = int(llMinFS + ((dwPosInBlock - dwPhase_samples) % dwBlockLen) * dSpan / dwBlockLen)
        else:
            print("Unknown signal shape selected, choose a different eShape\n")
            return False

        # invert sign for inverted waveforms
        if eShape in ["eInvertedSine", "eInvertedTriangle", "eInvertedSawtooth", "eInvertedRectangle"]:
            llValue *= -1

        # write value to array
        if llValue < llMinFS:
            llValue = llMinFS
        elif llValue > llMaxFS:
            llValue = llMaxFS
        pnData.append(int(llValue))
    return pnData

def residual_phase(frequency, SampleRate, LenSamples):
    return (frequency*LenSamples/SampleRate)%1


