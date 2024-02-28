import sys

import numpy as np 
from pyspcm import *
from spcm_tools import *
from helper_functions import *



class Spectrum_Instruments():
    def __init__(self):
        # open card
        self.hCard = spcm_hOpen(create_string_buffer(b'/dev/spcm0'))
        if not self.hCard:
            sys.stdout.write("no card found...\n")
            return
        
        performance_checks(self.hCard)

        # Default settings
        self.set_sampling_rate(MEGA(500)) # Sample rate to 500 MSa/S
        spcm_dwSetParam_i32(self.hCard, SPC_TRIG_TERM, int32(0)) # EXT0 input impedance -> 1 kOhm

        return
    
    def close(self):
        '''
        This MUST be run before spinning up another instance of Spectrum_Instruments
        Otherewise, the instrument tries to open the already-open card, and the
        only resolution I've found is to kill the kernel.
        '''
        spcm_vClose(self.hCard)
        print("card successfully closed")
        return

    def get_error(self):
        szErrorTextBuffer = create_string_buffer(ERRORTEXTLEN)
        spcm_dwGetErrorInfo_i32(self.hCard, None, None, szErrorTextBuffer)
        sys.stdout.write("{0}\n".format(szErrorTextBuffer.value))
        return


    def set_sampling_rate(self, rate):
        if rate < 50e6 or rate > 625e6:
            print("Invalid sampling rate")
            return
        spcm_dwSetParam_i64(self.hCard, SPC_SAMPLERATE, uint64(int(rate)))
        return


    def get_sampling_rate(self):
        lSampRate = int64(0)
        spcm_dwGetParam_i64(self.hCard, SPC_SAMPLERATE, byref(lSampRate))
        return lSampRate.value
        
    SampleRate = property(lambda self: self.get_sampling_rate, lambda self, x: self.set_sampling_rate(x))


    def set_output_mode(self, mode = 'SINGLE', loops = 0, MaxSegments = None):
        if mode.upper() == 'SINGLE':
            single_mode_init(self.hCard, loops)
        elif mode.upper() == 'SINGLERESTART':
            single_restart_mode_init(self.hCard, loops)
        elif mode.upper() == 'MULTI':
            multi_mode_init(self.hCard, loops)
        elif mode.upper() == "SEQUENCE":
            sequence_mode_init(self.hCard, MaxSegments)
        elif mode.upper() == "DOUBLEOUT":
            double_mode_init(self.hCard, loops)
        else:
            print("This mode isn't supported yet.")
        return

    def set_enable_ch(self, ch):
        '''
        This is important for memory allocation... I think?
        The name of the game here is to have as few channels as
        you need enabled (again... I think. It's a little unclear)
        '''
        print(ch)
        if ch == 1 or ch == 2:
            spcm_dwSetParam_i64(self.hCard, SPC_CHENABLE, uint64(ch))
        elif type(ch) == type('string') and ch.upper() == 'BOTH':
            spcm_dwSetParam_i64(self.hCard, SPC_CHENABLE, uint64(3))
        else:
            print("Provide valid channel")
        return

    def get_enable_ch(self):
        lCh = int64(0)
        spcm_dwGetParam_i64(self.hCard, SPC_CHENABLE, byref(lCh))
        if lCh.value == 1 or lCh.value == 2:
            return lCh.value
        elif lCh.value == 3:
            return 'both'
        else:
            print("???")
        return

    EnableChannel = property(lambda self: self.get_enable_ch(), lambda self, x: self.set_enable_ch(x))

    def get_enable_output(self, channel):
        lCh = int64(0)
        if channel == 1:
            spcm_dwGetParam_i64(self.hCard, SPC_ENABLEOUT0, byref(lCh))
        elif channel == 2:
            spcm_dwGetParam_i64(self.hCard, SPC_ENABLEOUT1, byref(lCh))
        else:
            print("invalid channel")
            return
        return bool(lCh.value)

    def set_enable_output(self, channel, onoff):
        if type(onoff) == type(True):
            cmmd = int32(int(onoff))
        else:
            print("You need to pass a boolean here")
        if channel == 1:
            spcm_dwSetParam_i64(self.hCard, SPC_ENABLEOUT0, cmmd)
        elif channel == 2:
            spcm_dwSetParam_i64(self.hCard, SPC_ENABLEOUT1, cmmd)
        else:
            print("invalid channel")
        return
    
    Out1 = property(lambda self: self.get_enable_output(1), lambda self, x: self.set_enable_output(1, x))
    Out2 = property(lambda self: self.get_enable_output(2), lambda self, x: self.set_enable_output(2, x))
        

    def set_amplitude(self, channel, amp):
        set_amplitude(self.hCard, channel, amp)
        return

    def get_amplitude(self, channel):
        return get_amplitude(self.hCard, channel)
    
    Amp1 = property(lambda self: self.get_amplitude(1), lambda self, amp: self.set_amplitude(1, amp))
    Amp2 = property(lambda self: self.get_amplitude(2), lambda self, amp: self.set_amplitude(2, amp))


    def configure_trigger(self, source = "SOFTWARE", level = 2.2, slope = "POS", level2 = 1):
        '''
        The card has a *lot* of trigger capabilities
        I'm going to keep it simple here, and only have a few options:
        SOFTWARE: immediately start playing waveforms after the you hit start. Default option.
        EXT0 or EXT1: trigger from one of the two external ports
        EXT0EXT1: EXT0 and EXT1 both active, either channel triggers (EXT0 OR EXT1)
        If using an external trigger, specify trigger level and slope.
        This should cover most of our use cases.
        '''
        if source.upper() == "SOFTWARE":
            # Default software trigger
            spcm_dwSetParam_i32(self.hCard, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE)
            spcm_dwSetParam_i32(self.hCard, SPC_TRIG_ANDMASK,     0)
            spcm_dwSetParam_i32(self.hCard, SPC_TRIG_CH_ORMASK0,  0)
            spcm_dwSetParam_i32(self.hCard, SPC_TRIG_CH_ORMASK1,  0)
            spcm_dwSetParam_i32(self.hCard, SPC_TRIG_CH_ANDMASK0, 0)
            spcm_dwSetParam_i32(self.hCard, SPC_TRIG_CH_ANDMASK1, 0)
            spcm_dwSetParam_i32(self.hCard, SPC_TRIGGEROUT,       0)
        elif source.upper() == "EXT0" or source.upper == "EXT1":
            configure_external_trigger_single(self.hCard, source, level, slope, level2)
        elif source.upper() == "EXT0EXT1":
            configure_external_trigger_double(self.hCard, level, slope, level2)
        
        # spcm_dwSetParam_i32(self.hCard, SPC_TRIG_ANDMASK,     0)
        # spcm_dwSetParam_i32(self.hCard, SPC_TRIG_CH_ORMASK0,  0)
        # spcm_dwSetParam_i32(self.hCard, SPC_TRIG_CH_ORMASK1,  0)
        # spcm_dwSetParam_i32(self.hCard, SPC_TRIG_CH_ANDMASK0, 0)
        # spcm_dwSetParam_i32(self.hCard, SPC_TRIG_CH_ANDMASK1, 0)
        # spcm_dwSetParam_i32(self.hCard, SPC_TRIGGEROUT,       0)
        return

    def load_waveform(self, waveform, channel, samplerate = None):
        '''
        Takes in a waveform, which is assumed to be a numpy array
        waveform is in volts into a 50 Ohm load. A high impedance load will see amp x 2
        loads that waveform into the AWG, and prep the enabled channel to play it
        Option to also set the sample rate here
        '''

        # Define the buffer
        self.wf_length = len(waveform)
        lSetChannels = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_CHCOUNT,     byref(lSetChannels))
        llMemSamples = int64(int(self.wf_length/lSetChannels.value))
        # dwError = spcm_dwSetParam_i64(self.hCard, SPC_MEMSIZE,     int64(int(llMemSamples/2)))
        dwError = spcm_dwSetParam_i64(self.hCard, SPC_MEMSIZE,     llMemSamples)
        if dwError != ERR_OK:
            self.get_error()
            return

        lBytesPerSample = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_MIINST_BYTESPERSAMPLE,  byref(lBytesPerSample))

        # self.BufferSize = llMemSamples.value * lBytesPerSample.value * lSetChannels.value
        self.BufferSize = llMemSamples.value * lBytesPerSample.value

        qwBufferSize = uint64(self.BufferSize)
        
        # Skip the continuous buffer if/else statement. If transfering data
        # somehow becomes the limiting factor, we can tinker with continuous memory.
        pvBuffer = c_void_p()
        pvBuffer = pvAllocMemPageAligned(qwBufferSize.value)

        amp, wf_16bit = to_AWG(waveform) # convert waveform to 16 bit array and get amplitude

        # move the data into a ctypes buffer, to be sent to the AWG.
        pnBuffer = cast(pvBuffer, ptr16)
        memmove(pnBuffer, wf_16bit.ctypes.data_as(ptr16), llMemSamples.value*sizeof(int16))
        # memmove(pnBuffer, wf_16bit.ctypes.data_as(ptr16), llMemSamples.value*sizeof(int16))

        del wf_16bit # paranorid memory management, probably not necessary

        # Set channel amplitude and sample rate
        set_amplitude(self.hCard, self.EnableChannel, amp)
        # set_amplitude(self.hCard, channel, amp)

        if type(samplerate) == type(12):
            self.samplerate = samplerate

        print("Starting the DMA transfer and waiting until data is in board memory\n")
        spcm_dwDefTransfer_i64(self.hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, int32(0), pvBuffer, uint64(0), qwBufferSize)
        spcm_dwSetParam_i32(self.hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)
        print("... data has been transferred to board memory\n")
        return
    

    def load_waveform_2_channels(self, Ch1waveform, Ch2waveform, samplerate = None):
        '''
        Takes in two waveform, which is assumed to be a numpy array
        waveform is in volts into a 50 Ohm load. A high impedance load will see amp x 2
        loads that waveform into the AWG, and prep the enabled channel to play it
        Option to also set the sample rate here
        '''
        # interleave the data into a single numpy array
        waveform = (np.column_stack((Ch1waveform,Ch2waveform))).flatten()
        # Define the buffer
        self.wf_length = len(waveform)
        lSetChannels = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_CHCOUNT,     byref(lSetChannels))
        llMemSamples = int64(int(self.wf_length/lSetChannels.value))
        # dwError = spcm_dwSetParam_i64(self.hCard, SPC_MEMSIZE,     int64(int(llMemSamples/2)))
        dwError = spcm_dwSetParam_i64(self.hCard, SPC_MEMSIZE,     llMemSamples)
        if dwError != ERR_OK:
            self.get_error()
            return

        lBytesPerSample = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_MIINST_BYTESPERSAMPLE,  byref(lBytesPerSample))

        # self.BufferSize = llMemSamples.value * lBytesPerSample.value * lSetChannels.value
        self.BufferSize = llMemSamples.value * lBytesPerSample.value

        qwBufferSize = uint64(self.BufferSize)
        
        # Skip the continuous buffer if/else statement. If transfering data
        # somehow becomes the limiting factor, we can tinker with continuous memory.
        pvBuffer = c_void_p()
        pvBuffer = pvAllocMemPageAligned(qwBufferSize.value)

        amp, wf_16bit = to_AWG(waveform) # convert waveform to 16 bit array and get amplitude

        # move the data into a ctypes buffer, to be sent to the AWG.
        pnBuffer = cast(pvBuffer, ptr16)
        memmove(pnBuffer, wf_16bit.ctypes.data_as(ptr16), llMemSamples.value*sizeof(int16))
        # memmove(pnBuffer, wf_16bit.ctypes.data_as(ptr16), llMemSamples.value*sizeof(int16))

        del wf_16bit # paranorid memory management, probably not necessary

        # Set channel amplitude and sample rate
        set_amplitude(self.hCard, self.EnableChannel, amp)
        # set_amplitude(self.hCard, channel, amp)

        if type(samplerate) == type(12):
            self.samplerate = samplerate

        print("Starting the DMA transfer and waiting until data is in board memory\n")
        spcm_dwDefTransfer_i64(self.hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, int32(0), pvBuffer, uint64(0), qwBufferSize)
        spcm_dwSetParam_i32(self.hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)
        print("... data has been transferred to board memory\n")
        return

    def load_waveform_multi(self, waveform, channel, ReplayNum = None, SegmentSize= None, samplerate = None):
        '''
        Since "MULTI" mode requires more restrictions on memsize this function
        should take care of necessary changes
        *** condisering defining output mode within load_waveform to minimize number of functions
        Takes in a waveform, which is assumed to be a numpy array
        waveform is in volts into a 50 Ohm load. A high impedance load will see amp x 2
        loads that waveform into the AWG, and prep the enabled channel to play it
        Option to also set the sample rate here
        '''
        self.wf_length = len(waveform)
        modetype = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_CARDMODE,     byref(modetype))
        lSetChannels = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_CHCOUNT,     byref(lSetChannels))

        if modetype.value == SPC_REP_STD_MULTI:
            if SegmentSize == None:
                SegmentSize = self.wf_length
            llMemSamples = int64(SegmentSize * ReplayNum)
            if lSetChannels.value == 1 and (int(SegmentSize) > int(GIGA(1)) or llMemSamples.value > int(GIGA(2))):
                print("Segment size or ReplayNum is too large ")
                return
            elif lSetChannels.value == 2 and (int(SegmentSize) > int(GIGA(2)/4) or llMemSamples.value > int(GIGA(1))):
                print("Segment size or ReplayNum is too large ")
                return
            waveform = np.tile(waveform, ReplayNum)
            spcm_dwSetParam_i64 (self.hCard, SPC_SEGMENTSIZE, int64(SegmentSize))
            spcm_dwSetParam_i64(self.hCard, SPC_MEMSIZE,     llMemSamples)
        if modetype.value == SPC_REP_STD_SINGLE:
            llMemSamples = int64(self.wf_length)
            spcm_dwSetParam_i64(self.hCard, SPC_MEMSIZE,     llMemSamples)
        # Define the buffer
        lBytesPerSample = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_MIINST_BYTESPERSAMPLE,  byref(lBytesPerSample))

        self.BufferSize = llMemSamples.value * lBytesPerSample.value * lSetChannels.value
        qwBufferSize = uint64(self.BufferSize)
        
        # Skip the continuous buffer if/else statement. If transfering data
        # somehow becomes the limiting factor, we can tinker with continuous memory.
        pvBuffer = c_void_p()
        pvBuffer = pvAllocMemPageAligned(qwBufferSize.value)

        amp, wf_16bit = to_AWG(waveform) # convert waveform to 16 bit array and get amplitude

        # move the data into a ctypes buffer, to be sent to the AWG.
        pnBuffer = cast(pvBuffer, ptr16)
        memmove(pnBuffer, wf_16bit.ctypes.data_as(ptr16), llMemSamples.value*sizeof(int16))
        del wf_16bit # paranorid memory management, probably not necessary

        # Set channel amplitude and sample rate

        set_amplitude(self.hCard, self.EnableChannel, amp)

        if type(samplerate) == type(12):
            self.samplerate = samplerate

        print("Starting the DMA transfer and waiting until data is in board memory\n")
        spcm_dwDefTransfer_i64(self.hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, int32(0), pvBuffer, uint64(0), qwBufferSize)
        spcm_dwSetParam_i32(self.hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)
        print("... data has been transferred to board memory\n")
        return
    
    def set_output(self, out):
        if type(out) != type(True):
            print("set_output takes a boolean.")
            return
        if out == True:
            # spcm_dwSetParam_i32(self.hCard, SPC_TIMEOUT, 10000)
            spcm_dwSetParam_i32(self.hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER) #| M2CMD_CARD_WAITREADY)
            # if dwError == ERR_TIMEOUT:
            #     spcm_dwSetParam_i32(self.hCard, SPC_M2CMD, M2CMD_CARD_STOP)
        if out == False:
            spcm_dwSetParam_i32(self.hCard, SPC_M2CMD,  M2CMD_CARD_STOP)
        return
    
    def use_ext_clock(self, freq = SPC_CM_EXTREFCLOCK, samplerate = SampleRate):
        set_ext_clock(self.hCard, freq, samplerate)
        #External clock source must be connected before using this function.
        #Setting freq to any other value will use the external clock as a reference for making
        # a gate.
        '''Would set output to true in order to read error type, consider checking external clock lock elsewhere:
        if (spcm_dwSetParam_i32 (self.hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER) == ERR_CLOCKNOTLOCKED):
             print("External clock not locked, check connection")
        '''
        return
    
    def ext_clock_output(self, output = False):
        if output == False:
            spcm_dwSetParam_i32 (self.hCard, SPC_CLOCKOUT, 0)
        elif output == True:
            spcm_dwSetParam_i32 (self.hCard, SPC_CLOCKOUT, 1)
        return

    
    def upload_sequence(self, Step, waveform, loops, channel, llCondition = SPCSEQ_ENDLOOPALWAYS , NextStep = None):
        if NextStep == None:
            NextStep = Step + 1
        # Condition can be SPCSEQ_ENDLOOPALWAYS, SPCSEQ_ENDLOOPONTRIG, SPCSEQ_END
        if Step == 0:
            spcm_dwSetParam_i32 (self.hCard, SPC_SEQMODE_STARTSTEP, 0)
        spcm_dwSetParam_i32 (self.hCard, SPC_SEQMODE_WRITESEGMENT, Step)
        wf_len = len(waveform)
        spcm_dwSetParam_i32 (self.hCard, SPC_SEQMODE_SEGMENTSIZE, int32(wf_len))
        self.load_waveform(waveform, channel)
        llValue = int64((llCondition<< 32) | (loops << 32) | (NextStep << 16) | (Step))
        spcm_dwSetParam_i64 (self.hCard, SPC_SEQMODE_STEPMEM0 + Step, llValue)
        return
    
    def multi_IO_line(self, line, mode = SPCM_XMODE_DISABLE):
        '''
        Available modes: SPCM_XMODE_ASYNCIN, SPCM_XMODE_ASYNCOUT, SPCM_XMODE_DIGOUT, SPCM_XMODE_TRIGOUT,
        SPCM_XMODE_RUNSTATE, SPCM_XMODE_ARMSTATE, SPCM_XMODE_REFCLKOUT, SPCM_XMODE_CONTOUTMARK, 
        SPCM_XMODE_SYSCLKOUT
        Mode will only be updated
        with the next call to either the M2CMD_CARD_START or M2CMD_CARD_WRITESETUP register.
        '''
        if line == "X0":
            spcm_dwSetParam_i32 (self.hCard, SPCM_X0_MODE, mode)
        elif line == "X1":
            spcm_dwSetParam_i32 (self.hCard, SPCM_X1_MODE, mode)
        elif line == "X2":
            spcm_dwSetParam_i32 (self.hCard, SPCM_X2_MODE, mode)
        else:
            print("This 'line' is not valid.")
        return

    # def load_simple_waveform(self, channel, amp, dwLenInSamples, eShape, dwLoops, frequency = None, dwPhase_degree = 0, samplerate = None):
    #     waveform = SpcMCalcSignal(self.hCard, dwLenInSamples, eShape, dwLoops, frequency, dwGainP = 1, dwPhase_degree = dwPhase_degree)
    #     # Define the buffer
    #     self.wf_length = len(waveform)
    #     llMemSamples = int64(self.wf_length)
    #     spcm_dwSetParam_i64(self.hCard, SPC_MEMSIZE,     llMemSamples)

    #     lSetChannels = int32(0)
    #     spcm_dwGetParam_i32(self.hCard, SPC_CHCOUNT,     byref(lSetChannels))
    #     lBytesPerSample = int32(0)
    #     spcm_dwGetParam_i32(self.hCard, SPC_MIINST_BYTESPERSAMPLE,  byref(lBytesPerSample))

    #     self.BufferSize = llMemSamples.value * lBytesPerSample.value * lSetChannels.value
    #     qwBufferSize = uint64(self.BufferSize)
        
    #     # Skip the continuous buffer if/else statement. If transfering data
    #     # somehow becomes the limiting factor, we can tinker with continuous memory.
    #     pvBuffer = c_void_p()
    #     pvBuffer = pvAllocMemPageAligned(qwBufferSize.value)

    #     _, wf_16bit = to_AWG(waveform)# convert waveform to 16 bit array

    #     # move the data into a ctypes buffer, to be sent to the AWG.
    #     pnBuffer = cast(pvBuffer, ptr16)
    #     memmove(pnBuffer, wf_16bit.ctypes.data_as(ptr16), llMemSamples.value*sizeof(int16))
    #     del wf_16bit # paranorid memory management, probably not necessary

    #     amp = int(amp*1000) #Convert amplitude from V to mV

    #     # Set channel amplitude and sample rate
    #     set_amplitude(self.hCard, channel, amp)
    #     if type(samplerate) == type(12):
    #         self.samplerate = samplerate

    #     print("Starting the DMA transfer and waiting until data is in board memory\n")
    #     spcm_dwDefTransfer_i64(self.hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, int32(0), pvBuffer, uint64(0), qwBufferSize)
    #     spcm_dwSetParam_i32(self.hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)
    #     print("... data has been transferred to board memory\n")
    #     return
    
    def load_simple_waveform(self, channel, amp, dwLenInSamples, eShape, dwLoops, frequency = None, dwPhase_degree = 0, samplerate = None):
        waveform = SpcMCalcSignal(self.hCard, dwLenInSamples, eShape, dwLoops, frequency, amp)
        # Define the buffer
        self.wf_length = len(waveform)
        llMemSamples = int64(self.wf_length)
        spcm_dwSetParam_i64(self.hCard, SPC_MEMSIZE,     llMemSamples)

        lSetChannels = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_CHCOUNT,     byref(lSetChannels))
        lBytesPerSample = int32(0)
        spcm_dwGetParam_i32(self.hCard, SPC_MIINST_BYTESPERSAMPLE,  byref(lBytesPerSample))

        self.BufferSize = llMemSamples.value * lBytesPerSample.value * lSetChannels.value
        qwBufferSize = uint64(self.BufferSize)
        
        # Skip the continuous buffer if/else statement. If transfering data
        # somehow becomes the limiting factor, we can tinker with continuous memory.
        pvBuffer = c_void_p()
        pvBuffer = pvAllocMemPageAligned(qwBufferSize.value)

        wf_16bit = (np.array(waveform)).astype(np.int16) # convert waveform to 16 bit array
        # move the data into a ctypes buffer, to be sent to the AWG.
        pnBuffer = cast(pvBuffer, ptr16)
        memmove(pnBuffer, wf_16bit.ctypes.data_as(ptr16), llMemSamples.value*sizeof(int16))
        del wf_16bit # paranorid memory management, probably not necessary

        amp = int(amp*1000) #Convert amplitude from V to mV
        # print (amp)

        # Set channel amplitude and sample rate
        set_amplitude(self.hCard, channel, amp)
        if type(samplerate) == type(12):
            self.samplerate = samplerate

        # print("Starting the DMA transfer and waiting until data is in board memory\n")
        spcm_dwDefTransfer_i64(self.hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, int32(0), pvBuffer, uint64(0), qwBufferSize)
        dwError = spcm_dwSetParam_i32(self.hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)
        if dwError != ERR_OK:
            self.get_error()
            self.close()
        # print("... data has been transferred to board memory\n")
        return 




        
    
    
