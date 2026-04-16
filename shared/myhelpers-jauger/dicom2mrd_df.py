# Modified 09/2023 for handling fMRI datasets saved as enhanced dicoms
# In these datasets, the slices are stored as frames, but should be separated
# before sending to the server (better replicating the streaming scenario while scanning)

import pydicom
import argparse
import ismrmrd
import numpy as np
import os
import ctypes
import re
import base64

# Defaults for input arguments
defaults = {
    'outGroup':       'dataset',
}

# Lookup table between DICOM and MRD image types
imtype_map = {'M': ismrmrd.IMTYPE_MAGNITUDE,
              'P': ismrmrd.IMTYPE_PHASE,
              'R': ismrmrd.IMTYPE_REAL,
              'I': ismrmrd.IMTYPE_IMAG}

# Lookup table between DICOM and Siemens flow directions
venc_dir_map = {'rl'  : 'FLOW_DIR_R_TO_L',
                'lr'  : 'FLOW_DIR_L_TO_R',
                'ap'  : 'FLOW_DIR_A_TO_P',
                'pa'  : 'FLOW_DIR_P_TO_A',
                'fh'  : 'FLOW_DIR_F_TO_H',
                'hf'  : 'FLOW_DIR_H_TO_F',
                'in'  : 'FLOW_DIR_TP_IN',
                'out' : 'FLOW_DIR_TP_OUT'}

def CreateMrdHeader(dset):
    """Create MRD XML header from a DICOM file"""

    mrdHead = ismrmrd.xsd.ismrmrdHeader()

    mrdHead.measurementInformation                             = ismrmrd.xsd.measurementInformationType()
    mrdHead.measurementInformation.measurementID               = dset.SeriesInstanceUID
    mrdHead.measurementInformation.patientPosition             = dset.PatientPosition
    mrdHead.measurementInformation.protocolName                = dset.SeriesDescription
    mrdHead.measurementInformation.frameOfReferenceUID         = dset.FrameOfReferenceUID

    mrdHead.acquisitionSystemInformation                       = ismrmrd.xsd.acquisitionSystemInformationType()
    mrdHead.acquisitionSystemInformation.systemVendor          = dset.Manufacturer
    mrdHead.acquisitionSystemInformation.systemModel           = dset.ManufacturerModelName
    mrdHead.acquisitionSystemInformation.systemFieldStrength_T = float(dset.MagneticFieldStrength)
    mrdHead.acquisitionSystemInformation.institutionName       = dset.InstitutionName

    mrdHead.userParameters                                     = ismrmrd.xsd.userParametersType()
    mrdHead.userParameters.userParameterLong                   = [ismrmrd.xsd.userParameterLongType(name="MultiBandFactor", value=args.sms_factor)] 
    try:
        mrdHead.acquisitionSystemInformation.stationName       = dset.StationName
    except:
        pass

    mrdHead.experimentalConditions                             = ismrmrd.xsd.experimentalConditionsType()
    mrdHead.experimentalConditions.H1resonanceFrequency_Hz     = int(dset.MagneticFieldStrength*4258e4)

    enc = ismrmrd.xsd.encodingType()
    enc.trajectory                                              = ismrmrd.xsd.trajectoryType('cartesian')
    encSpace                                                    = ismrmrd.xsd.encodingSpaceType()
    encSpace.matrixSize                                         = ismrmrd.xsd.matrixSizeType()
    encSpace.matrixSize.x                                       = dset.Columns
    encSpace.matrixSize.y                                       = dset.Rows
    encSpace.matrixSize.z                                       = 1
    encSpace.fieldOfView_mm                                     = ismrmrd.xsd.fieldOfViewMm()
    if dset.SOPClassUID.name == 'Enhanced MR Image Storage':
        encSpace.fieldOfView_mm.x                               =       dset.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing[0]*dset.Rows
        encSpace.fieldOfView_mm.y                               =       dset.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing[1]*dset.Columns
        encSpace.fieldOfView_mm.z                               = float(dset.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].SpacingBetweenSlices)
    else:
        encSpace.fieldOfView_mm.x                               =       dset.PixelSpacing[0]*dset.Rows
        encSpace.fieldOfView_mm.y                               =       dset.PixelSpacing[1]*dset.Columns
        encSpace.fieldOfView_mm.z                               = float(dset.SpacingBetweenSlices)
    enc.encodedSpace                                            = encSpace
    enc.reconSpace                                              = encSpace
    enc.encodingLimits                                          = ismrmrd.xsd.encodingLimitsType()
    enc.parallelImaging                                         = ismrmrd.xsd.parallelImagingType()

    enc.encodingLimits.slice = ismrmrd.xsd.limitType(minimum=0, maximum=args.max_slice_number, center=0)

    enc.parallelImaging.accelerationFactor                      = ismrmrd.xsd.accelerationFactorType()
    if dset.SOPClassUID.name == 'Enhanced MR Image Storage':
        # DNF: check for existence of parallel reduction factor fields
        if hasattr(dset.SharedFunctionalGroupsSequence[0].MRModifierSequence[0], 'ParallelReductionFactorInPlane'):
            enc.parallelImaging.accelerationFactor.kspace_encoding_step_1 = dset.SharedFunctionalGroupsSequence[0].MRModifierSequence[0].ParallelReductionFactorInPlane
        else:
            enc.parallelImaging.accelerationFactor.kspace_encoding_step_1 = 1
        if hasattr(dset.SharedFunctionalGroupsSequence[0].MRModifierSequence[0], 'ParallelReductionFactorOutOfPlane'):
            enc.parallelImaging.accelerationFactor.kspace_encoding_step_2 = dset.SharedFunctionalGroupsSequence[0].MRModifierSequence[0].ParallelReductionFactorOutOfPlane
        else:
            enc.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 
    else:
        enc.parallelImaging.accelerationFactor.kspace_encoding_step_1 = 1
        enc.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1

    mrdHead.encoding.append(enc)

    mrdHead.sequenceParameters                                  = ismrmrd.xsd.sequenceParametersType()

    return mrdHead

def GetDicomFiles(directory):
    """Get path to all DICOMs in a directory and its sub-directories"""
    for entry in os.scandir(directory):
        # DNF: handle dicoms saved wtih CRL script without file extension
        #if entry.is_file() and (entry.path.lower().endswith(".dcm") or entry.path.lower().endswith(".ima")):
        if entry.is_file():
            yield entry.path
        elif entry.is_dir():
            yield from GetDicomFiles(entry.path)


def main(args):
    dsetsAll = []
    for entryPath in GetDicomFiles(args.folder):
        dsetsAll.append(pydicom.dcmread(entryPath))

    if args.time:
        print("Sorting slices by acquisition time")

    # Group by series number
    uSeriesNum = np.unique([dset.SeriesNumber for dset in dsetsAll])
    print("Found %d unique series from %d files in folder %s" % (len(uSeriesNum), len(dsetsAll), args.folder))

    print("Creating MRD XML header from file %s" % dsetsAll[0].filename)
    mrdHead = CreateMrdHeader(dsetsAll[0])
    print(mrdHead.toXML())
    
    imgAll = [None]*len(uSeriesNum)

    for iSer in range(len(uSeriesNum)):
        dsets = [dset for dset in dsetsAll if dset.SeriesNumber == uSeriesNum[iSer]]

        imgAll[iSer] = [None]*len(dsets)

        # Sort images by instance number, as they may be read out of order
        def get_instance_number(item):
            return item.InstanceNumber
        dsets = sorted(dsets, key=get_instance_number)

        # Build a list of unique SliceLocation and TriggerTimes, as the MRD
        # slice and phase counters index into these
        # uSliceLoc = np.unique([dset.SliceLocation for dset in dsets])
        # if dsets[0].SliceLocation != uSliceLoc[0]:
        #     uSliceLoc = uSliceLoc[::-1]
        
        # DNF: slice locations stored per frame, assume consistent over all datasets
        uSliceLoc = np.empty(len(dsets[0].PerFrameFunctionalGroupsSequence))
        for iSlice in range(len(dsets[0].PerFrameFunctionalGroupsSequence)):
            uSliceLoc[iSlice] = dsets[0].PerFrameFunctionalGroupsSequence[iSlice].PlanePositionSequence[0].ImagePositionPatient[2]

        try:
            # This field may not exist for non-gated sequences
            uTrigTime = np.unique([dset.TriggerTime for dset in dsets])
            if dsets[0].TriggerTime != uTrigTime[0]:
                uTrigTime = uTrigTime[::-1]
        except:
            uTrigTime = np.zeros_like(uSliceLoc)

        print("Series %d has %d volumes with %d slices and %d phases" % (uSeriesNum[iSer], len(dsets), len(uSliceLoc), len(uTrigTime)))
        
        # DNF: !!! Hack to add repetition number for fetal_neural_segmentation_volume_reg_df script
        currentRep = 0
        
        for iImg in range(len(dsets)):
            tmpDset = dsets[iImg]

            ###### Get all slice acquisition times and return list of indices of sorted list
            ###### EXAMPLE: [0.5, 0.1, 1.2, 1.5, 0.5, 0.1, 1.2, 1.5] -> [1, 5, 2, 6, 3, 7, 4, 8] (SMS factor of 2) (original list sorted by slice location so the slices are correctly ordered for volume)
            uSliceTime = np.empty(len(tmpDset.PerFrameFunctionalGroupsSequence))
            for iSlice in range(len(tmpDset.PerFrameFunctionalGroupsSequence)):
                uSliceTime[iSlice] = float(tmpDset.PerFrameFunctionalGroupsSequence[iSlice].FrameContentSequence[0].FrameAcquisitionDateTime) 

            sortedSliceTime = [i for i, x in sorted(enumerate(uSliceTime), key=lambda x: x[1])]
            
            tmpDset_pixel_array = tmpDset.pixel_array
            # Remove pixel data from pydicom class
            del tmpDset['PixelData']
            
            imgAll[iSer][iImg] = [None]*len(uSliceLoc)

            # DNF loop over frames (slices) in dataset
            for SliceNum in range(len(dsets[iImg].PerFrameFunctionalGroupsSequence)):
                # Create new MRD image instance.
                # pixel_array data has shape [row col], i.e. [y x].
                # from_array() should be called with 'transpose=False' to avoid warnings, and when called
                # with this option, can take input as: [cha z y x], [z y x], or [y x]
                if args.time: #if sorting by time, then index into array of indices sorted by acquisition time
                    iSlice = sortedSliceTime[SliceNum] #iSlice is set to whatever index value is stored at SliceNum index, and iSlice is used to index into array of slices throughout the loop
                else: #otherwise continue linearly (i.e 0->31) 
                    iSlice = SliceNum
                tmpMrdImg = ismrmrd.Image.from_array(np.squeeze(tmpDset_pixel_array[iSlice,:,:]), transpose=False)
                tmpMeta   = ismrmrd.Meta()
    
                try:
                    tmpMrdImg.image_type                = imtype_map[tmpDset.ImageType[2]]
                except:
                    #print("Unsupported ImageType %s -- defaulting to IMTYPE_MAGNITUDE" % tmpDset.ImageType[2])
                    tmpMrdImg.image_type                = ismrmrd.IMTYPE_MAGNITUDE
                
                perFrameInfo = tmpDset.PerFrameFunctionalGroupsSequence[iSlice]

                tmpMrdImg.field_of_view            = (perFrameInfo.PixelMeasuresSequence[0].PixelSpacing[0]*tmpDset.Rows, perFrameInfo.PixelMeasuresSequence[0].PixelSpacing[1]*tmpDset.Columns, perFrameInfo.PixelMeasuresSequence[0].SpacingBetweenSlices)
                tmpMrdImg.position                 = tuple(np.stack(perFrameInfo.PlanePositionSequence[0].ImagePositionPatient))
                tmpMrdImg.read_dir                 = tuple(np.stack(perFrameInfo.PlaneOrientationSequence[0].ImageOrientationPatient[0:3]))
                tmpMrdImg.phase_dir                = tuple(np.stack(perFrameInfo.PlaneOrientationSequence[0].ImageOrientationPatient[3:7]))
                tmpMrdImg.slice_dir                = tuple(np.cross(np.stack(perFrameInfo.PlaneOrientationSequence[0].ImageOrientationPatient[0:3]), np.stack(perFrameInfo.PlaneOrientationSequence[0].ImageOrientationPatient[3:7])))
                acqDateTime                        = perFrameInfo.FrameContentSequence[0].FrameAcquisitionDateTime
                tmpMrdImg.acquisition_time_stamp   = round((int(acqDateTime[4:6])*3600 + int(acqDateTime[6:8])*60 + int(acqDateTime[8:10]) + float(acqDateTime[10:]))*1000/2.5)
                try:
                    tmpMrdImg.physiology_time_stamp[0] = round(int(tmpDset.TriggerTime/2.5))
                except:
                    pass

                try:
                    ImaAbsTablePosition = tmpDset.get_private_item(0x0019, 0x13, 'SIEMENS MR HEADER').value
                    tmpMrdImg.patient_table_position = (ctypes.c_float(ImaAbsTablePosition[0]), ctypes.c_float(ImaAbsTablePosition[1]), ctypes.c_float(ImaAbsTablePosition[2]))
                except:
                    pass
    
                tmpMrdImg.image_series_index     = uSeriesNum.tolist().index(tmpDset.SeriesNumber)
                tmpMrdImg.image_index            = tmpDset.get('InstanceNumber', 0)
                tmpMrdImg.slice                  = uSliceLoc.tolist().index(perFrameInfo.PlanePositionSequence[0].ImagePositionPatient[2])
                
                
                # DNF: !!! Hack to add repetition number for fetal_neural_segmentation_volume_reg_df script
                tmpMrdImg.repetition = currentRep
                if(SliceNum == (len(uSliceLoc)-1)):
                    currentRep = currentRep + 1
                #print("Current slice: %d, current rep: %d" % (tmpMrdImg.slice, tmpMrdImg.repetition))
                
                try:
                    tmpMrdImg.phase                  = uTrigTime.tolist().index(tmpDset.TriggerTime)
                except:
                    pass
    
                try:
                    res  = re.search(r'(?<=_v).*$',     tmpDset.SequenceName)
                    venc = re.search(r'^\d+',           res.group(0))
                    dir  = re.search(r'(?<=\d)[^\d]*$', res.group(0))
    
                    tmpMeta['FlowVelocity']   = float(venc.group(0))
                    tmpMeta['FlowDirDisplay'] = venc_dir_map[dir.group(0)]
                except:
                    pass
    
                tmpMeta['SequenceDescription'] = "slimmon"
                tmpMeta['EchoTime'] = perFrameInfo.MREchoSequence[0].EffectiveEchoTime
                tmpMeta['RepetitionTime'] = tmpDset.SharedFunctionalGroupsSequence[0].MRTimingAndRelatedParametersSequence[0].RepetitionTime

    
                # Store the complete base64, json-formatted DICOM header so that non-MRD fields can be
                # recapitulated when generating DICOMs from MRD images
                #tmpMeta['DicomJson'] = base64.b64encode(tmpDset.to_json().encode('utf-8')).decode('utf-8')

    
                tmpMrdImg.attribute_string = tmpMeta.serialize()
                imgAll[iSer][iImg][SliceNum] = tmpMrdImg 

    # Create an MRD file
    print("Creating MRD file %s with group %s" % (args.outFile, args.outGroup))
    mrdDset = ismrmrd.Dataset(args.outFile, args.outGroup)
    mrdDset._file.require_group(args.outGroup)

    # Write MRD Header
    mrdDset.write_xml_header(bytes(mrdHead.toXML(), 'utf-8'))

    # Write all images
    for iSer in range(len(imgAll)):
        for iImg in range(len(imgAll[iSer])):
            for iSli in range(len(imgAll[iSer][iImg])):
                mrdDset.append_image("images_%d" % imgAll[iSer][iImg][iSli].image_series_index, imgAll[iSer][iImg][iSli])

    mrdDset.close()

if __name__ == '__main__':
    """Basic conversion of a folder of DICOM files to MRD .h5 format"""

    parser = argparse.ArgumentParser(description='Convert DICOMs to MRD file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('folder',            help='Input folder of DICOMs')
    parser.add_argument('-o', '--outFile',  help='Output MRD file')
    parser.add_argument('-g', '--outGroup', help='Group name in output MRD file')
    parser.add_argument('-sms', '--sms_factor', help='SMS Factor of Dataset')
    parser.add_argument('-max', '--max_slice_number', help='Maximum slice number per volume')
    parser.add_argument('--time', action='store_true', help='Sort by acquisition time') ##add --time to command to sort by acquisition time

    parser.set_defaults(**defaults)

    args = parser.parse_args()

    if args.outFile is None:
        args.outFile = os.path.basename(args.folder) + '.h5'

    main(args)
