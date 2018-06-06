# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
PROJECT( PhaseSymmetryImageFilter)

SET(Required_Packages
ITK
)


ADD_EXECUTABLE( PhaseSymmetryImageFilter examplePS.cxx)

TARGET_LINK_LIBRARIES( PhaseSymmetryImageFilter ITKCommon ITKIO ITKBasicFilters ITKAlgorithms)