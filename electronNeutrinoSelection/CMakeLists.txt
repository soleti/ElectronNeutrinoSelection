#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB})

art_make(
          LIBRARY_NAME
                        ElectronFilter
          LIB_LIBRARIES
                        larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        art_Framework_Core
                        art_Framework_Services_Registry
                        art_Framework_Principal
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
#                        uboone_SpaceChargeServices_SpaceChargeServiceMicroBooNE_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
                        ${SIMULATIONBASE}
                   	larsim_MCCheater_BackTracker_service
			uboonecode_uboone_UBXSec_DataTypes
                        uboonecode_uboone_UBXSec_Algorithms
                        uboonecode_uboone_SelectionTool_OpT0FinderBase
                        uboonecode_uboone_SelectionTool_OpT0FinderAlgorithms
                        uboonecode_uboone_BasicTool_GeoAlgo
                        uboone_UBFlashFinder
                        larevt_CalibrationDBI_IOVData
                        larevt_CalibrationDBI_Providers
          MODULE_LIBRARIES
                        ElectronFilter
                        larpandora_LArPandoraInterface
			larsim_MCCheater_BackTracker_service
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        art_Framework_Core
                        art_Framework_Services_Registry
                        art_Framework_Principal
                        art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
                        ${SIMULATIONBASE}
			uboonecode_uboone_UBXSec_DataTypes
                        uboonecode_uboone_BasicTool_GeoAlgo
                        uboonecode_uboone_SelectionTool_OpT0FinderBase
                        uboonecode_uboone_SelectionTool_OpT0FinderAlgorithms
                        uboone_UBFlashFinder
                        uboonecode_uboone_UBXSec_Algorithms
                        EventWeight
                        larevt_CalibrationDBI_IOVData
                        larevt_CalibrationDBI_Providers
)

add_subdirectory(job)

install_headers()
install_source()
install_fhicl()
