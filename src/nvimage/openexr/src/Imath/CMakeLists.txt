PROJECT (Imath)

SET(IMATH_SRCS	
#	ImathBox.cpp
	ImathBox.h
	ImathBoxAlgo.h
	ImathColor.h
	ImathColorAlgo.cpp
	ImathColorAlgo.h
	ImathEuler.h
	ImathExc.h
	ImathFrame.h
	ImathFrustum.h
	ImathFun.cpp
	ImathFun.h
	ImathGL.h
	ImathGLU.h
	ImathHalfLimits.h
	ImathInt64.h
	ImathInterval.h
	ImathLimits.h
	ImathLine.h
	ImathLineAlgo.h
	ImathMath.h
	ImathMatrix.h
	ImathMatrixAlgo.cpp
	ImathMatrixAlgo.h
	ImathPlane.h
	ImathPlatform.h
	ImathQuat.h
	ImathRandom.cpp
	ImathRandom.h
	ImathRoots.h
#	ImathShear.cpp
	ImathShear.h
	ImathSphere.h
	ImathVec.cpp
	ImathVec.h
	ImathVecAlgo.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(ILMBASE_SHARED)
	ADD_LIBRARY(Imath SHARED ${IMATH_SRCS})
ELSE(ILMBASE_SHARED)
	ADD_LIBRARY(Imath ${IMATH_SRCS})
ENDIF(ILMBASE_SHARED)

#TARGET_LINK_LIBRARIES(Imath ${LIBS})

INSTALL(TARGETS Imath
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib/static)

INSTALL(FILES 
	ImathBox.h
	ImathBoxAlgo.h
	ImathColor.h
	ImathColorAlgo.h
	ImathEuler.h
	ImathExc.h
	ImathFrame.h
	ImathFrustum.h
	ImathFun.cpp
	ImathFun.h
	ImathGL.h
	ImathGLU.h
	ImathHalfLimits.h
	ImathInt64.h
	ImathInterval.h
	ImathLimits.h
	ImathLine.h
	ImathLineAlgo.h
	ImathMath.h
	ImathMatrix.h
	ImathMatrixAlgo.h
	ImathPlane.h
	ImathPlatform.h
	ImathQuat.h
	ImathRandom.h
	ImathRoots.h
	ImathShear.h
	ImathSphere.h
	ImathVec.h
	ImathVecAlgo.h
	DESTINATION include/OpenEXR)
