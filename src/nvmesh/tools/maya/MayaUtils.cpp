// Copyright NVIDIA Corporation 2007 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmesh/MeshBuilder.h>
#include <nvmesh/halfedge/halfedgemesh.h>

// Required by maya!
#include <iostream>

#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MFnMesh.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MPointArray.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MItDependencyGraph.h>

#include "MayaUtils.h"
#include "MayaMeshBuilder.h"

// References:
//  http://florian.loitsch.com/gpExport/index-4.html
//  http://www.greggman.com/pages/mayastuff.htm#skin
//  http://www.robthebloke.org/research/index.htm
//
// tired of MFnSkinCluster::getWeights?? Try this:
//  http://www.highend3d.com/boards/index.php?showtopic=210482
//
// MTransformationMatrix
//  http://www.gamedev.net/community/forums/topic.asp?topic_id=253944


using namespace nv;





HalfEdge::Mesh * MayaUtils::buildHalfEdgeMesh(const MDagPath & dagPath)
{
	MayaMeshBuilderOptions options;
	MayaMeshBuilder builder(options);

	builder.addNode(dagPath);

	return builder.buildHalfEdgeMesh();
}


bool MayaUtils::getMeshPositions(const MDagPath & dagPath, Array<Vector3> * pointArray)
{
	nvDebugCheck(pointArray != NULL);

	MStatus status;
	MFnMesh meshFn(dagPath, &status);

	MItMeshPolygon polyIt(dagPath, MObject::kNullObj, &status);
	if (MS::kSuccess != status) return false;

	// Add positions.
	MPointArray positionArray;
	status = meshFn.getPoints(positionArray, MSpace::kObject);
	if (MS::kSuccess != status) return false;

	const uint positionCount = positionArray.length();
	pointArray->reserve(positionCount);
	
	for (uint i = 0; i < positionCount; i++)
	{
		MPoint point = positionArray[i];
		pointArray->append(Vector3(point.x, point.y, point.z));
	}

	return true;
}



MStatus MayaUtils::getSkinAndMesh(/*ref*/ MFnMesh & meshFn, /*out*/MFnSkinCluster & skinCluster)
{
	MStatus status;
	MObject meshObj = meshFn.object(&status);

	if (MS::kSuccess != status)
	{
		return status;
	}

	MItDependencyGraph dgIter(meshObj,
		MFn::kSkinClusterFilter, 
		MItDependencyGraph::kUpstream, 
		MItDependencyGraph::kBreadthFirst,
		MItDependencyGraph::kNodeLevel,
		&status);

	if (MS::kSuccess != status)
	{
		return status;
	}
	if(dgIter.isDone())
	{
		return MS::kNotFound;
	}

	return skinCluster.setObject(dgIter.thisNode());
}


// Implementation described here:
// http://www.greggman.com/pages/mayastuff.htm
MStatus MayaUtils::getSkinAndMeshGMAN(/*ref*/ MFnMesh & meshFn, /*out*/MFnSkinCluster & skinCluster)
{
	MStatus status;
	MFnDagNode dagNode(meshFn.dagPath(), &status);

	// the deformed mesh comes into the visible mesh
	// through its "inmesh" plug
	MPlug inMeshPlug = dagNode.findPlug("inMesh", &status);

	if (status != MS::kSuccess)
	{
		return status;
	}

	if (!inMeshPlug.isConnected())
	{
		return MS::kNotFound;
	}

	// walk the tree of stuff upstream from this plug
	MItDependencyGraph dgIt(inMeshPlug,
		MFn::kInvalid,
		MItDependencyGraph::kUpstream,
		MItDependencyGraph::kDepthFirst,
		MItDependencyGraph::kPlugLevel,
		&status);

	if (status != MS::kSuccess)
	{
		return status;
	}

	dgIt.disablePruningOnFilter();

	for (; !dgIt.isDone(); dgIt.next())
	{
		MObject thisNode = dgIt.thisNode();

		// go until we find a skinCluster
		if (thisNode.apiType() == MFn::kSkinClusterFilter)	// @@ Why does he filter explicitely?
		{
			//MFnSkinCluster skinCluster(thisNode);
			skinCluster.setObject(thisNode);

			// get the mesh coming into the skinCluster.  This
			// is the mesh before being deformed but after
			// being edited/tweaked/etc.

			MPlug inputPlug = skinCluster.findPlug("input", &status);
			if (MS::kSuccess == status)
			{
				MPlug childPlug = inputPlug.elementByLogicalIndex(0);
				MPlug geomPlug = childPlug.child(0);

				MObject dataObj1;
				geomPlug.getValue(dataObj1);

				// let use this mesh instead of the visible one
				meshFn.setObject(dataObj1);

				return MS::kSuccess;
			}
		}
	}

	return MS::kNotFound;
}




Skeleton * buildSkeleton(const MDagPath & dagPath)
{
	MStatus status;
	MFnMesh meshFn(dagPath, &status);
	//MObject meshObj = meshFn.object(&status);

	MFnSkinCluster skinCluster;
	status = MayaUtils::getSkinAndMesh(/*ref*/meshFn, /*out*/skinCluster);

    MDagPathArray influencePaths;
    int numInfluencePaths = skinCluster.influenceObjects(influencePaths, &status);


	Array<MString> usedTransforms(numInfluencePaths);
    for (int i = 0; i < numInfluencePaths; i++)
    {
		MString name = influencePaths[i].partialPathName();
		nvDebug("Influence %d : %s\n", i, name.asChar());
        usedTransforms.append(name);
    }

	/*MItGeometry gIter(meshObj, &status);
	for (; !gIter.isDone(); gIter.next())
	{

	}*/

	// @@ TODO!



	return NULL;
}



#if 0 

MStatus ObjTranslator::ManageSkinning(MObject meshobj)
{
    MFnDagNode curDagNode(meshobj);
    MStatus status;
    MFnMesh mesh(meshobj);

    /* EASIER:
    MItDependencyGraph dgIter(m_mesh.object(),
                              MFn::kSkinClusterFilter,
                              MItDependencyGraph::kUpstream,
                              MItDependencyGraph::kBreadthFirst,
                              MItDependencyGraph::kNodeLevel,
                              &g_status);

    if(dgIter.isDone())
    {        
        return false;
    }
    // We found a skin cluster!
    MFnSkinCluster skinCluster(dgIter.thisNode(), &g_status); assert(g_status);
    */
    MObject skinclusterobj = getConnectedObjectType(meshobj, "inMesh", MFn::kSkinClusterFilter, true, false, 4, &status);
    MFnSkinCluster skinCluster(skinclusterobj, &status);
    if(!status)
    return status;
    DPF(("====> Found SkinCluster...\n"));
    MDagPathArray influencePaths;
    int numInfluencePaths;
    numInfluencePaths = skinCluster.influenceObjects( influencePaths, &status );
    m_usedTransforms.clear(); // we don't need any other reference but Bones. So let's clear
    for(int i=0; i<numInfluencePaths; i++)
    {
        DPF(("Influence %d : %s\n", i, influencePaths[i].partialPathName().asChar()));
        m_usedTransforms.push_back(influencePaths[i].partialPathName());
    }
    //
    // Weights and indices
    //
    m_inattr_bonesoffsets = m_inattr.size();
    m_inattr.push_back(VertexAttribute(MESH_BONESOFFSETS,0));
    m_outattr.push_back(VertexAttribute(MESH_BONESOFFSETS,0));
    m_inattr_bonesweights = m_inattr.size();
    m_inattr.push_back(VertexAttribute(MESH_BONESWEIGHTS,0));
    m_outattr.push_back(VertexAttribute(MESH_BONESWEIGHTS,0));

    std::vector<std::map<int, float> > weightArrayMap;
    
    MItGeometry gIter(meshobj, &status);
    m_numbones = 0;
    DPF(("\n"));
    for (; !gIter.isDone(); gIter.next())
    {
        // Get the weights of the ith coordinate
        // If there aren't any, then the skin cluster is old and not used
        uint weightCount = 0;
        MDoubleArray weights;
        MObject component = gIter.component(&status);
        status = skinCluster.getWeights(m_currentInfo.top().m_dagPath, component, weights, weightCount);
        if (status)
        {
            std::map<int, float> weightmap;
            float sum=0;
            for(unsigned int i=0; i<weights.length(); i++)
            {
                if(weights[i] > 0.051) 
                {
                    sum += (float)weights[i];
                    weightmap[i] = (float)weights[i];
                    if((int)weightmap.size() > m_numbones) 
                        m_numbones = weightmap.size();
                }
                /*else
                {
                    DPF(("Skipping weight %d (%f)\n", i, weights[i]))
                }*/
            }
            weightArrayMap.push_back(weightmap);
            //DPF(("%f | ", sum));
        }
        else
        {
            weightArrayMap.push_back(std::map<int, float>());
        }
    }
    DPF(("\nSkinning : %d bones for %d vertices\n",m_numbones, weightArrayMap.size()));
    //
    // Now create the serie of offsets/weights made of N x m_numbones in in_attr
    //
    // Set the # of components.
    // Warning : these components may be >4... will have to turn this into many sets of 4 components
    // OR: later, using bytes in attribs would be cool...
    m_inattr[m_inattr_bonesoffsets].m_numComponents = m_numbones;
    m_inattr[m_inattr_bonesweights].m_numComponents = m_numbones;
    m_outattr[m_inattr_bonesoffsets].m_numComponents = m_numbones;
    m_outattr[m_inattr_bonesweights].m_numComponents = m_numbones;
    //DPF(("Weight sums\n"));
    for(unsigned int i=0; i<weightArrayMap.size(); i++)
    {
        std::map<int, float>::iterator it = weightArrayMap[i].begin();
        //DPF((">>>> Weights : "))
        int nb = 0;
        float sum = 0;
        for(int j=0; j<m_numbones; j++)
        {
            if(it == weightArrayMap[i].end())
            {
                m_inattr[m_inattr_bonesoffsets].m_floatVector.push_back(0.0f);
                m_inattr[m_inattr_bonesweights].m_floatVector.push_back(0.0f);
                //DPF(("0 "));
            } else {
                if(nb >= 4) DPF(("WARNING: %d exceeding 4 weights for the model\n", i));
                // Note : we may want ot use 'int' for the offset... TODO
                m_inattr[m_inattr_bonesoffsets].m_floatVector.push_back((float)it->first);
                m_inattr[m_inattr_bonesweights].m_floatVector.push_back(it->second);
                //DPF(("(%d:%f) ", it->first, (float)it->second))
                sum += it->second;
                nb++;
                ++it;
            }
        }
        //DPF((" - "));
    }
    DPF(("\n"));
    return status;
}

#endif
