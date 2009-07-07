// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#include <nvmesh/TriMesh.h>
#include <nvmesh/QuadTriMesh.h>

#include <nvmesh/kdtree/MeshKDTree.h>
#include <nvmesh/kdtree/KDTreeBuilder.h>

using namespace nv;



KDTree * nv::buildKDTree(const TriMesh * mesh)
{
	KDTreeBuilder builder;

	// Allocate indices.
	const uint faceCount = mesh->faceCount();

	Array<uint> indexArray;
	indexArray.resize(faceCount * 3);

	// Copy face indices to temporary array.
	for(uint f = 0; f < faceCount; f++)
	{
		const TriMesh::Face & face = mesh->faceAt(f);

		indexArray[3 * f + 0] = face.v[0];
		indexArray[3 * f + 1] = face.v[1];
		indexArray[3 * f + 2] = face.v[2];
	}


	// Allocate vertices.
	const uint vertexCount = mesh->vertexCount();

	Array<Vector3> vertexArray;
	vertexArray.resize(vertexCount);

	// Copy vertex position to temporary array.
	for(uint v = 0; v < vertexCount; v++)
	{
		const TriMesh::Vertex & vertex = mesh->vertexAt(v);

		vertexArray[v] = vertex.pos;
	}

	
	// Create the tree.
	return builder.createTree(indexArray, vertexArray);
}


KDTree * nv::buildKDTree(const QuadTriMesh * mesh)
{
	KDTreeBuilder builder;

	// Allocate indices.
	const uint faceCount = mesh->faceCount();

	Array<uint> indexArray;
	indexArray.reserve(faceCount * 3);

	// Copy face indices to temporary array.
	for(uint f = 0; f < faceCount; f++)
	{
		const QuadTriMesh::Face & face = mesh->faceAt(f);

		indexArray.append(face.v[0]);
		indexArray.append(face.v[1]);
		indexArray.append(face.v[2]);

		if (face.isQuadFace())
		{
			indexArray.append(face.v[0]);
			indexArray.append(face.v[2]);
			indexArray.append(face.v[3]);
		}
	}


	// Allocate vertices.
	const uint vertexCount = mesh->vertexCount();

	Array<Vector3> vertexArray;
	vertexArray.resize(vertexCount);

	// Copy vertex position to temporary array.
	for(uint v = 0; v < vertexCount; v++)
	{
		const TriMesh::Vertex & vertex = mesh->vertexAt(v);

		vertexArray[v] = vertex.pos;
	}

	
	// Create the tree.
	return builder.createTree(indexArray, vertexArray);
}




#if 0

/** A set of rays to test the performance of the KD-Tree. */
struct RaySet {

	/** Implement Nick Chapman test to compare performance.
	 * http://homepages.paradise.net.nz/nickamy/benchmark.html
	 */
	void GenerateChapmanTest(int num) {

		ray_array.Resize(num);

		foreach(r, ray_array) {
			PiKDTree::Ray & ray = ray_array[r];

			// Generate origin.
			RandomUnitVector(&ray.origin);
		//	ray.origin *= 0.2;
		//	ray.origin += Vec3(-0.016840, 0.110154, -0.001537);

			// Generate direction.				
			Vec3 end;
			RandomUnitVector(&end);
		//	end *= 0.2;
		//	end += Vec3(-0.016840, 0.110154, -0.001537);

			ray.dir.Sub(end, ray.origin);
			ray.maxt = ray.dir.NormalizeSlow();
			ray.idir.Set(1.0f / ray.dir.x, 1.0f / ray.dir.y, 1.0f / ray.dir.z);
		}	
	}


	/** Generate rays inside the box with a random direction. */
	void GenerateCastanoTest(uint num) {

		ray_array.Resize(num);

		foreach(r, ray_array) {
			PiKDTree::Ray & ray = ray_array[r];

			// Generate origin.
			RandomUnitVector(&ray.origin);
			//ray.origin = Vec3::Origin;

			// Generate direction.
			ray.dir.Sub(Vec3::Origin, ray.origin);
			ray.maxt = ray.dir.NormalizeSlow();
			ray.idir.Set(1.0f / ray.dir.x, 1.0f / ray.dir.y, 1.0f / ray.dir.z);

			// Generate direction.						
			//RandomUnitVector(&ray.dir);
			//ray.maxt = 1;
		}	
	}

	/** Generate montecarlo test. */
	void GenerateMontecarloTest(uint num) {

		ray_array.Resize(num);

		for(uint i = 0; i < num; i += 128) {

			// Generate origin.
			Vec3 origin;
			RandomUnitVector(&origin);

			for(uint e = 0; e < 128; e++) {

				if( i + e == num ) return;

				PiKDTree::Ray & ray = ray_array[i+e];

				ray.origin = origin;

				ray.dir.Sub(Vec3::Origin, ray.origin);
				ray.maxt = ray.dir.NormalizeSlow() * 0.5f;
				ray.idir.Set(1.0f / ray.dir.x, 1.0f / ray.dir.y, 1.0f / ray.dir.z);
			}
		}
	}




	/** Test the performance of the KD-Tree. */
	void TestKDTree(const PiKDTree * tree) {

		PiKDTree::Hit hit;
		intersection_num = 0;
		node_test_num = 0;
		face_test_num = 0;
		leaf_test_num = 0;
		uint64 begin = piGetClockCounter();

		foreach(r, ray_array) {
			if( tree->TestRay(ray_array[r], &hit) ) {
				intersection_num++;
			}
			node_test_num += tree->node_test_num;
			leaf_test_num += tree->leaf_test_num;
			face_test_num += tree->face_test_num;
		}

		elapsed = piGetClockCounter() - begin;
	}


	/** Test the performance of the KD-Tree for montecarlo evaluation. */
	void TestMontecarlo(const PiKDTree * tree) {

		PiKDTree::Hit hit;
		intersection_num = 0;
		node_test_num = 0;
		face_test_num = 0;
		leaf_test_num = 0;
		uint64 begin = piGetClockCounter();

		const uint num = ray_array.Size();
		for(uint i = 0; i < num; i += 128) {

			PiKDTree::Ray & ray = ray_array[i];

			PiKDTree::Cache cache;
			tree->CreateSphericalCache(ray, &cache);

			for(uint e = 0; e < 128; e++) {

				if( i + e == num ) return;

				if( tree->TestRay(cache, ray, &hit) ) {
					intersection_num++;
				}
				node_test_num += tree->node_test_num;
				leaf_test_num += tree->leaf_test_num;
				face_test_num += tree->face_test_num;
			}
		}

		elapsed = piGetClockCounter() - begin;
	}


	/** Test the performance of the KD-Tree with shadow rays. */
	void TestShadowRays(const PiKDTree * tree) {
		intersection_num = 0;
		node_test_num = 0;
		face_test_num = 0;
		leaf_test_num = 0;
		uint64 begin = piGetClockCounter();

		foreach(r, ray_array) {
			if( tree->TestRay(ray_array[r]) ) {
				intersection_num++;
			}
			node_test_num += tree->node_test_num;
			leaf_test_num += tree->leaf_test_num;
			face_test_num += tree->face_test_num;
		}

		elapsed = piGetClockCounter() - begin;
	}


	/** Test the performance of the KD-Tree. */
	void TestBruteForce(const PiKDTree * tree) {

		PiKDTree::Hit hit;
		intersection_num = 0;
		node_test_num = 0;
		face_test_num = 0;
		leaf_test_num = 0;
		uint64 begin = piGetClockCounter();

		foreach(r, ray_array) {
			if( tree->TestAllFaces(ray_array[r], &hit) ) {
				intersection_num++;
			}
			face_test_num += tree->face_test_num;
		}

		elapsed = piGetClockCounter() - begin;
	}


	/** Get a random unit vector. */
	void RandomUnitVector(Vec3 * v) {
		float length;
		do {
			v->x = 2 * rnd.GetReal() - 1;
			v->y = 2 * rnd.GetReal() - 1;
			v->z = 2 * rnd.GetReal() - 1;
			length = v->Length();
		} while( length > 1.0f || length < PI_EPSILON );

		*v /= length;
	}


	/** Print test results. */
	void PrintResult() {

		const uint ray_num = ray_array.Size();
		float percentage = 100 * float(intersection_num) / ray_num;

		piDebug("  %d rays traced in %.6fM clocks.\n", ray_num, double(elapsed) / (1024*1024) );
		piDebug("  %.3fK clocks per ray.\n", double(elapsed) / (1024 * ray_num) );
		piDebug("  %d rays intersected (%.2f%%).\n", intersection_num, percentage);
		piDebug("  %.2f nodes per ray.\n", double(node_test_num) / ray_num);
		piDebug("  %.2f leafs per ray.\n", double(leaf_test_num) / ray_num);
		piDebug("  %.2f faces per ray.\n", double(face_test_num) / ray_num);
	}

	Rand48 rnd;
	PiArray<PiKDTree::Ray> ray_array;

	uint64 elapsed;
	uint intersection_num;
	uint node_test_num;
	uint face_test_num;
	uint leaf_test_num;
};



PI_DECLARE_TEST(MeshKDTreeTest) {

#if !_DEBUG
	// Try to increase the priority for better tunning.
	piSetRealtimePriority();
#endif

	// Allocation.
	PiMeshPtr mesh( new PiMesh() );

	// Importing.
	//if( mesh->Load("meshes/standford/bunny.ply") == NULL ) {
	if( mesh->Load("meshes/hoppe/cathead.m") == NULL ) {
		return Failed;
	}

	// Transform.
	PiMeshTransform::Ptr mesh_transform( new PiMeshTransform(mesh.GetPtr()) );
	mesh_transform->FitBox(Box(-1, -1, -1, 1, 1, 1));
	//mesh_transform->FlipAxis(0, 2, 1);

	Box bounds;
	PiMeshBoundsPtr mesh_bounds( new PiMeshBounds(mesh.GetPtr()) );
	mesh_bounds->GetBoundingBox(&bounds);
//	bounds.mins.Print();
//	bounds.maxs.Print();

	Vec3 center;
	bounds.GetCenter(&center);
//	center.Print();


	// Weld mesh vertices.
//	PiMeshVertexWeldPtr mesh_vertex_weld( mesh->GetPlugin<PiMeshVertexWeld>() );
//	mesh_vertex_weld->WeldVertices();

	// AABB Tree.
	//PiMeshAABBTreePtr mesh_aabb_tree( mesh->GetPlugin<PiMeshAABBTree>() );
	//mesh_aabb_tree->BuildTree();
	
	PiMeshKDTreePtr mesh_kd_tree( new PiMeshKDTree(mesh.GetPtr()) );
	mesh_kd_tree->BuildTree();

	// Implement Nick Chapman test to compare performance.
	// http://homepages.paradise.net.nz/nickamy/benchmark.html

	RaySet ray_set;
//	ray_set.GenerateChapmanTest(512 * 1024);
//	ray_set.TestKDTree(mesh_kd_tree->GetTree());
//	ray_set.PrintResult();

//	ray_set.TestShadowRays(mesh_kd_tree->GetTree());
//	ray_set.PrintResult();

	ray_set.GenerateMontecarloTest(512 * 1024 * 2);
	ray_set.TestKDTree(mesh_kd_tree->GetTree());
	ray_set.PrintResult();

	ray_set.TestMontecarlo(mesh_kd_tree->GetTree());
	ray_set.PrintResult();

//	ray_set.TestBruteForce(mesh_kd_tree->GetTree());
//	ray_set.PrintResult();
	
//	ray_set.GenerateCastanoTest(512 * 1024);
//	ray_set.TestKDTree(mesh_kd_tree->GetTree());
//	ray_set.PrintResult();

	return Succeed;
}




PI_DECLARE_TEST(MeshKDTreeBuilderTest) {

#if !_DEBUG
	// Try to increase the priority for better tunning.
	piSetRealtimePriority();
#endif

		
	// Allocation.
	PiMeshPtr mesh( new PiMesh() );
	
	// Importing.
	//if( mesh->Load("meshes/hoppe/cathead.m") == NULL ) {
	if( mesh->Load("meshes/standford/bunny.ply") == NULL ) {
	//if( mesh->Load("meshes/hires/face.obj") == NULL ) {
	//if( mesh->Load("meshes/hires/triceratops.obj") == NULL ) {
	//if( mesh->Load("meshes/hires/nefertiti tired.obj") == NULL ) {
		return Failed;
	}
	
	PiMeshKDTreePtr mesh_kd_tree( new PiMeshKDTree(mesh.GetPtr()) );

	uint64 begin = piGetClockCounter();

	mesh_kd_tree->BuildTree();

	uint64 elapsed = piGetClockCounter() - begin;
	
	piDebug("KD Tree built in %.6fM clocks.\n", double(elapsed) / (1024*1024) );
		
	return Succeed;
}

#endif // 0
