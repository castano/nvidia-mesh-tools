// todo

// Emil Person (Humus) also has code to load T3D files.

#if 0

// ==================================================================
// Copyright (c) 2003 Bart Sekura 
//					  sekurb@yahoo.com
//
// $Id: T3DMap.cpp,v 1.7 2003/06/13 06:50:16 bs Exp $
//
// ==================================================================

#include "app/T3DMap.h"
#include "app/Light.h"
#include "app/Surface.h"
#include "app/Mesh.h"
#include "app/ASE.h"
#include "app/LWO.h"
#include "app/T3D.h"
#include "base/Array.h"
#include "base/StringObject.h"
#include "base/Parser.h"
#include "base/Global.h"
#include "math/Vec3.h"
#include "math/Vec2.h"
#include "math/Matrix.h"
#include "math/Polygon.h"
#include "render/extgl.h"
#include "render/Shader.h"
#include "render/RenderObject.h"
#include "render/RenderDevice.h"

#define MAX_T3D_POLY_VERTNUM    (64)
#define MAX_T3D_MATERIAL_NAME   (256)

class T3DPolygon : public ConvexPolygon {
public:
	T3DPolygon() : ConvexPolygon(), next(0), material("$notex") {}

	T3DPolygon*		next;
    int				flags;
    int				link;
	Vec3			origin;
	Vec3			normal;
	Vec3			tu,tv;
    Vec2			pan;
	String			material;
};

// -------------------------
// misc util local functions
// -------------------------

static const char* vec3Names[] = { "X","Y","Z" };
static const char* rot3Names[] = { "Pitch","Yaw","Roll" };

static inline void swapVec3(Vec3& v)
{
    float a;

    a = v[1];
    v[1] = v[2];
    v[2] = a;
}

static inline Vec3 fixPos(Vec3& p)
{
	swapVec3(p);
	return p;
}

static inline Vec3 fixDir(Vec3& p)
{
	swapVec3(p);
	return -p;
}

static int parseName(Parser& parser, char* name)
{
	if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		return 0;
	if (!parser.nextToken())
		return 0;

	strcpy(name,parser.getTokenString());
	return 1;
}

static int parseVec3(Parser& parser, Vec3& v)
{
	for (int i=0; i<3; ++i) {
		if (!parser.nextToken())
			return 0;

		v[i] = (float) atof(parser.getTokenString());
		if (i<2) {
			if (!parser.nextToken() || strcmp(parser.getTokenString(),","))
				return 0;
		}
	}

	return 1;
}

static int parseTriple(Parser& parser, const char* names[], Vec3& v)
{
    // opening '('
	if (!parser.nextToken() || strcmp(parser.getTokenString(),"("))
		return 0;

    while (parser.nextToken()) {
        int i;
        const char* token;

        token = parser.getTokenString();
        if (!strcmp(token,")"))
            break;

        if (!strcmp(token,names[0]))
            i = 0;
        else if (!strcmp(token,names[1]))
            i = 1;
        else if (!strcmp(token,names[2]))
            i = 2;
        else
            continue;

        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
            return 0;
        if (!parser.nextToken())
            return 0;

        v[i] = (float) atof(parser.getTokenString());
    }

    return 1;
}

static void convertName(const char* s, const char* prefix, String& out)
{
	const char* dot = strrchr(s,'.');
	out = prefix;
	out += dot ? dot + 1 : s;
}

// ---------------------
// T3DMap implementation
// ---------------------

T3DMap::T3DMap()
:	Scene()
,	polyList(0), surfList(0)
,	worldMesh(0)
{
}

T3DMap::~T3DMap()
{
}

int T3DMap::parsePolygon(Parser& parser, const Vec3& origin)
{
	T3DPolygon poly;
	poly.pan = Vec2(0,0);

    if (!parser.nextToken()) 
		return 0;

	while (1) {
		const char* token;

		token = parser.getTokenString();
		if (!strcmp("End",token)) {
			if (!parser.nextToken())
				return 0;

			token = parser.getTokenString();
			if (!strcmp("Polygon",token))
				break;
		}
		else if (!strcmp("Vertex",token)) {

			Vec3 p;
			if (!parseVec3(parser,p))
				return 0;

			poly.verts.add(fixPos(p+origin));
		}
		else if (!strcmp("Origin",token)) {
			Vec3 p;
			if (!parseVec3(parser,p))
				return 0;

            poly.origin = fixPos(p+origin);
		}
		else if (!strcmp("Normal",token)) {
			Vec3 n;
			if (!parseVec3(parser,n))
				return 0;

            poly.normal = fixDir(n);
		}
		else if (!strcmp("TextureU",token)) {
			Vec3 t;
			if (!parseVec3(parser,t))
				return 0;

			poly.tu = fixDir(t);
		}
		else if (!strcmp("TextureV",token)) {
			Vec3 t;
			if (!parseVec3(parser,t))
				return 0;

            poly.tv = fixDir(t);
		}
        else if (!strcmp("Pan",token)) {

            for (int i=0; i<2; ++i) {
                int j;

                if (!parser.nextToken())
                    return 0;

                token = parser.getTokenString();
                if (!strcmp("U",token)) {
                    j = 0;
                } else if (!strcmp("V",token)) {
                    j = 1;
                } else {
                    break;
                }

	            if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		            return 0;
                if (!parser.nextToken())
                    return 0;

				poly.pan[j] = (float) atof(parser.getTokenString());
            }

            // no need to get next token
            continue;
        }
        else if (!strcmp("Texture",token)) {
	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		        return 0;
            if (!parser.nextToken())
                return 0;

			convertName(parser.getTokenString(),"base/",poly.material);
		}
        else if (!strcmp("Link",token)) {
	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		        return 0;
            if (!parser.nextToken())
                return 0;

            poly.link = atoi(parser.getTokenString());
        }

        if (!parser.nextToken())
            break;
	}

	if (poly.verts.getSize()) {

		T3DPolygon* p = new T3DPolygon(poly);
		p->next = polyList;
		polyList = p;
	}

	return 1;	
}

int T3DMap::parsePolyList(Parser& parser, const Vec3& origin)
{
	while (parser.nextToken()) {
		const char* token;

		token = parser.getTokenString();
		if (!strcmp("Begin",token)) {
			if (!parser.nextToken())
				return 0;

			token = parser.getTokenString();
			if (!strcmp("Polygon",token)) {
				if (!parsePolygon(parser,origin))
					return 0;
			}
		}
		else if (!strcmp("End",token)) {
			if (!parser.nextToken())
				return 0;

			token = parser.getTokenString();
			if (!strcmp("PolyList",token))
				break;
		}
	}

	return 1;
}

int T3DMap::parseBrush(Parser& parser)
{
	int brushDef;
	int buildBrush;
	Vec3 location(0,0,0);
	char name[256] = "<NULL>";

	brushDef = 0;
	buildBrush = 0;

	while (parser.nextToken()) {
		const char* token;

		token = parser.getTokenString();

		if (!strcmp("Begin",token)) {
			if (!parser.nextToken())
				return 0;

			token = parser.getTokenString();
			if (!strcmp("Brush",token))
				brushDef = 1;
			else if (!strcmp("PolyList",token)) {
				if (!brushDef)
					return 0;

				if (!buildBrush)
					continue;

				if (!parsePolyList(parser,location))
					return 0;
			}
		} else if (!strcmp("End",token)) {
			if (!parser.nextToken())
				return 0;

			token = parser.getTokenString();
			if (!strcmp("Actor",token))
				break;
			else if (!strcmp("Brush",token))
				brushDef = 0;

		} else if (!strcmp("Location",token)) {
			if (buildBrush) {
				if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
					return 0;
				if (!parseTriple(parser,vec3Names,location))
					return 0;

				location = fixPos(location);
			}
		} else if (!strcmp("Name",token)) {
			if (brushDef) {
				if (!parseName(parser,name))
					return 0;			

				if (!strcmp("Brush",name))
					buildBrush = 1;
			}
		}
	}

	GLog.print("Brush: %s\n",name);
	return 1;
}

int T3DMap::parseLight(Parser& parser)
{
	char name[256] = "<NULL>";
    Vec3 pos,rot;
    float fov,range;
    Matrix spotMat;

    pos = Vec3(0,0,0);
    rot = Vec3(0,0,0);
    fov = 90;
    range = 1024;

	while (parser.nextToken()) {
		const char* token;

		token = parser.getTokenString();
		if (!strcmp("End",token)) {
			if (!parser.nextToken())
				return 0;
			if (!strcmp("Actor",parser.getTokenString()))
				break;
		} else if (!strcmp("Name",token)) {
			if (!parseName(parser,name))
				return 0;			

		} else if (!strcmp("Location",token)) {
	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		        return 0;
            if (!parseTriple(parser,vec3Names,pos))
                return 0;

            pos = fixPos(pos);

        } else if (!strcmp("Rotation",token)) {

	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		        return 0;
            if (!parseTriple(parser,rot3Names,rot))
                return 0;

            for (int i=0; i<3; ++i)
                rot[i] = (float)((int)rot[i] % 65536) * 0.0054931640625f;
        }
	}

	spotMat.setIdentity();
	spotMat.rotate(-rot[0],Vec3(1,0,0));
	spotMat.rotate(rot[1]+90,Vec3(0,1,0));
	spotMat.rotate(rot[2],Vec3(0,0,1));
	spotMat.translate(-pos);

	Light* light = new Light();
	light->init(spotMat,fov,range);

	light->next = lightList;
	lightList = light;

    return 1;
}

int T3DMap::parseStaticMesh(Parser& parser)
{
	String name;
    Vec3 pos(0,0,0),rot(0,0,0),scale(1,1,1);

	while (parser.nextToken()) {
		const char* token;

		token = parser.getTokenString();
		if (!strcmp("End",token)) {
			if (!parser.nextToken())
				return 0;
			if (!strcmp("Actor",parser.getTokenString()))
				break;
		}
		else if (!strcmp("StaticMesh",token)) {
	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		        return 0;
	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"StaticMesh"))
		        return 0;
            if (!parser.nextToken())
                return 0;

			convertName(parser.getTokenString(),"data/models/",name);

		}
        else if (!strcmp("ColLocation",token)) {
	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		        return 0;
            if (!parseTriple(parser,vec3Names,pos))
                return 0;

            pos = fixPos(pos);
        }
        else if (!strcmp("Rotation",token)) {

	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		        return 0;
            if (!parseTriple(parser,rot3Names,rot))
                return 0;

            for (int i=0; i<3; ++i)
                rot[i] = (float)((int)rot[i] % 65536) * 0.0054931640625f;

			//rot = fixDir(rot);
        }
        else if (!strcmp("DrawScale3D",token)) {

	        if (!parser.nextToken() || strcmp(parser.getTokenString(),"="))
		        return 0;
            if (!parseTriple(parser,vec3Names,scale))
                return 0;
		}
	}

	GLog.print("staticMesh: %s (%.2f %.2f %.2f)\n",name.asChar(),
		pos.x,pos.y,pos.z);

	Matrix m;
	m.setIdentity();
	/*
	m.rotate(rot[0],Vec3(1,0,0));
	m.rotate(rot[1],Vec3(0,1,0));
	m.rotate(rot[2],Vec3(0,0,1));
	m.translate(-pos);
	m.scale(scale);
	*/

	m.translate(pos);
	m.rotate(90,Vec3(1,0,0));
	m.rotate(180,Vec3(0,0,1));
	m.rotate(rot[0],Vec3(0,1,0));
	m.rotate(rot[1],Vec3(0,0,1));
	m.rotate(rot[2],Vec3(1,0,0));
	m.scale(scale);

	//m.rotate(90,Vec3(1,0,0));
	//m.rotate(180,Vec3(0,0,1));

	/*
	MeshFileASE ase(m);
	name += ".ase";
	ase.load(meshList,name.asChar());

	return 1;
	*/

	MeshFileASE ase(m);
	if (ase.load(meshList,(name+".ase").asChar()))
		return 1;

	MeshFileLWO lwo(m);
	if (lwo.load(meshList,(name+".lwo").asChar()))
		return 1;

	GLog.print("ERROR: static mesh \"%s\" NOT FOUND\n",name.asChar());
	return 0;
}

int T3DMap::parseActor(Parser& parser)
{
    if (!parser.expectToken("Class"))
        return 0;
    if (!parser.expectToken("="))
        return 0;
    if (!parser.nextToken())
        return 0;

    const char* token = parser.getTokenString();
	if (!strcmp("Brush",token))
		return parseBrush(parser);
	else if (!strcmp("Light",token) || !strcmp("Spotlight",token))
		return parseLight(parser);
	else if (!strcmp("StaticMeshActor",token))
		return parseStaticMesh(parser);

	// parse dummy actor
	while (parser.nextToken()) {
		if (!strcmp("End",parser.getTokenString())) {
			if (!parser.nextToken())
				return 0;
			if (!strcmp("Actor",parser.getTokenString()))
				break;
		}
	}
    return 1;
}

int T3DMap::load(const char* fileName)
{
    // parse map file
    Parser parser;
    parser.setDelims("()=,");
    if (parser.parseFile(fileName)) {

		close();
		GLog.print("parsing map: \"%s\"\n",fileName);
        int inMap = 0;

        while (parser.nextToken()) {
			const char* token = parser.getTokenString();			
            if (!strcmp("Begin",token)) {
				if (!parser.nextToken()) {
					GLog.print("ERROR: token expected after 'Begin'\n");
					break;
				}

				token = parser.getTokenString();
				if (!strcmp("Map",token)) {
					inMap = 1;
					continue;
				}

				if (!strcmp("Actor",token)) {
					if (!inMap) {
						GLog.print("WARNING: Actor outside map\n");
						continue;
					}

					if (!parseActor(parser))
						break;					
				}
			}
			else if (!strcmp("End",token)) {
				if (!parser.nextToken()) {
					GLog.print("ERROR: token expected after 'End'\n");
					break;
				}

				if (!strcmp("Map",parser.getTokenString())) {
					if (!inMap)
						GLog.print("ERROR: 'End' without matching 'Begin'\n");

					break;
				}
			}
		}

		GLog.print("done!\n\n");
		for (Mesh* mesh = meshList; mesh; mesh=mesh->getNext())
			GLog.print("mesh: %s\n",mesh->getName());
    }

    return 1;
}

void T3DMap::createSurfs(ShaderManager* shaderMan)
{
    T3DPolygon* poly = polyList;
    while (poly) {
		Surface* s = new Surface();
		s->verts.resize(poly->verts.getSize());

		Shader* shader = shaderMan->createShader(poly->material.asChar());

		int j;
		for (j=0; j<matGroups.getSize(); ++j) {
			MaterialGroup& g = matGroups[j];
			if (g.shader == shader) {
				s->nextGroup = g.surfList;
				g.surfList = s;
				break;
			}
		}

		if (j == matGroups.getSize()) {
			MaterialGroup g;
			g.name = poly->material;
			g.surfList = s;
			g.shader = shader;

			matGroups.add(g);
		}

		float width = (float) (shader->getDiffuseMap()->getWidth()),
			  height = (float) (shader->getDiffuseMap()->getHeight());

		for (int i=0; i<poly->verts.getSize(); ++i) {
			s->verts[i].xyz = poly->verts[i];
			s->verts[i].normal = poly->normal;

			Vec3 a = poly->verts[i] - poly->origin;
			float u = a * poly->tu,
				  v = a * poly->tv;

			s->verts[i].st[0] = (u + poly->pan[0]) / width;
			s->verts[i].st[1] = (v + poly->pan[1]) / height;
		}

		s->next = surfList;
		surfList = s;

		poly = poly->next;
	}

	Surface::makeTris(surfList);
	worldMesh = new MeshT3D(matGroups);

	for (Mesh* mesh=meshList; mesh; mesh=mesh->getNext())
		mesh->setupShaders(shaderMan,"models/");

	worldMesh->next = meshList;
	meshList = worldMesh;

	cleanup();
}

void T3DMap::cleanup()
{
	// delete surfs
	for (Surface* s=surfList; s;) {
		Surface* surf=s;
		s=s->next;
		delete surf;
	}

	surfList = 0;

	// delete polys
	for (T3DPolygon* p=polyList; p;) {
		T3DPolygon* poly = p;
		p = p->next;
		delete poly;
	}

	polyList = 0;

	matGroups.clear();
}

void T3DMap::draw(RenderDevice* rdev)
{
	/*
    glColor3f(1,0,0);
    glBegin(GL_LINES);

    T3DPolygon* poly = polyList;
    while (poly) {
        for (int i=0; i<poly->verts.getSize(); ++i) {
            glVertex3fv(poly->verts[i]);
            glVertex3fv(poly->verts[(i+1)%poly->verts.getSize()]);
        }

        poly = poly->next;
    }
    glEnd();
	*/

	/*
	glColor3f(1,1,1);
	Surface* s = surfList;
	while (s) {
		rdev->setTexture(0,s->shader->getDiffuseMap());
		glBegin(GL_POLYGON);
		for (int i=0; i<s->verts.getSize(); ++i) {
			const Vertex& v = s->verts[i];
			glTexCoord2fv(v.st);
			glVertex3fv(v.xyz);
		}
		glEnd();

		s = s->next;
	}
	*/

	glColor3f(1,1,1);

	/*
	for (int i=0; i<matGroups.getSize(); ++i) {
		MaterialGroup& g = matGroups[i];
		rdev->setTexture(0,g.shader->getDiffuseMap());

		for (Surface* s=g.surfList; s; s=s->nextGroup) {
			glBegin(GL_POLYGON);
			for (int j=0; j<s->verts.getSize(); ++j) {
				const Vertex& v = s->verts[j];
				glTexCoord2fv(v.st);
				glVertex3fv(v.xyz);
			}
			glEnd();
		}
	}
	*/
	//worldMesh->drawDebug(rdev,Mesh::EDrawLines);

	rdev->disableAllStages();
	//glMatrixMode(GL_MODELVIEW);
	for (Mesh* mesh = meshList; mesh; mesh=mesh->getNext()) {
		//glPushMatrix();
		//glMultMatrixf(mesh->getTransform().x[0]);

		mesh->drawDebug(rdev,Mesh::EDrawDiffuse);

		//glPopMatrix();
	}

	glColor3f(1,0,0);
	for (Light* light = lightList; light; light=light->getNext()) {
		light->volume.draw();
	}
}

#endif
