// This code is in the public domain -- castanyo@yahoo.es

#include <nvcore/StrLib.h>
#include <nvcore/Tokenizer.h>
//#include "Pi/Core/LinkList.h"

#include <nvmath/Vector.h>

#include "MeshImportPLY.h"

using namespace nv;

namespace {

	// File formats.
	enum PlyFormat {
		PLY_ASCII     = 1,		// ascii PLY file
		PLY_BINARY_BE = 2,		// binary PLY file, big endian
		PLY_BINARY_LE = 3		// binary PLY file, little endian
	};
	
	
	// Ply number types.
	enum PlyType {
		PLY_INT8	= 1,
		PLY_INT16	= 2,
		PLY_INT32	= 3,
		PLY_UINT8	= 4,
		PLY_UINT16	= 5,
		PLY_UINT32	= 6,
		PLY_FLOAT32	= 7,
		PLY_FLOAT64	= 8
	};
	
	
	// Ply property type.
	enum PlyElementType {
		PLET_VERTEX	= 0,
		PLET_FACE	= 1,
		PLET_OTHER	= 2
	};
	
	
	// Ply property type.
	enum PlyPropertyMode {
		PLPM_SCALAR	= 0,
		PLPM_LIST	= 1,
		PLPM_STRING	= 2
	};
	
	
	// Ply element.
	struct PlyProperty {
		char *			name;
		PlyPropertyMode	mode;
		PlyType			type;
		PlyType			ctype;		// Type of the count.
		PlyProperty *	next;
		PlyProperty *	prev;
	};
	
	
	// Ply element.
	struct PlyElement {
		char *			name;
		PlyElementType	type;
		uint32			num;
		PlyProperty		props;
		PlyElement *	next;
		PlyElement *	prev;
	};
	
	
	// Ply file header.
	struct PlyHeader {
		PlyFormat	format;
		uint32		version;
		PlyElement  elems;
	};



	// Read and parse a type name.
	static bool ReadType( PiParser &parser, int &type ) {
		if( !strCmp( parser.token, "Int8" ) )			type = PLY_INT8;
		else if( !strCmp( parser.token, "char" ) )	type = PLY_INT8;
		else if( !strCmp( parser.token, "Int16" ) )	type = PLY_INT16;
		else if( !strCmp( parser.token, "short" ) )	type = PLY_INT16;
		else if( !strCmp( parser.token, "Int32" ) )	type = PLY_INT32;
		else if( !strCmp( parser.token, "int" ) )		type = PLY_INT32;
		else if( !strCmp( parser.token, "Uint8" ) )	type = PLY_UINT8;
		else if( !strCmp( parser.token, "uchar" ) )	type = PLY_UINT8;
		else if( !strCmp( parser.token, "Uint16" ) )	type = PLY_UINT16;
		else if( !strCmp( parser.token, "ushort" ) )	type = PLY_UINT16;
		else if( !strCmp( parser.token, "Uint32" ) )	type = PLY_UINT32;
		else if( !strCmp( parser.token, "uint" ) )	type = PLY_UINT32;
		else if( !strCmp( parser.token, "Float32" ) )	type = PLY_FLOAT32;
		else if( !strCmp( parser.token, "float" ) )	type = PLY_FLOAT32;
		else if( !strCmp( parser.token, "Float64" ) )	type = PLY_FLOAT64;
		else if( !strCmp( parser.token, "double" ) )	type = PLY_FLOAT64;
		else return false;	
		return true;
	}
	
	
	// Read an element of the body of the PLY file.
	static bool ReadFloatValue( PiParser &parser, int format, int type, float &val ) {
	
		if( format == PLY_ASCII ) {
			if( !parser.GetToken( true ) ) return false;
			sscanf( parser.token, "%f", &val );
		}
		else {
			switch( type ) {
				case PLY_INT8: {
					sint8 i;
					if( !parser.SkipBytes( &i, 1 ) ) return false;
					val = (float) i;
				} break;
				case PLY_INT16: {
					sint16 i;
					if( !parser.SkipBytes( &i, 2 ) ) return false;
					else if( format==PLY_BINARY_BE ) i = piSwapBE16( i );
					else if( format==PLY_BINARY_LE ) i = piSwapLE16( i );
					val = (float) i;
				} break;
				case PLY_UINT8: {
					uint8 u;
					if( !parser.SkipBytes( &u, 1 ) ) return false;
					val = (float) u;
				} break;
				case PLY_UINT16: {
					uint16 u;
					if( !parser.SkipBytes( &u, 2 ) ) return false;
					else if( format==PLY_BINARY_BE ) u = piSwapBE16( u );
					else if( format==PLY_BINARY_LE ) u = piSwapLE16( u );
					val = (float) u;
				} break;
				case PLY_FLOAT32: {
					uint32 u;
					if( !parser.SkipBytes( &u, 4 ) ) return false;
					else if( format==PLY_BINARY_BE ) u = piSwapBE32( u );
					else if( format==PLY_BINARY_LE ) u = piSwapLE32( u );
					val = (float &)u;
				} break;
				case PLY_FLOAT64: {
					uint64 u;
					if( !parser.SkipBytes( &u, 8 ) ) return false;
					else if( format==PLY_BINARY_BE ) u = piSwapBE64( u );
					else if( format==PLY_BINARY_LE ) u = piSwapLE64( u );
					val = (float) (double &)u;
				} break;
				default:
					val = 0.0f;
					return false;
			}
		}
	
		return true;
	}
	
	
	/** Read an element of the body of the PLY file. */
	static bool ReadIntValue( PiParser &parser, int format, int type, uint &val ) {
	
		if( format == PLY_ASCII ) {
			if( !parser.GetToken( true ) ) return false;
			sscanf( parser.token, "%i", &val );
		}
		else {
			switch( type ) {
				case PLY_INT8: {
					sint8 i;
					if( !parser.SkipBytes( &i, 1 ) ) return false;
					val = uint(i);
				} break;
				case PLY_INT16: {
					uint16 u;
					if( !parser.SkipBytes( &u, 2 ) ) return false;
					else if( format==PLY_BINARY_BE ) u = piSwapBE16( u );
					else if( format==PLY_BINARY_LE ) u = piSwapLE16( u );
					val = (uint) (sint16 &)u;
				} break;
				case PLY_INT32: {
					uint32 u;
					if( !parser.SkipBytes( &u, 4 ) ) return false;
					else if( format==PLY_BINARY_BE ) u = piSwapBE32( u );
					else if( format==PLY_BINARY_LE ) u = piSwapLE32( u );
					val = (uint) (sint32 &)u;
				} break;
				case PLY_UINT8: {
					uint8 i;
					if( !parser.SkipBytes( &i, 1 ) ) return false;
					val = (uint)i;
				} break;
				case PLY_UINT16: {
					uint16 u;
					if( !parser.SkipBytes( &u, 2 ) ) return false;
					else if( format==PLY_BINARY_BE ) u = piSwapBE16( u );
					else if( format==PLY_BINARY_LE ) u = piSwapLE16( u );
					val = (uint) (uint16 &)u;
				} break;
				case PLY_UINT32: {
					uint32 u;
					if( !parser.SkipBytes( &u, 4 ) ) return false;
					else if( format==PLY_BINARY_BE ) u = piSwapBE32( u );
					else if( format==PLY_BINARY_LE ) u = piSwapLE32( u );
					val = (uint) (uint32 &)u;
				} break;
				default:
					val = 0;
					return false;
			}
		}
	
		return true;
	}
	

} // namespace


/*----------------------------------------------------------------------------
	Methods:
----------------------------------------------------------------------------*/

/** Load the PLY mesh. */
bool PiMeshImportPLY::Import(const char * name) {

	// Reset mesh.
	mesh->Reset();

	// Open file.
	PiAutoPtr<PiVirtualFile> vf( PiFileSystem::OpenFile( name ) );
	if( vf == NULL ) {
		piDebug( "*** Cannot open PLY file: '%s'\n", name );
		return false;
	}

	piDebug("--- Loading PLY file '%s'.\n", name );	
	
	PiParser parser;
	parser.StartParseBuffer( vf->GetMem(), vf->GetSize() );

	if( !parser.GetToken( true ) || strCmp(parser.token, "ply") ) {
		piDebug( "*** File '%s' does not appear to be a ply file.\n", name );
		return false;
	}

	PlyHeader head;
	PlyElement * current = NULL;
	LL_Reset( &head.elems, next, prev );


	bool result = false;

	while(true) {

		if( !parser.GetToken( true ) ) {
			piDebug( "*** Unexpected 'end of file' found while reading ply file '%s'.\n", name );
			break;
		}

		if( !strCmp( parser.token, "format" ) ) {
			// 'format' <type> <version>
			if( !parser.GetToken( false ) ) {
				piDebug( "*** Type expected while reading ply file '%s'.\n", name );
				break;
			}

			if( !strCmp( parser.token, "ascii" ) ) head.format = PLY_ASCII;
			else if( !strCmp( parser.token, "binary_big_endian" ) ) head.format = PLY_BINARY_BE;
			else if( !strCmp( parser.token, "binary_little_endian" ) ) head.format = PLY_BINARY_LE;

			if( !parser.GetToken( false ) ) {
				piDebug( "*** Version expected while reading ply file '%s'.\n", name );
				break;
			}

			head.version = 0;
		}
		else if( !strCmp( parser.token, "element" ) ) {
			// 'element' <type> <num>

			current = new PlyElement;
			LL_Reset( current, next, prev );
			LL_AddLast( &head.elems, current, next, prev );
			LL_Reset( &current->props, next, prev );

			if( !parser.GetToken( false ) ) {
				piDebug( "*** Element type expected while reading ply file '%s'.\n", name );
				break;
			}
			current->name = piStrDup( parser.token );
			if( !strCmp( parser.token, "vertex" ) ) current->type = PLET_VERTEX;
			else if( !strCmp( parser.token, "face" ) ) current->type = PLET_FACE;
			else current->type = PLET_OTHER;

			if( !parser.GetToken( false ) ) {
				piDebug( "*** Number of elements expected while reading ply file '%s'.\n", name );
				break;
			}
			current->num = atoi( parser.token );

		}
		else if( !strCmp( parser.token, "property" ) ) {
			// 'property' <type> <name>

			PlyProperty * prop = new PlyProperty;
			LL_Reset( prop, next, prev );
			LL_AddLast( &current->props, prop, next, prev );

			if( !parser.GetToken( false ) ) {
				piDebug( "*** Property type expected while reading ply file '%s'.\n", name );
				break;
			}

			if( !strCmp( parser.token, "list" ) ) {
				prop->mode = PLPM_LIST;
				if( !parser.GetToken( false ) ) {
					piDebug( "*** Property type expected while reading ply file '%s'.\n", name );
					break;
				}
				if( !ReadType( parser, (int &)prop->ctype ) ) {
					piDebug( "*** Unknown property type while reading ply file '%s'.\n", name );
					break;
				}

				if( !parser.GetToken( false ) ) {
					piDebug( "*** Property type expected while reading ply file '%s'.\n", name );
					break;
				}
			}
			else {
				prop->mode = PLPM_SCALAR;
			}

			if( !ReadType( parser, (int &)prop->type ) ) {
				piDebug( "*** Unknown property type while reading ply file '%s'.\n", name );
				break;
			}

			// read name
			if( !parser.GetToken( false ) ) {
				piDebug( "*** Property name expected while reading ply file '%s'.\n", name );
				break;
			}
			prop->name = piStrDup( parser.token );
		}
		else if( !strCmp( parser.token, "comment" ) ) {
			AutoString comment;
			while( parser.GetToken( false ) ) {
				comment.Append(" ");				
				comment.Append(parser.token);
			};
			piDebug( "---   Comment:%s\n", comment.GetStr() );
		}
		else if( !strCmp( parser.token, "end_header" ) ) {
			if( parser.NextLine() )
				result = true;
			break;
		}
	};


	if( result == false ) {
		piDebug("***   Load failed.\n");
	}
	else {

		// Analize what we have:
		PlyElement * e;
		for( e = head.elems.next; e != &head.elems; e = e->next ) {
			if( strCmp( e->name, "vertex" ) == 0 ) {
				mesh->SetVertexNum( e->num );
				//pos_array.Resize( e->num );
				//col_array.Resize( e->num );
			}
			else if( strCmp( e->name, "face" ) == 0 ) {
				// Number of faces is not the number of triangles.
				//face_array.Resize( e->num );
			}
		}


		PiMesh::Channel * pos_channel = mesh->GetChannel(mesh->CreateChannel("position", VS_POS, VF_FLOAT, 3));

		
		
		uint pos_num = 0;
		Vec2 tmp(0, 0);

		// Extract data.
		for( e = head.elems.next; e != &head.elems; e = e->next ) {

			for( uint i = 0; i < e->num; i++ ) {

				PlyProperty * p;
				for( p = e->props.next; p != &e->props; p = p->next ) {

					if( p->mode == PLPM_SCALAR ) {
						float value;
						ReadFloatValue( parser, head.format, p->type, value );

						if( strCmp( p->name, "x" ) == 0 ) {
							tmp.x = value;
						}
						else if( strCmp( p->name, "y" ) == 0 ) {
							tmp.y = value;
						}
						else if( strCmp( p->name, "z" ) == 0 ) {
							pos_channel->data[pos_num++] = Vec4(tmp, value, 1.0f);
						}
					}
					else if( p->mode == PLPM_LIST ) {
						uint count;
						ReadIntValue( parser, head.format, p->ctype, count );

						uint v0, v1, v2;
						ReadIntValue( parser, head.format, p->type, v0 );
						ReadIntValue( parser, head.format, p->type, v1 );

						// conver poly to fan
						for( uint f = 0; f < count-2; f++ ) {
							ReadIntValue( parser, head.format, p->type, v2 );
							mesh->AddFace(v0, v1, v2);
							v1 = v2;
						}
					}
				}
			}
		}

		piDebug("---   Load succeed.\n");
		piDebug( "---   %d vertices\n", mesh->GetVertexNum() );
		piDebug( "---   %d faces\n", mesh->GetFaceNum() );
	}

	// clean header
	PlyElement * e, * enext;
	for( e=head.elems.next; e!=&head.elems; e=enext ) {
		enext = e->next;
		LL_Reset( e, next, prev );

		PlyProperty * p, * pnext;
		for( p=e->props.next; p!=&e->props; p=pnext ) {
			pnext = p->next;
			LL_Reset( p, next, prev );

			mem::free( p->name );
			delete p;
		}

		mem::free( e->name );
		delete e;
	}

	return result;
}

