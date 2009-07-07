// Copyright NVIDIA Corporation 2006 -- Denis Kovacs <dkovacs@nvidia.com>

#ifndef NV_GEOMETRYIMAGE_H
#define NV_GEOMETRYIMAGE_H

#include <nvimage/FloatImage.h>
#include <nvimage/HoleFilling.h>

#include <nvmath/Vector.h>

namespace nv
{

	class GeometryImage
	{
	public:

		enum Flags
		{
			PositionFlag = 0x01,
			NormalFlag = 0x02,
			DisplacementFlag = 0x04,
			VectorDisplacementFlag = 0x08,
			CurvatureFlag = 0x10,
			OcclusionFlag = 0x20,
			BentNormalFlag = 0x40
		};

		GeometryImage(uint flags, Vector2 extents)
		{
			setupFlags(flags);

			m_geometryImage = new FloatImage();
			m_geometryImage->allocate(m_componentCount, (uint)extents.x(), (uint)extents.y());
			m_geometryImage->clear(0.0f);
		}

		// GeometryImage takes ownership of the given image.
		GeometryImage(uint flags, FloatImage * img)
		{
			setupFlags(flags);
			m_geometryImage = img;
		}

		~GeometryImage()
		{
			delete m_geometryImage;
		}

		inline FloatImage * img()
		{
			return m_geometryImage;
		}

		inline int width() const { return m_geometryImage->width(); }
		inline int height() const { return m_geometryImage->height(); }
		inline int componentNum() const { 
			return m_geometryImage->componentNum(); 
		}
		inline int coverageChannel() const { return m_coverageChannel; }
		inline int positionChannel(uint channel = 0) const { 
			nvDebugCheck(channel < 3); return m_positionChannel + channel; 
		}
		inline int normalChannel(uint channel = 0) const { 
			nvDebugCheck(channel < 3); return m_normalChannel + channel; 
		}
		inline int displacementChannel(uint channel = 0) const { 
			nvDebugCheck(channel < 3); return m_displacementChannel + channel; 
		}
		inline int curvatureChannel() const { return m_curvatureChannel; }
		inline int occlusionChannel() const { return m_occlusionChannel; }
		inline int bentNormalChannel(uint channel = 0) const { 
			nvDebugCheck(channel < 3); return m_bentNormalChannel + channel; 
		}

		inline const float pixel(int x, int y, int c) const
		{
			return m_geometryImage->pixel(x, y, c);
		}

		inline void setPixel(float f, int x, int y, int c)
		{
			m_geometryImage->setPixel(f, x, y, c);
		}

		inline void addPixel(float f, int x, int y, int c)
		{
			m_geometryImage->addPixel(f, x, y, c);
		}

		inline void setPixel(Vector3 v, int x, int y, int c)
		{
			m_geometryImage->setPixel(v.x(), x, y, c+0);
			m_geometryImage->setPixel(v.y(), x, y, c+1);
			m_geometryImage->setPixel(v.z(), x, y, c+2);
		}

		inline void addPixel(Vector3 v, int x, int y, int c)
		{
			m_geometryImage->addPixel(v.x(), x, y, c+0);
			m_geometryImage->addPixel(v.y(), x, y, c+1);
			m_geometryImage->addPixel(v.z(), x, y, c+2);
		}

		Vector3 position(int x, int y) const
		{
			return Vector3(
				pixel(x, y, positionChannel(0)),
				pixel(x, y, positionChannel(1)),
				pixel(x, y, positionChannel(2)));
		}

		Vector3 normal(int x, int y) const
		{
			return Vector3(
				pixel(x, y, normalChannel(0)),
				pixel(x, y, normalChannel(1)),
				pixel(x, y, normalChannel(2)));
		}

		Vector3 bentNormal(int x, int y) const
		{
			return Vector3(
				pixel(x, y, bentNormalChannel(0)),
				pixel(x, y, bentNormalChannel(1)),
				pixel(x, y, bentNormalChannel(2)));
		}

		GeometryImage * fastDownSample() const
		{
			return new GeometryImage( m_flags, 
				m_geometryImage->fastDownSample() );
		}

		// normalize pixels by coverage 
		void fillSubpixels(BitMap * bmap, float threshold)
		{
			nvDebugCheck(bmap != NULL);
			nvDebugCheck(threshold >= 0.0f && threshold < 1.0f);

			uint w = width();
			uint h = height();
			uint c = componentNum(); // do coverage last.
			uint coverageCh = coverageChannel();

			float * coverage = img()->channel(coverageCh);

			// clear insignificant coverage
			for (uint i = 0; i < h*w; i++) 
			{
				if (coverage[i] < threshold) {
					coverage[i] = 0.0f;
					bmap->clearBitAt(i);
				}
			}

			// normalize pixels (including coverage)
			for (uint k = 0; k <= c; k++) 
			{
				if (k == coverageCh) continue; // skip coverage until the end

				uint ch = k;
				if (ch == c) ch = coverageCh;      // normalize coverage a the end

				float * channel = img()->channel(ch);
				for (uint i = 0; i < h*w; i++) 
				{
					if (coverage[i] > 0.0f) {
						//nvDebugCheck(!isZero(channel[i]));
						channel[i] /= coverage[i];
					}
				}
			}
		}

		BitMap * getBitMaskFromCoverage() const
		{
			uint w = width();
			uint h = height();

			BitMap * bm = new BitMap(w, h);
			bm->clearAll();

			float * coverage = m_geometryImage->channel(coverageChannel());

			// clear insignificant coverage
			for (uint i = 0; i< h*w; i++) 
			{
				if (coverage[i] != 0.0) bm->setBitAt(i);
			}

			return bm;
		}

		void setNormalMap(FloatImage * nmap)
		{
			uint w = width();
			uint h = height();
			uint size = w * h;

			nvCheck(nmap->width() == w && nmap->height() == h);

			uint xoffset = size * normalChannel(0);
			uint yoffset = size * normalChannel(1);
			uint zoffset = size * normalChannel(2);
			uint coffset = size * coverageChannel();

			for (uint i = 0; i < size; i++)
			{
				Vector3 normal = Vector3(
					nmap->pixel(i),
					nmap->pixel(i + size),
					nmap->pixel(i + size * 2));

				float len = length(normal);

				if (len > 0.1f)
				{
					m_geometryImage->setPixel(normal.x(), i + xoffset);
					m_geometryImage->setPixel(normal.y(), i + yoffset);
					m_geometryImage->setPixel(normal.z(), i + zoffset);
					m_geometryImage->setPixel(clamp(len, 0.0f, 1.0f), i + coffset);
				}
			}
		}

	private:

		void setupFlags(uint flags)
		{
			m_flags = flags;
			m_componentCount = 1;

			m_coverageChannel = 0;
			m_positionChannel = -1;
			m_normalChannel = -1;
			m_displacementChannel = -1;
			m_curvatureChannel = -1;
			m_occlusionChannel = -1;
			m_bentNormalChannel = -1;

			if (flags & PositionFlag) 
			{
				m_positionChannel = m_componentCount;
				m_componentCount += 3;
			}

			if (flags & NormalFlag) 
			{
				m_normalChannel = m_componentCount;
				m_componentCount += 3;
			}

			if (flags & VectorDisplacementFlag) 
			{
				m_displacementChannel = m_componentCount; 
				m_componentCount += 3;
			}
			else if (flags & DisplacementFlag) 
			{
				m_displacementChannel = m_componentCount; 
				m_componentCount++;
			} 

			if (flags & CurvatureFlag)
			{
				m_curvatureChannel = m_componentCount;
				m_componentCount++;
			}

			if (flags & OcclusionFlag)
			{
				m_occlusionChannel = m_componentCount;
				m_componentCount++;
			}

			if (flags & NormalFlag) 
			{
				m_bentNormalChannel = m_componentCount;
				m_componentCount += 3;
			}
		}

	private:

		FloatImage * m_geometryImage;

		uint m_flags;
		uint m_componentCount;

		int m_coverageChannel;
		int m_positionChannel;
		int m_normalChannel;
		int m_displacementChannel;
		int m_curvatureChannel;
		int m_occlusionChannel;
		int m_bentNormalChannel;

	};

} // nv

#endif // NV_GEOMETRYIMAGE_H
