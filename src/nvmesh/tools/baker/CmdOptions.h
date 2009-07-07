// Copyright NVIDIA Corporation 2008 -- Edgar Velazquez-Armendariz <edgarv@nvidia.com>

#ifndef NV_BAKER_CMDOPTIONS_H
#define NV_BAKER_CMDOPTIONS_H

#include <tclap/CmdLine.h>

namespace nv {

	/// Basic common protocol for command line options using TCLAP: the classes will provide
	/// a method which returns a fresh opaque object with the options and which provides
	/// a single method to add its options to a command line.
	struct CmdOptions {
		virtual void add(TCLAP::CmdLine & cmd) = 0;
	};

	class CmdOptionsProvider {
		// The responsible of handling the options object is the caller: the Options
		// provider only creates it.
		virtual CmdOptions * getCmdOptions() = 0;
		virtual void setCmdOptions(CmdOptions * opt) = 0;
	};


}	// namespace nv

#endif // NV_BAKER_CMDOPTIONS_H
