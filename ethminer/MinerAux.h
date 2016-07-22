#pragma once

/*
	This file is part of cpp-ethereum.

	cpp-ethereum is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	cpp-ethereum is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with cpp-ethereum.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @file MinerAux.cpp
 * @author Gav Wood <i@gavwood.com>
 * @date 2014
 * CLI module for mining.
 */

#include <boost/algorithm/string/case_conv.hpp>

#include <libdevcore/CommonJS.h>
#include <libevmcore/Instruction.h>
#include <libethcore/BasicAuthority.h>
#include <libethcore/Exceptions.h>
#include <libethashseal/EthashAux.h>
#include <libethashseal/EthashGPUMiner.h>
#include <libethashseal/EthashCPUMiner.h>
#include <libethashseal/Ethash.h>

#if ETH_ETHASHCL
#include <libethash-cl/ethash_cl_miner.h>
#endif // ETH_ETHASHCL

#if ETH_JSONRPC
#include <jsonrpccpp/server/connectors/httpserver.h>
#include <jsonrpccpp/client/connectors/httpclient.h>
#endif // ETH_JSONRPC

#include "cpp-ethereum/BuildInfo.h"

#if ETH_JSONRPC
#include "PhoneHome.h"
#include "FarmClient.h"
#endif // ETH_JSONRPC

#undef RETURN

bool isTrue(std::string const& _m)
{
	return _m == "on" || _m == "yes" || _m == "true" || _m == "1";
}

bool isFalse(std::string const& _m)
{
	return _m == "off" || _m == "no" || _m == "false" || _m == "0";
}

inline std::string credits()
{
	std::ostringstream out;
	out
		<< "cpp-ethereum " << dev::Version << std::endl
		<< "  By cpp-ethereum contributors, (c) 2013-2016." << std::endl
		<< "  See the README for contributors and credits." << std::endl;
	return out.str();
}

class BadArgument: public dev::Exception {};
struct MiningChannel: public dev::LogChannel
{
	static const char* name() { return EthGreen "miner"; }
	static const int verbosity = 2;
	static const bool debug = false;
};
#define minelog clog(MiningChannel)

class MinerCLI
{
public:
	enum class OperationMode
	{
		None,
		DAGInit,
		Benchmark,
		Farm
	};


	MinerCLI(OperationMode _mode = OperationMode::None): mode(_mode) {
		dev::eth::Ethash::init();
		dev::eth::NoProof::init();
		dev::eth::BasicAuthority::init();
	}

	bool interpretOption(int& i, int argc, char** argv)
	{
		std::string arg = argv[i];
		if ((arg == "-F" || arg == "--farm") && i + 1 < argc)
		{
			mode = OperationMode::Farm;
			m_farmURL = argv[++i];
		}
		else if (arg == "--farm-recheck" && i + 1 < argc)
			try {
				m_farmRecheckPeriod = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		else if (arg == "--opencl-platform" && i + 1 < argc)
			try {
				m_openclPlatform = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		else if (arg == "--opencl-device" && i + 1 < argc)
			try {
				m_openclDevice = std::stol(argv[++i]);
				m_miningThreads = 1;
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
#if ETH_ETHASHCL
		else if (arg == "--cl-global-work" && i + 1 < argc)
			try {
				m_globalWorkSizeMultiplier = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		else if (arg == "--cl-local-work" && i + 1 < argc)
			try {
				m_localWorkSize = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		else if (arg == "--cl-ms-per-batch" && i + 1 < argc)
			try {
				m_msPerBatch = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
#endif // ETH_ETHASHCL
		else if (arg == "--list-devices")
			m_shouldListDevices = true;
		else if (arg == "--allow-opencl-cpu")
			m_clAllowCPU = true;
		else if (arg == "--cl-extragpu-mem" && i + 1 < argc)
			m_extraGPUMemory = 1000000 * std::stol(argv[++i]);
		else if (arg == "--phone-home" && i + 1 < argc)
		{
			std::string m = argv[++i];
			if (isTrue(m))
				m_phoneHome = true;
			else if (isFalse(m))
				m_phoneHome = false;
			else
			{
				std::cerr << "Bad " << arg << " option: " << m << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		}
		else if (arg == "--benchmark-warmup" && i + 1 < argc)
			try {
				m_benchmarkWarmup = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		else if (arg == "--benchmark-trial" && i + 1 < argc)
			try {
				m_benchmarkTrial = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		else if (arg == "--benchmark-trials" && i + 1 < argc)
			try {
				m_benchmarkTrials = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		else if (arg == "-C" || arg == "--cpu")
			m_minerType = "cpu";
		else if (arg == "-G" || arg == "--opencl")
			m_minerType = "opencl";
		else if (arg == "--current-block" && i + 1 < argc)
			m_currentBlock = std::stol(argv[++i]);
		else if (arg == "--no-precompute")
		{
			m_precompute = false;
		}
		else if ((arg == "-D" || arg == "--create-dag") && i + 1 < argc)
		{
			std::string m = boost::algorithm::to_lower_copy(std::string(argv[++i]));
			mode = OperationMode::DAGInit;
			try
			{
				m_initDAG = std::stol(m);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << m << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		}
		else if ((arg == "-w" || arg == "--check-pow") && i + 4 < argc)
		{
			std::string m;
			try
			{
				dev::eth::BlockHeader bi;
				m = boost::algorithm::to_lower_copy(std::string(argv[++i]));
				dev::h256 powHash(m);
				m = boost::algorithm::to_lower_copy(std::string(argv[++i]));
				dev::h256 seedHash;
				if (m.size() == 64 || m.size() == 66)
					seedHash = dev::h256(m);
				else
					seedHash = dev::eth::EthashAux::seedHash(std::stol(m));
				m = boost::algorithm::to_lower_copy(std::string(argv[++i]));
				bi.setDifficulty(dev::u256(m));
				auto boundary = dev::eth::Ethash::boundary(bi);
				m = boost::algorithm::to_lower_copy(std::string(argv[++i]));
				dev::eth::Ethash::setNonce(bi, dev::h64(m));
				auto r = dev::eth::EthashAux::eval(seedHash, powHash, dev::h64(m));
				bool valid = r.value < boundary;
				std::cout << (valid ? "VALID :-)" : "INVALID :-(") << std::endl;
				std::cout << r.value << (valid ? " < " : " >= ") << boundary << std::endl;
				std::cout << "  where " << boundary << " = 2^256 / " << bi.difficulty() << std::endl;
				std::cout << "  and " << r.value << " = ethash(" << powHash << ", " << dev::h64(m) << ")" << std::endl;
				std::cout << "  with seed as " << seedHash << std::endl;
				if (valid)
					std::cout << "(mixHash = " << r.mixHash << ")" << std::endl;
				std::cout << "SHA3( light(seed) ) = " << sha3(dev::eth::EthashAux::light(dev::eth::Ethash::seedHash(bi))->data()) << std::endl;
				exit(0);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << m << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		}
		else if (arg == "-M" || arg == "--benchmark")
			mode = OperationMode::Benchmark;
		else if ((arg == "-t" || arg == "--mining-threads") && i + 1 < argc)
		{
			try {
				m_miningThreads = std::stol(argv[++i]);
			}
			catch (...)
			{
				std::cerr << "Bad " << arg << " option: " << argv[i] << std::endl;
				BOOST_THROW_EXCEPTION(BadArgument());
			}
		}
                else if (arg == "--disable-submit-hashrate")
                        m_submitHashrate = false;
		else
			return false;
		return true;
	}

	void execute()
	{
		if (m_shouldListDevices)
		{
#if ETH_ETHASHCL
			dev::eth::EthashGPUMiner::listDevices();
#endif // ETH_ETHASHCL

			exit(0);
		}

		if (m_minerType == "cpu")
			dev::eth::EthashCPUMiner::setNumInstances(m_miningThreads);
		else if (m_minerType == "opencl")
		{
#if ETH_ETHASHCL
			if (!dev::eth::EthashGPUMiner::configureGPU(
					m_localWorkSize,
					m_globalWorkSizeMultiplier,
					m_msPerBatch,
					m_openclPlatform,
					m_openclDevice,
					m_clAllowCPU,
					m_extraGPUMemory,
					m_currentBlock
				))
				exit(1);
			dev::eth::EthashGPUMiner::setNumInstances(m_miningThreads);
#else
			cerr << "Selected GPU mining without having compiled with -DETHASHCL=1" << endl;
			exit(1);
#endif // ETH_ETHASHCL
		}
		if (mode == OperationMode::DAGInit)
			doInitDAG(m_initDAG);
		else if (mode == OperationMode::Benchmark)
			doBenchmark(m_minerType, m_phoneHome, m_benchmarkWarmup, m_benchmarkTrial, m_benchmarkTrials);
		else if (mode == OperationMode::Farm)
			doFarm(m_minerType, m_farmURL, m_farmRecheckPeriod);
	}

	static void streamHelp(std::ostream& _out)
	{
		_out
#if ETH_JSONRPC
			<< "Work farming mode:" << std::endl
			<< "    -F,--farm <url>  Put into mining farm mode with the work server at URL (default: http://127.0.0.1:8545)" << std::endl
			<< "    --farm-recheck <n>  Leave n ms between checks for changed work (default: 500)." << std::endl
			<< "    --no-precompute  Don't precompute the next epoch's DAG." << std::endl
#endif // ETH_JSONRPC
			<< "Ethash verify mode:" << std::endl
			<< "    -w,--check-pow <headerHash> <seedHash> <difficulty> <nonce>  Check PoW credentials for validity." << std::endl
			<< std::endl
			<< "Benchmarking mode:" << std::endl
			<< "    -M,--benchmark  Benchmark for mining and exit; use with --cpu and --opencl." << std::endl
			<< "    --benchmark-warmup <seconds>  Set the duration of warmup for the benchmark tests (default: 3)." << std::endl
			<< "    --benchmark-trial <seconds>  Set the duration for each trial for the benchmark tests (default: 3)." << std::endl
			<< "    --benchmark-trials <n>  Set the number of trials for the benchmark tests (default: 5)." << std::endl
#if ETH_JSONRPC
			<< "    --phone-home <on/off>  When benchmarking, publish results (default: on)" << std::endl
#endif // ETH_JSONRPC
			<< "DAG creation mode:" << std::endl
			<< "    -D,--create-dag <number>  Create the DAG in preparation for mining on given block and exit." << std::endl
			<< "Mining configuration:" << std::endl
			<< "    -C,--cpu  When mining, use the CPU." << std::endl
			<< "    -G,--opencl  When mining use the GPU via OpenCL." << std::endl
			<< "    --opencl-platform <n>  When mining using -G/--opencl use OpenCL platform n (default: 0)." << std::endl
			<< "    --opencl-device <n>  When mining using -G/--opencl use OpenCL device n (default: 0)." << std::endl
			<< "    -t, --mining-threads <n> Limit number of CPU/GPU miners to n (default: use everything available on selected platform)" << std::endl
			<< "    --allow-opencl-cpu Allows CPU to be considered as an OpenCL device if the OpenCL platform supports it." << std::endl
			<< "    --list-devices List the detected OpenCL devices and exit." << std::endl
			<< "    --current-block Let the miner know the current block number at configuration time. Will help determine DAG size and required GPU memory." << std::endl
			<< "    --disable-submit-hashrate  When mining, don't submit hashrate to node." << std::endl
#if ETH_ETHASHCL
			<< "    --cl-extragpu-mem Set the memory (in MB) you believe your GPU requires for stuff other than mining. Windows rendering e.t.c.." << std::endl
			<< "    --cl-local-work Set the OpenCL local work size. Default is " << dev::toString(ethash_cl_miner::c_defaultLocalWorkSize) << std::endl
			<< "    --cl-global-work Set the OpenCL global work size as a multiple of the local work size. Default is " << dev::toString(ethash_cl_miner::c_defaultGlobalWorkSizeMultiplier) << " * " << dev::toString(ethash_cl_miner::c_defaultLocalWorkSize) << std::endl
			<< "    --cl-ms-per-batch Set the OpenCL target milliseconds per batch (global workgroup size). Default is " << dev::toString(ethash_cl_miner::c_defaultMSPerBatch) << ". If 0 is given then no autoadjustment of global work size will happen" << std::endl
#endif // ETH_ETHASHCL
			;
	}

	enum class MinerType
	{
		CPU,
		GPU
	};

	std::string minerType() const { return m_minerType; }
	bool shouldPrecompute() const { return m_precompute; }

private:
	void doInitDAG(unsigned _n)
	{
		dev::h256 seedHash = dev::eth::EthashAux::seedHash(_n);
		std::cout << "Initializing DAG for epoch beginning #" << (_n / 30000 * 30000) << " (seedhash " << seedHash.abridged() << "). This will take a while." << std::endl;
		dev::eth::EthashAux::full(seedHash, true);
		exit(0);
	}

	void doBenchmark(std::string _m, bool _phoneHome, unsigned _warmupDuration = 15, unsigned _trialDuration = 3, unsigned _trials = 5)
	{
		dev::eth::BlockHeader genesis;
		genesis.setDifficulty(1 << 18);
		cdebug << dev::eth::Ethash::boundary(genesis);

		dev::eth::GenericFarm<dev::eth::EthashProofOfWork> f;
		std::map<std::string, dev::eth::GenericFarm<dev::eth::EthashProofOfWork>::SealerDescriptor> sealers;
		sealers["cpu"] = dev::eth::GenericFarm<dev::eth::EthashProofOfWork>::SealerDescriptor{&dev::eth::EthashCPUMiner::instances, [](dev::eth::GenericMiner<dev::eth::EthashProofOfWork>::ConstructionInfo ci){ return new dev::eth::EthashCPUMiner(ci); }};
#if ETH_ETHASHCL
		sealers["opencl"] = dev::eth::GenericFarm<dev::eth::EthashProofOfWork>::SealerDescriptor{&dev::eth::EthashGPUMiner::instances, [](dev::eth::GenericMiner<dev::eth::EthashProofOfWork>::ConstructionInfo ci){ return new dev::eth::EthashGPUMiner(ci); }};
#endif // ETH_ETHASHCL
		f.setSealers(sealers);
		f.onSolutionFound([&](dev::eth::EthashProofOfWork::Solution) { return false; });

		std::string platformInfo =
			_m == "cpu" ? dev::eth::EthashCPUMiner::platformInfo() :
#if ETH_ETHASHCL
			_m == "opencl" ? dev::eth::EthashGPUMiner::platformInfo() :
#endif // ETH_ETHASHCL
			"";
		std::cout << "Benchmarking on platform: " << platformInfo << std::endl;

		std::cout << "Preparing DAG..." << std::endl;
		dev::eth::Ethash::ensurePrecomputed(0);

		genesis.setDifficulty(dev::u256(1) << 63);
		f.setWork(genesis);
		f.start(_m);

		std::map<dev::u256, dev::eth::WorkingProgress> results;
		dev::u256 mean = 0;
		dev::u256 innerMean = 0;
		for (unsigned i = 0; i <= _trials; ++i)
		{
			if (!i)
				std::cout << "Warming up..." << std::endl;
			else
				std::cout << "Trial " << i << "... " << std::flush;
			boost::this_thread::sleep_for(boost::chrono::seconds(i ? _trialDuration : _warmupDuration));

			auto mp = f.miningProgress();
			f.resetMiningProgress();
			if (!i)
				continue;
			auto rate = mp.rate();

			std::cout << rate << std::endl;
			results[rate] = mp;
			mean += rate;
		}
		f.stop();
		int j = -1;
		for (auto const& r: results)
			if (++j > 0 && j < (int)_trials - 1)
				innerMean += r.second.rate();
		innerMean /= (_trials - 2);
		std::cout << "min/mean/max: " << results.begin()->second.rate() << "/" << (mean / _trials) << "/" << results.rbegin()->second.rate() << " H/s" << std::endl;
		std::cout << "inner mean: " << innerMean << " H/s" << std::endl;

		(void)_phoneHome;
#if ETH_JSONRPC
		if (_phoneHome)
		{
			std::cout << "Phoning home to find world ranking..." << std::endl;
			jsonrpc::HttpClient client("http://gav.ethdev.com:3000");
			PhoneHome rpc(client);
			try
			{
				unsigned ranking = rpc.report_benchmark(platformInfo, (int)innerMean);
				std::cout << "Ranked: " << ranking << " of all benchmarks." << std::endl;
			}
			catch (...)
			{
			}
		}
#endif // ETH_JSONRPC
		exit(0);
	}

	// dummy struct for special exception.
	struct NoWork {};
	void doFarm(std::string _m, std::string const& _remote, unsigned _recheckPeriod)
	{
		std::map<std::string, dev::eth::GenericFarm<dev::eth::EthashProofOfWork>::SealerDescriptor> sealers;
		sealers["cpu"] = dev::eth::GenericFarm<dev::eth::EthashProofOfWork>::SealerDescriptor{&dev::eth::EthashCPUMiner::instances, [](dev::eth::GenericMiner<dev::eth::EthashProofOfWork>::ConstructionInfo ci){ return new dev::eth::EthashCPUMiner(ci); }};
#if ETH_ETHASHCL
		sealers["opencl"] = dev::eth::GenericFarm<dev::eth::EthashProofOfWork>::SealerDescriptor{&dev::eth::EthashGPUMiner::instances, [](dev::eth::GenericMiner<dev::eth::EthashProofOfWork>::ConstructionInfo ci){ return new dev::eth::EthashGPUMiner(ci); }};
#endif // ETH_ETHASHCL
		(void)_m;
		(void)_remote;
		(void)_recheckPeriod;
#if ETH_JSONRPC
		jsonrpc::HttpClient client(_remote);

		dev::h256 id = dev::h256::random();
		::FarmClient rpc(client);
		dev::eth::GenericFarm<dev::eth::EthashProofOfWork> f;
		f.setSealers(sealers);
		f.start(_m);

		dev::eth::EthashProofOfWork::WorkPackage current;
		dev::eth::EthashAux::FullType dag;
		while (true)
			try
			{
				bool completed = false;
				dev::eth::EthashProofOfWork::Solution solution;
				f.onSolutionFound([&](dev::eth::EthashProofOfWork::Solution sol)
				{
					solution = sol;
					completed = true;
					return true;
				});
				
				while (!completed)
				{
					auto mp = f.miningProgress();
					f.resetMiningProgress();
					if (current)
						minelog << "Mining on PoWhash" << current.headerHash << ": " << mp;
					else
						minelog << "Getting work package...";

					if (m_submitHashrate)
					{
						auto rate = mp.rate();
						try
						{
							rpc.eth_submitHashrate(dev::toJS((dev::u256)rate), "0x" + id.hex());
						}
						catch (jsonrpc::JsonRpcException const& _e)
						{
							cwarn << "Failed to submit hashrate.";
							cwarn << boost::diagnostic_information(_e);
						}
					}

					Json::Value v = rpc.eth_getWork();
					if (v[0].asString().empty())
						throw NoWork();
					dev::h256 hh(v[0].asString());
					dev::h256 newSeedHash(v[1].asString());
					if (current.seedHash != newSeedHash)
						minelog << "Grabbing DAG for" << newSeedHash;
					if (!(dag = dev::eth::EthashAux::full(newSeedHash, true, [&](unsigned _pc){ std::cout << "\rCreating DAG. " << _pc << "% done..." << std::flush; return 0; })))
						BOOST_THROW_EXCEPTION(dev::eth::DAGCreationFailure());
					if (m_precompute)
						dev::eth::EthashAux::computeFull(sha3(newSeedHash), true);
					if (hh != current.headerHash)
					{
						current.headerHash = hh;
						current.seedHash = newSeedHash;
						current.boundary = dev::h256(dev::fromHex(v[2].asString()), dev::h256::AlignRight);
						minelog << "Got work package:";
						minelog << "  Header-hash:" << current.headerHash.hex();
						minelog << "  Seedhash:" << current.seedHash.hex();
						minelog << "  Target: " << dev::h256(current.boundary).hex();
						f.setWork(current);
					}
					boost::this_thread::sleep_for(boost::chrono::milliseconds(_recheckPeriod));
				}
				cnote << "Solution found; Submitting to" << _remote << "...";
				cnote << "  Nonce:" << solution.nonce.hex();
				cnote << "  Mixhash:" << solution.mixHash.hex();
				cnote << "  Header-hash:" << current.headerHash.hex();
				cnote << "  Seedhash:" << current.seedHash.hex();
				cnote << "  Target: " << dev::h256(current.boundary).hex();
				cnote << "  Ethash: " << dev::h256(dev::eth::EthashAux::eval(current.seedHash, current.headerHash, solution.nonce).value).hex();
				if (dev::eth::EthashAux::eval(current.seedHash, current.headerHash, solution.nonce).value < current.boundary)
				{
					bool ok = rpc.eth_submitWork("0x" + toString(solution.nonce), "0x" + toString(current.headerHash), "0x" + toString(solution.mixHash));
					if (ok)
						cnote << "B-) Submitted and accepted.";
					else
						cwarn << ":-( Not accepted.";
				}
				else
					cwarn << "FAILURE: GPU gave incorrect result!";
				current.reset();
			}
			catch (jsonrpc::JsonRpcException&)
			{
				for (auto i = 3; --i; boost::this_thread::sleep_for(boost::chrono::seconds(1)))
					std::cerr << "JSON-RPC problem. Probably couldn't connect. Retrying in " << i << "... \r";
				std::cerr << std::endl;
			}
			catch (NoWork&)
			{
				boost::this_thread::sleep_for(boost::chrono::milliseconds(100));
			}

#endif // ETH_JSONRPC
		exit(0);
	}

	/// Operating mode.
	OperationMode mode;

	/// Mining options
	std::string m_minerType = "cpu";
	unsigned m_openclPlatform = 0;
	unsigned m_openclDevice = 0;
	unsigned m_miningThreads = UINT_MAX;
	bool m_shouldListDevices = false;
	bool m_clAllowCPU = false;
#if ETH_ETHASHCL
	unsigned m_globalWorkSizeMultiplier = ethash_cl_miner::c_defaultGlobalWorkSizeMultiplier;
	unsigned m_localWorkSize = ethash_cl_miner::c_defaultLocalWorkSize;
	unsigned m_msPerBatch = ethash_cl_miner::c_defaultMSPerBatch;
#endif // ETH_ETHASHCL
	uint64_t m_currentBlock = 0;
	// default value is 350MB of GPU memory for other stuff (windows system rendering, e.t.c.)
	unsigned m_extraGPUMemory = 350000000;

	/// DAG initialisation param.
	unsigned m_initDAG = 0;

	/// Benchmarking params
	bool m_phoneHome = true;
	unsigned m_benchmarkWarmup = 3;
	unsigned m_benchmarkTrial = 3;
	unsigned m_benchmarkTrials = 5;

	/// Farm params
	std::string m_farmURL = "http://127.0.0.1:8545";
	unsigned m_farmRecheckPeriod = 500;
	bool m_precompute = true;
	bool m_submitHashrate = true;
};
