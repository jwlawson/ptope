#include <benchmark/benchmark.h>

#include <polytope_check.h>

#include "calc.h"

using ptope::calc::min_cos_angle;
static void EsselmanPolytopeCheck(benchmark::State& state) {
	ptope::PolytopeCandidate p({ { 1, -.5, 0, 0 }, 
												{ -.5, 1, min_cos_angle(4), 0 }, 
												{ 0, min_cos_angle(4), 1, -.5 }, 
												{ 0, 0, -.5, 1 } });
	auto q = p.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	auto r = q.extend_by_inner_products({ min_cos_angle(8), 0, 0, 0 });

	ptope::PolytopeCheck pcheck;
  while (state.KeepRunning()) {
		pcheck(p);
	}
}
BENCHMARK(EsselmanPolytopeCheck);

static void TumarkinPolytopeCheck(benchmark::State& state) {
	ptope::PolytopeCandidate p( { { 1, min_cos_angle(4), 0, 0 },
			{ min_cos_angle(4), 1, -.5, 0 },
			{0, -.5, 1, -.5 },
			{ 0, 0, -.5, 1 } });
	auto q = p.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });
	auto r = q.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	r.rebase_vectors({ 1, 2, 3, 4 });
	auto s = r.extend_by_inner_products({ 0, 0, min_cos_angle(4), 0 });

	ptope::PolytopeCheck pcheck;
  while (state.KeepRunning()) {
		pcheck(s);
	}
}
BENCHMARK(TumarkinPolytopeCheck);

static void NotPolytopeCheck(benchmark::State& state) {
	ptope::PolytopeCandidate p( { { 1, min_cos_angle(4), 0, 0 },
			{ min_cos_angle(4), 1, -.5, 0 },
			{0, -.5, 1, -.5 },
			{ 0, 0, -.5, 1 } });
	auto q = p.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });
	auto r = q.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });

	ptope::PolytopeCheck pcheck;
  while (state.KeepRunning()) {
		pcheck(r);
	}
}
BENCHMARK(NotPolytopeCheck);

static void NotBigPolytopeCheck(benchmark::State& state) {
	ptope::PolytopeCandidate p( { { 1, min_cos_angle(4), 0, 0 },
			{ min_cos_angle(4), 1, -.5, 0 },
			{0, -.5, 1, -.5 },
			{ 0, 0, -.5, 1 } });
	auto q = p.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });
	auto r = q.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	r.rebase_vectors({ 1, 2, 3, 4 });
	auto s = r.extend_by_inner_products({ min_cos_angle(4), 0, min_cos_angle(3), min_cos_angle(8) }); 

	ptope::PolytopeCheck pcheck;
  while (state.KeepRunning()) {
		pcheck(s);
	}
}
BENCHMARK(NotBigPolytopeCheck);

BENCHMARK_MAIN();
