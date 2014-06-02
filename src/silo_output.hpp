/*
 * silo_output.hpp
 *
 *  Created on: May 31, 2014
 *      Author: dmarce1
 */

#ifndef SILO_OUTPUT_HPP_
#define SILO_OUTPUT_HPP_

namespace xtree {

template<int Ndim, int NFields>
class silo_output: public hpx::components::managed_component_base<silo_output<Ndim, NFields>> {
public:
	static constexpr double precision = 1.0e-10;
	static constexpr int Nchild = 1 << Ndim;
	struct zone {
		std::vector<double> fields;
		std::vector<double> position;
		double span;
		zone() :
				fields(NFields), position(Ndim) {
		}
		zone(const zone& z) :
				fields(NFields), position(Ndim) {
			fields = z.fields;
			position = z.position;
			span = z.span;
		}
		zone(zone&& z) :
				fields(NFields), position(Ndim) {
			fields = std::move(z.fields);
			position = std::move(z.position);
			span = std::move(z.span);
		}
	};
	struct silo_zone {
		std::vector<double> fields;
		std::vector<int> vertices;
		silo_zone() :
				fields(NFields), vertices(Nchild) {
		}
		silo_zone(const silo_zone& s) :
				fields(NFields), vertices(Nchild) {
			fields = s.fields;
			vertices = s.vertices;
		}
		silo_zone(silo_zone&& s) :
				fields(NFields), vertices(Nchild) {
			fields = std::move(s.fields);
			vertices = std::move(s.vertices);
		}
	};
	struct vertex: public std::vector<double> {
		int index;
		vertex() :
				std::vector<double>(Ndim) {
		}
		vertex(const vertex& v) :
				std::vector<double>(Ndim, v) {
			index = v.index;
		}
		vertex(vertex&& v) :
				std::vector<double>(Ndim) {
			std::vector<double>::operator=(std::move(*(static_cast<std::vector<double>*>(&v))));
			index = v.index;
		}
	};
	struct vertex_less_functor: std::binary_function<vertex, vertex, bool> {
		bool operator()(const vertex& x, const vertex& y) const {
			for (int i = 0; i < Ndim; i++) {
				if (y[i] - x[i] < precision * 0.5) {
					return false;
				}
			}
			return true;
		}
	};
	using vertex_dir_type = std::set<vertex,vertex_less_functor>;
	using silo_zone_dir_type = std::list<silo_zone>;
private:
	std::vector<std::string> names;
	int current_index;
	vertex_dir_type nodedir;
	silo_zone_dir_type zonedir;
	std::vector<bool> received;
	mutable hpx::lcos::local::mutex mutex0;
public:
	silo_output() :
			names(NFields) {
		received.resize((hpx::find_all_localities()).size());
		reset();

	}
	virtual ~silo_output() {
	}
	bool all_received() const {
		mutex0.lock();
		for (int i = 0; i < received.size(); i++) {
			if (!received[i]) {
				return false;
			}
		}
		mutex0.unlock();
		return true;
	}

	void do_output() {
		mutex0.lock();
		constexpr int nshapes = 1;
		const int nnodes = nodedir.size();
		const int nzones = zonedir.size();
		int shapesize[1] = { Nchild };
		int shapetype[1] = { DB_ZONETYPE_HEX };
		int shapecnt[1] = { nzones };
		DBfile* db;
		DBoptlist* olist;
		double* coords[Ndim];
		char* coordnames[Ndim];

		for (int di = 0; di < Ndim; di++) {
			coords[di] = new double[nnodes];
			coordnames[di] = new char[2];
			coordnames[di][0] = 'x' + char(di);
			coordnames[di][1] = '\0';
		}

		for (auto ni = nodedir.begin(); ni != nodedir.end(); ni++) {
			for (int di = 0; di < Ndim; di++) {
				coords[di][ni->index] = (*ni)[di];
			}
		}
		nodedir.clear();

		int* zone_nodes = new int[zonedir.size() * Nchild];
		int zni = 0;
		for (auto zi = zonedir.begin(); zi != zonedir.end(); zi++) {
			for (int ci0 = 0; ci0 < Nchild; ci0++) {
				zone_nodes[zni] = zi->vertices[ci0];
				zni++;
			}
		}

		olist = DBMakeOptlist(1);
		db = DBCreate("X.silo", DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
		DBPutZonelist2(db, "zones", nzones, Ndim, zone_nodes, Nchild * nzones, 0, 0, 0, shapetype, shapesize, shapecnt, nshapes, olist);
		DBPutUcdmesh(db, "mesh", Ndim, coordnames, reinterpret_cast<float **>(coords), nnodes, nzones, "zones", NULL, DB_DOUBLE, olist);

		double* data = new double[nzones];
		char fname[2];
		fname[1] = '\0';
		for (int fi = 0; fi < NFields; fi++) {
			int i = 0;
			for (auto zi = zonedir.begin(); zi != zonedir.end(); zi++) {
				data[i] = zi->fields[fi];
				i++;
			}
			fname[0] = '1' + char(fi);
			DBPutUcdvar1(db, fname, "mesh", data, nzones, 0, 0, DB_DOUBLE, DB_ZONECENT, olist);
		}

		zonedir.clear();

		delete[] zone_nodes;
		for (int di = 0; di < Ndim; di++) {
			delete[] coords[di];
			delete[] coordnames[di];
		}
		DBClose(db);
		mutex0.unlock();
		reset();
	}

	void reset() {
		mutex0.lock();
		current_index = 0;
		for (int i = 0; i < received.size(); i++) {
			received[i] = false;
		}
		mutex0.unlock();
	}

	void send_zones_to_silo(int proc_num_from, std::vector<zone> zones) {
		constexpr int vertex_order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
		const int sz = zones.size();
		assert(!received[proc_num_from]);
		for (int i = 0; i < sz; i++) {
			silo_zone s;
			int j;
			s.fields = std::move(zones[i].fields);
			for (int ci0; ci0 < Nchild; ci0++) {
				vertex v;
				child_index_type<Ndim> ci;
				ci = vertex_order[ci0];
				for (int k = 0; k < Nchild; k++) {
					v[k] = zones[i].position[k] + double(2 * ci[k] - 1) * zones[i].span * 0.5;
				}
				mutex0.lock();
				auto iter = nodedir.find(v);
				if (iter == nodedir.end()) {
					j = current_index;
					v.index = j;
					current_index++;
					nodedir.insert(std::move(v));
				} else {
					j = iter->index;
				}
				mutex0.unlock();
				s.vertices[ci] = j;
			}
			mutex0.lock();
			zonedir.push_back(std::move(s));
			mutex0.unlock();
		}
		mutex0.lock();
		received[proc_num_from] = true;
		mutex0.unlock();
		if (all_received()) {
			hpx::apply([this]() {do_output();});
		}

	}
	XTREE_MAKE_ACTION( action_send_zones_to_siloe, silo_output::send_zones_to_silo );
}
;

} /* namespace xtree */

#endif /* SILO_OUTPUT_HPP_ */
