/*
 * silo_output.hpp
 *
 *  Created on: May 31, 2014
 *      Author: dmarce1
 */

#ifndef SILO_OUTPUT_HPP_
#define SILO_OUTPUT_HPP_
#ifndef KILL_SILO_DEP

namespace xtree {

template<int Ndim>
class silo_output: public hpx::components::managed_component_base<silo_output<Ndim>> {
public:
	static const double precision = 1.0e-10;
	static const int Nchild = 1 << Ndim;
	struct zone {
		std::vector<double> fields;
		std::array<double, Ndim> position;
		std::array<double, Ndim> span;
		zone() {
		}
		zone(const zone& z) {
			*this = z;
		}
		zone(zone&& z) {
			*this = z;
		}
		zone& operator=(const zone& z) {
			fields = z.fields;
			position = z.position;
			span = z.span;
			return *this;
		}
		zone& operator=(zone&& z) {
			fields = std::move(z.fields);
			position = std::move(z.position);
			span = std::move(z.span);
			return *this;
		}
		template<typename Archive>
		void serialize(Archive& ar, const int v) {
			ar & fields;
			ar & position;
			ar & span;
		}
	};
	struct silo_zone {
		std::vector<double> fields;
		std::vector<int> vertices;
		silo_zone() :
				vertices(Nchild) {
		}
		silo_zone(const silo_zone& s) :
				vertices(Nchild) {
			fields = s.fields;
			vertices = s.vertices;
		}
		silo_zone(silo_zone&& s) :
				vertices(Nchild) {
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
				if (x[i] - y[i] > precision) {
					return false;
				} else if (x[i] - y[i] < -precision) {
					return true;
				}
			}
			return false;
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
	mutable hpx::lcos::local::spinlock mutex0;
public:
	silo_output() {
	//	printf("Silo in\n");
		received.resize((hpx::find_all_localities()).size());
		reset();
	//	printf("Silo out\n");

	}
	virtual ~silo_output() {
	}

	void do_output() {
		mutex0.lock();
		const int nshapes = 1;
		const int nnodes = nodedir.size();
		const int nzones = zonedir.size();
		int shapesize[1] = { Nchild };
		int shapetype[1] = { DB_ZONETYPE_HEX };
		int shapecnt[1] = { nzones };
		DBfile* db;
		DBoptlist* olist;
		std::vector<double> coord_vectors[Ndim];
		std::string coordname_strs[Ndim];
		double* coords[Ndim];
		const char* coordnames[Ndim];
		const int nfields = (zonedir.begin())->fields.size();

		for (int di = 0; di < Ndim; di++) {
			coord_vectors[di].resize(nnodes);
			coords[di] = coord_vectors[di].data();
			coordname_strs[di] = ('x' + char(di));
			coordnames[di] = coordname_strs[di].c_str();
		}

		for (auto ni = nodedir.begin(); ni != nodedir.end(); ni++) {
			for (int di = 0; di < Ndim; di++) {
				coords[di][ni->index] = (*ni)[di];
				//	printf("%e ", (*ni)[di]);
			}
			//	printf("\n");
		}
		nodedir.clear();

		std::vector<int> zone_nodes(zonedir.size() * Nchild);
		int zni = 0;
		for (auto zi = zonedir.begin(); zi != zonedir.end(); zi++) {
			for (int ci0 = 0; ci0 < Nchild; ci0++) {
				zone_nodes[zni] = zi->vertices[ci0];
				//		printf("%i ", zone_nodes[zni]);
				zni++;
			}
			//		printf("\n");
		}

		olist = DBMakeOptlist(1);
		db = DBCreate("X.silo", DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
		DBPutZonelist2(db, "zones", nzones, Ndim, zone_nodes.data(), Nchild * nzones, 0, 0, 0, shapetype, shapesize,
				shapecnt, nshapes, olist);
		DBPutUcdmesh(db, "mesh", Ndim, const_cast<char**>(coordnames), reinterpret_cast<float **>(coords), nnodes,
				nzones, "zones", NULL, DB_DOUBLE, olist);

		std::vector<double> data(nzones);
		char fname[2];
		fname[1] = '\0';
		for (int fi = 0; fi != nfields; ++fi) {
			int i = 0;
			for (auto zi = zonedir.begin(); zi != zonedir.end(); zi++) {
				data[i] = zi->fields[fi];
				i++;
			}
			fname[0] = 'A' + char(fi);
			DBPutUcdvar1(db, fname, "mesh", reinterpret_cast<float*>(data.data()), nzones, 0, 0, DB_DOUBLE, DB_ZONECENT, olist);
		}

		zonedir.clear();

		DBClose(db);
		mutex0.unlock();
		reset();
	}

	void reset() {
		mutex0.lock();
		current_index = 0;
		std::fill(received.begin(), received.end(), false);
		mutex0.unlock();
	}

	void send_zones_to_silo(int proc_num_from, std::vector<zone> zones) {
		const int vertex_order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
		const int sz = zones.size();
		assert(!received[proc_num_from]);
		for (int i = 0; i < sz; i++) {
			silo_zone s;
			int j;
			s.fields = std::move(zones[i].fields);
			//	printf("%e %e %e %e %e %e\n", zones[i].position[0], zones[i].position[1], zones[i].position[2], zones[i].span[0], zones[i].span[1],
			//			zones[i].span[2]);
			for (int ci0 = 0; ci0 < Nchild; ci0++) {
				vertex v;
				int ci = vertex_order[ci0];
				for (int k = 0; k < Ndim; k++) {
					const double factor = (0.5 * double(2 * ((ci >> k) & 1) - 1));
					v[k] = zones[i].position[k] + zones[i].span[k] * factor;
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
				s.vertices[ci0] = j;
			}
			mutex0.lock();
			zonedir.push_back(std::move(s));
			mutex0.unlock();
		}
		mutex0.lock();
		received[proc_num_from] = true;
		printf( "Receive output from %i\n", proc_num_from );
		if (std::all_of(received.begin(), received.end(), [](bool b) {return b;})) {
			printf( "Doing output\n");
			mutex0.unlock();
			/*hpx::apply([this]() {*/do_output();/*});*/
			mutex0.lock();
		}
		mutex0.unlock();

	}
	HPX_DEFINE_COMPONENT_ACTION_TPL( silo_output,send_zones_to_silo,action_send_zones_to_silo );
}
;

} /* namespace xtree */

#endif

#endif /* SILO_OUTPUT_HPP_ */
